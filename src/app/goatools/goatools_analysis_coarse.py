from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from csv import reader, DictReader, DictWriter
from goatools.anno.genetogo_reader import Gene2GoReader # type: ignore
from goatools.base import dnld_file # type: ignore
from goatools.go_enrichment import GOEnrichmentStudy # type: ignore
from goatools.mapslim import mapslim # type: ignore
from goatools.obo_parser import GODag # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.colors import LinearSegmentedColormap # type: ignore
from numpy import where # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from typing import Any, Tuple

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class GeneCountContrastManager:
    def __init__(self, gene_count_contrasts_path: Path) -> None:
        self.gene_count_contrasts = self.extract_count_contrasts(gene_count_contrasts_path)
        self.cell_line = gene_count_contrasts_path.name.split("_")[0]

    @staticmethod
    def extract_count_contrasts(count_contrasts: Path) -> pd.DataFrame:
        df = pd.read_csv(count_contrasts, index_col=0)
        columns = [eval(column) for column in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        return df

    def run(self, contrast_of_interest: str, expression_change_direction: str) -> set[str]:
        return self.extract_differentially_expressed_genes(self.gene_count_contrasts, contrast_of_interest, expression_change_direction)

    @staticmethod
    def extract_differentially_expressed_genes(
        gene_count_contrasts: pd.DataFrame, contrast_of_interest: str, expression_change_direction: str,
        log2foldchange_threshold: float = 1.5) -> set[str]:

        with suppress_stdout_stderr():
            groupby_data = gene_count_contrasts.groupby(level=0, axis=1)
        for contrast, deseq_results in groupby_data:
            if contrast != contrast_of_interest:
                continue
            
            significant_genes = set()

            print(f"\nStarting on contrast: {contrast}")
            deseq_results = deseq_results.droplevel(0, axis=1)
            deseq_results.columns = deseq_results.columns.set_levels(deseq_results.columns.levels[0].astype(float), level=0)

            with suppress_stdout_stderr():
                groupby_data = deseq_results.groupby(level=0, axis=1)
            for time_point, deseq_results_by_time in groupby_data:
                print(f"Time point: {time_point}")
                deseq_results_by_time = deseq_results_by_time.droplevel(0, axis=1)

                if expression_change_direction == "up":
                    deseq_results_by_time["Significant"] = where(
                        (deseq_results_by_time["log2FoldChange"] > log2foldchange_threshold) &
                        (deseq_results_by_time["padj"] < 0.05), "yes", "no")
                elif expression_change_direction == "down":
                    deseq_results_by_time["Significant"] = where(
                        (deseq_results_by_time["log2FoldChange"] < (-1*log2foldchange_threshold)) &
                        (deseq_results_by_time["padj"] < 0.05), "yes", "no")
                significant_genes_by_time = set(deseq_results_by_time.loc[deseq_results_by_time["Significant"] == "yes"].index)

                significant_genes.update(significant_genes_by_time)

            return significant_genes

class GoatoolsManager:
    def __init__(self, taxon_id: int, background_genes_path: Path, goea_dir: Path, outpath: Path) -> None:
        self.taxon_id = taxon_id

        # Background genes for enrichment analysis
        self.background_genes_path = background_genes_path

        # Full and slim GO terms in obo format
        self.obo_ftp_path = "http://current.geneontology.org/ontology/go.obo"
        self.obo_path = goea_dir / "reference_files" / "go.obo"
        self.slim_obo_ftp_path = "https://current.geneontology.org/ontology/subsets/goslim_pir.obo"
        self.slim_obo_path = goea_dir / "reference_files" / "goslim-pir.obo"

        # gene2go and taxon specific gene2go info
        self.gene2go_ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
        self.gene2go_path = goea_dir / "reference_files" / "gene2go.txt"
        self.gene2go_taxon_path = goea_dir / "reference_files" / f"gene2go_{self.taxon_id}.txt"

        # Final results path
        self.outpath = outpath

    def setup(self) -> None:
        if not self.gene2go_taxon_path.is_file():
            dnld_file(self.gene2go_ftp_path, str(self.gene2go_path))
            self.filter_gene2go_by_taxon_id(self.gene2go_path, self.gene2go_taxon_path, self.taxon_id)
        dnld_file(self.obo_ftp_path, str(self.obo_path))
        dnld_file(self.slim_obo_ftp_path, str(self.slim_obo_path))
        self.get_background_genes(self.background_genes_path)
        
    @staticmethod
    def filter_gene2go_by_taxon_id(gene2go_path: Path, gene2go_taxon_path: Path, taxon_id: int) -> None:
        # Creating a taxon specific gene2go file dramatically reduces Gene2GoReader load time
        with gene2go_path.open() as inhandle, gene2go_taxon_path.open("w") as outhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            outhandle.write("\t".join(header) + "\n")
            for line in reader_iterator:
                tax_id = int(line[0])
                if tax_id != taxon_id:
                    continue
                outhandle.write("\t".join(line) + "\n")

    @staticmethod
    def get_background_genes(background_genes_path: Path) -> None:
        # NOTE: This is a WIP
        # Currenty instructions are:
        # Go to: https://www.ncbi.nlm.nih.gov/gene
        # For Rousettus search: "9407"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]
        # Send to file (must be tabular)
        pass

    def run(self, study_genes: set[str]) -> list:
        # Returns list of goatools result objects
        background_genes = self.extract_background_genes(self.background_genes_path)
        symbol_id_mapper = self.set_symbol_id_mapper(self.background_genes_path)
        study_genes = self.convert_study_gene_symbols(study_genes, symbol_id_mapper)

        godag = GODag(str(self.obo_path), optional_attrs=["relationship"])

        namespace_associations = Gene2GoReader(str(self.gene2go_taxon_path), taxids=[self.taxon_id]).get_ns2assc()
        gene_go_terms = self.extract_gene_go_terms(namespace_associations)

        go_enrichment_analysis = GOEnrichmentStudy(background_genes, gene_go_terms, godag,
                                    methods=['bonferroni'], pvalcalc='fisher_scipy_stats')
        return go_enrichment_analysis.run_study_nts(study_genes)
    
    @staticmethod
    def extract_background_genes(background_genes_path: Path) -> set[str]:
        background_genes = set()
        with background_genes_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            for line in reader_iterator:
                gene_id = int(line[2])
                background_genes.add(gene_id)
        return background_genes

    @staticmethod
    def set_symbol_id_mapper(background_genes_path: Path) -> dict[int]:
        symbol_id_mapper = {}
        with background_genes_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            for line in reader_iterator:
                gene_id = int(line[2])
                symbol = line[5]
                symbol_id_mapper[symbol] = gene_id
        return symbol_id_mapper
    
    @staticmethod
    def convert_study_gene_symbols(study_genes: set[str], symbol_id_mapper: dict[str]) -> set[str]:
        converted_study_genes = set()
        for gene in study_genes:
            try:
                converted_study_genes.add(symbol_id_mapper[gene])
            except KeyError:
                continue
        return converted_study_genes

    @staticmethod
    def extract_gene_go_terms(namespace_associations: dict[dict[set]]) -> dict[set]:
        gene_go_terms = defaultdict(set)
        for _, associations in namespace_associations.items():
            for gene_id, go_terms in associations.items():
                gene_go_terms[gene_id].update(go_terms)
        return gene_go_terms

class GoatoolsResultsManager(GoatoolsManager):
    def __init__(self, taxon_id: int, background_genes_path: Path, goea_dir: Path, outpath: Path) -> None:
        super().__init__(taxon_id, background_genes_path, goea_dir, outpath)

    def run(self, goatools_results: list) -> None:
        symbol_id_mapper = self.set_symbol_id_mapper(self.background_genes_path)
        id_symbol_mapper = {id: symbol for symbol, id in symbol_id_mapper.items()}

        godag = GODag(str(self.obo_path), optional_attrs=["relationship"])
        goslim_dag = GODag(str(self.slim_obo_path), optional_attrs=["relationship"])

        significant_results = self.extract_significant_results(goatools_results)
        grouped_results = self.group_go_terms(significant_results, id_symbol_mapper, godag, goslim_dag)
        grouped_results_highlights = self.extract_highlights_from_grouped_go_terms(grouped_results, godag)
        
        with self.outpath.open("w") as outhandle:
            header = grouped_results_highlights[0].keys()
            writer = DictWriter(outhandle, fieldnames=header)
            writer.writeheader()
            writer.writerows(grouped_results_highlights)

    @staticmethod
    def extract_significant_results(results: list) -> list:
        return [r for r in results if r.p_uncorrected < 0.05]

    @classmethod
    def group_go_terms(cls, results: list, id_symbol_mapper: dict[str], godag: GODag, goslim_dag: GODag) -> dict[list[dict]]:
        grouped_results = defaultdict(list)
        for r in results:
            go_genes = [id_symbol_mapper[gene] for gene in r.study_items]
            if len(go_genes) < 2: # NOTE: May want to modulate threshold
                continue

            group_go_term = cls.calculate_group_go_term(r.GO, godag, goslim_dag)
            
            info = {"go_id": r.GO,
                    "go_name": f'"{r.goterm.name}"',
                    "p-val": round(r.p_uncorrected, 10),
                    "genes": go_genes}
            grouped_results[group_go_term].append(info)
        return grouped_results

    @classmethod
    def calculate_group_go_term(cls, go_term: str, godag: GODag, goslim_dag: GODag) -> str:
        direct_ancestors, all_ancestors = mapslim(go_term, godag, goslim_dag)
        if len(all_ancestors) == 0:
            return "NONE"
        return cls.get_lowest_level_ancestor(godag, all_ancestors)

    @staticmethod
    def get_lowest_level_ancestor(godag: GODag, ancestors: set[str]) -> str:
        ancestor_levels = dict()
        for ancestor in ancestors:
            ancestor_levels[ancestor] = godag[ancestor].level
        return max(ancestor_levels, key= lambda x: ancestor_levels[x])

    @classmethod
    def extract_highlights_from_grouped_go_terms(cls, grouped_results: dict[list[dict]], godag: GODag) -> list[dict]:
        grouped_results_highlights = list()
        for group_go_id, results in grouped_results.items():
            group_info = {"group_go_id": group_go_id,
                          "group_go_name": f'"{godag[group_go_id].name}"'}

            most_significant_result = cls.extract_most_significant_result(group_go_id, results)
            msr = most_significant_result
            most_significant_result_info = {"most_significant_go_id": msr["go_id"],
                                            "most_significant_go_name": msr["go_name"],
                                            "most_significant_p_val": msr["p-val"],
                                            "most_significant_study_genes": "|".join(sorted(msr["genes"]))}
            group_info.update(most_significant_result_info)

            largest_result = sorted(results, key=lambda result: len(result["genes"]), reverse=True)[0]
            lr = largest_result
            largest_result_info = {"largest_go_id": lr["go_id"],
                                   "largest_go_name": lr["go_name"],
                                   "largest_p_val": lr["p-val"],
                                   "largest_study_genes": "|".join(sorted(lr["genes"]))}
            group_info.update(largest_result_info)
            
            grouped_results_highlights.append(group_info)
        return sorted(grouped_results_highlights, key=lambda result: result["most_significant_p_val"])
    
    @staticmethod
    def extract_most_significant_result(group_go_id: str, results: list[dict]) -> dict:
        significant_results = sorted(results, key=lambda result: result["p-val"])
        if (significant_results[0]["go_id"] != group_go_id) or (len(significant_results) == 1):
            return significant_results[0]
        return significant_results[1] #NOTE: May not want this, but does help remove very broad terms

class GoatoolsResultsGrapher():
    def run(self, results_path: Path, title: str) -> None:
        results = self.extract_results(results_path)
        results_top_10 = results[:10][::-1] # NOTE: Need to reverse order so they plot correctly to bar chart

        go_names = self.extract_transform_go_names(results_top_10)
        gene_counts = self.extract_gene_counts(results_top_10)
        p_vals = self.extract_pvals(results_top_10)
        
        cmap_colors = ["crimson", "purple", "blue"]
        cmap = LinearSegmentedColormap.from_list("red_to_purple", cmap_colors)
        barchart_colors = self.set_colors(cmap, p_vals)
        
        _, ax = self.make_barchart(go_names, gene_counts, barchart_colors, title)
        self.add_colorbar(cmap, p_vals, ax)
        
        figure_outpath = results_path.with_suffix(".svg")
        plt.savefig(figure_outpath, bbox_inches="tight")
        plt.close()

    @staticmethod
    def extract_results(results_path: Path) -> list[dict]:
        with results_path.open() as inhandle:
            reader_iterator = DictReader(inhandle)
            return list(reader_iterator)

    @staticmethod
    def extract_transform_go_names(results: list[dict]) -> list[str]:
        go_names = [r["most_significant_go_name"][1:-1] for r in results]
        go_ids = [r["most_significant_go_id"] for r in results]

        transformed_go_names = []
        for go_name, go_id in zip(go_names, go_ids):
            if len(go_name) > 25:
                print(f"{go_id}={go_name}")
                transformed_go_names.append(go_id)
            else:
                transformed_go_names.append(go_name)
        return transformed_go_names

    @staticmethod
    def extract_pvals(results: list[dict]) -> list[float]:
        return [float(r["most_significant_p_val"]) for r in results]
    
    @staticmethod
    def extract_gene_counts(results: list[dict]) -> list[int]:
        return [len(r["most_significant_study_genes"].split("|")) for r in results]
    
    @staticmethod
    def set_colors(cmap, p_vals: list[float]) -> list[list[float]]:
        if (max(p_vals) - min(p_vals)) == 0:
            p_vals[0] = 0.0000000000000001
        rescale = [(p_val - min(p_vals)) / (max(p_vals) - min(p_vals)) for p_val in p_vals]
        return cmap(rescale)

    @staticmethod
    def make_barchart(go_names: list[str], gene_counts: list[int], colors: list[list[float]], title: str) -> Tuple[Any]:
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["font.weight"] = "bold"
        fig, ax = plt.subplots(figsize=(4, 3), dpi=1200)
        ax.barh(go_names, gene_counts, color=colors)
        ax.tick_params(axis="y", labelsize=10)
        ax.set_title(title.replace("_", " "), fontweight="bold")
        ax.set_xlabel("Genes per GO Term", fontweight="bold")
        if max(gene_counts) > 100:
            ax.set_xscale("log")
        return fig, ax

    @staticmethod
    def add_colorbar(cmap, p_vals: list[float], ax) -> None:
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_clim(vmin=min(p_vals), vmax=max(p_vals))
        cbar = plt.colorbar(sm, ax=ax)
        cbar.ax.invert_yaxis()
        cbar.ax.set_title("p-val", fontweight="bold")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_contrast_path", type=str, required=True)
    parser.add_argument("-contrasts", nargs="+", required=True)
    parser.add_argument("-expression_direction", type=str, required=True)
    parser.add_argument("-taxon_id", type=int, required=True)
    parser.add_argument("-background_genes", type=str, required=True)
    parser.add_argument("-goea_dir", type=str, required=True)
    args = parser.parse_args()

    gccm = GeneCountContrastManager(Path(args.gene_contrast_path))
    for contrast in args.contrasts:
        degs = gccm.run(contrast, args.expression_direction)

        outpath = Path(args.goea_dir) / "coarse_results" / f"{gccm.cell_line}_{contrast}_{args.expression_direction}_goea_coarse.csv"

        gm = GoatoolsManager(args.taxon_id, Path(args.background_genes), Path(args.goea_dir), outpath)
        gm.setup()
        results = gm.run(degs)
        
        grm = GoatoolsResultsManager(args.taxon_id, Path(args.background_genes), Path(args.goea_dir), outpath)
        grm.run(results)
        
        grg = GoatoolsResultsGrapher()
        grg.run(grm.outpath, contrast)