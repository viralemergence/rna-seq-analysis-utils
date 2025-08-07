from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from gseapy import enrichr # type: ignore
import matplotlib.pyplot as plt # type: ignore
from numpy import log10, where # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
import seaborn as sns # type: ignore
from typing import Any

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class BackgroundGeneManager:
    @staticmethod
    def extract_background_genes_from_contrasts(count_contrasts: Path) -> set[str]:
        df = pd.read_csv(count_contrasts, index_col=0)
        columns = [eval(column) for column in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        return set(df.index)

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

    def run_coarse(self, contrast_of_interest: str, expression_change_direction: str) -> set[str]:
        return self.extract_differentially_expressed_genes(self.gene_count_contrasts, contrast_of_interest, expression_change_direction)

    def run_fine(self, contrast_of_interest: str, expression_change_direction: str) -> set[str]:
        return self.extract_differentially_expressed_genes_by_timepoint(
            self.gene_count_contrasts, contrast_of_interest, expression_change_direction
            )

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

            print(f"\nExtracting coarse DEGs for contrast: {contrast}")
            deseq_results = deseq_results.droplevel(0, axis=1)
            deseq_results.columns = deseq_results.columns.set_levels(deseq_results.columns.levels[0].astype(float), level=0)

            with suppress_stdout_stderr():
                groupby_data = deseq_results.groupby(level=0, axis=1)
            for time_point, deseq_results_by_time in groupby_data:
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

    @staticmethod
    def extract_differentially_expressed_genes_by_timepoint(
        gene_count_contrasts: pd.DataFrame, contrast_of_interest: str, expression_change_direction: str,
        log2foldchange_threshold: float = 1.5) -> set[str]:

        with suppress_stdout_stderr():
            groupby_data = gene_count_contrasts.groupby(level=0, axis=1)
        for contrast, deseq_results in groupby_data:
            if contrast != contrast_of_interest:
                continue
            
            significant_genes = defaultdict()

            print(f"\nExtracting fine DEGs for contrast: {contrast}")
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

                significant_genes[time_point] = significant_genes_by_time

            return significant_genes

def perform_gene_over_representation_analysis(
        degs: set[str],
        background_genes: set[str],
        gene_set: str = "Reactome_Pathways_2024") -> Any:
    enr_bg = enrichr(
        gene_sets = [gene_set],
        gene_list = list(degs),
        background = background_genes,
        outdir = None
    )
    return enr_bg

def generate_heatmap(heatmap_dataframe: pd.DataFrame, contrast: str, figure_outpath: Path) -> None:
    plt.rcParams["svg.fonttype"] = "none"
    sns.set_theme(font="sans-serif", font_scale=0.6, rc={"font.weight": "bold"})
    fig, ax = plt.subplots(figsize=(3, 3), dpi=1200)

    heatmap = sns.heatmap(heatmap_dataframe, ax=ax, cmap="Reds", square=True, cbar_kws={"shrink": 0.75})

    # Title
    title = contrast
    title = title.split("_vs_")
    title.insert(1, "vs")
    title = "\n".join(title)
    heatmap.set_title(title, fontweight="bold")

    # Colorbar and null values
    cbar = heatmap.collections[0].colorbar
    cbar.set_label("log10(Odds Ratio)", labelpad=10, fontweight="bold")
    heatmap.collections[0].cmap.set_bad("grey")

    # Truncating and setting labels
    heatmap.set_xlabel("HPI", fontweight="bold")
    heatmap.set_ylabel("")
    
    labels = []
    for label in heatmap.get_yticklabels():
        text = label.get_text()
        if len(text) > 25:
            text = text[:25] + "..." # cls.kegg_pathway_id_mapper(text)
        labels.append(text)
    heatmap.set_yticklabels(labels)
    plt.tight_layout()
    
    heatmap.tick_params(left=False, bottom=False)

    # Saving graph
    fig.savefig(figure_outpath, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_contrast_path", type=str, required=True)
    parser.add_argument("-contrasts", nargs="+", required=True)
    parser.add_argument("-expression_direction", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    args = parser.parse_args()

    background_genes = BackgroundGeneManager.extract_background_genes_from_contrasts(Path(args.gene_contrast_path))

    gccm = GeneCountContrastManager(Path(args.gene_contrast_path))
    for contrast in args.contrasts:
        degs = gccm.run_coarse(contrast, args.expression_direction)

        gora = perform_gene_over_representation_analysis(degs, background_genes)
        top_pathways_coarse = list(gora.results.head(10)["Term"])

        collated_results = {}
        degs_by_time_point = gccm.run_fine(contrast, args.expression_direction)
        for time_point, degs in degs_by_time_point.items():
            print(f"Beginning over representation analysis for {time_point}")
            gora = perform_gene_over_representation_analysis(degs, background_genes)
            results = gora.results.set_index("Term")

            present_pathways = [pathway for pathway in top_pathways_coarse if pathway in results.index]
            pathways_of_interest_stat = results.loc[present_pathways]["Odds Ratio"]
            collated_results[time_point] = pathways_of_interest_stat

        ora_results = pd.DataFrame(collated_results,
                                   columns=degs_by_time_point.keys()).sort_values(list(degs_by_time_point.keys())[0], ascending=False)
        ora_results = ora_results.apply(lambda x: log10(x))

        figure_outpath = Path(args.outdir) / f"{gccm.cell_line}_{contrast}_{args.expression_direction}_gora_heatmap.png"
        generate_heatmap(ora_results, contrast, figure_outpath)