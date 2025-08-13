from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from numpy import where # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.dds import DeseqDataSet # type: ignore
import subprocess
from typing import Union
from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class ContrastFoldChangeManager:
    def run(self, count_contrasts_path: Path) -> dict[set[str]]:
        count_contrasts = self.extract_count_contrasts(count_contrasts_path)
        return self.extract_significant_genes(count_contrasts)

    @staticmethod
    def extract_count_contrasts(count_contrasts: Path) -> pd.DataFrame:
        df = pd.read_csv(count_contrasts, index_col=0)
        columns = [eval(column) for column in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        return df

    @staticmethod
    def extract_significant_genes(count_contrasts: pd.DataFrame) -> dict[set[str]]:
        contrasts_significant_genes = dict()
        for contrast, deseq_results in count_contrasts.T.groupby(level=0):
            deseq_results = deseq_results.T.droplevel(0, axis=1)
            significant_genes = set()
            for time, deseq_results_by_time in deseq_results.T.groupby(level=0):
                deseq_results_by_time = deseq_results_by_time.T.droplevel(0, axis=1)
                deseq_results_by_time["Significant"] = where((abs(deseq_results_by_time["log2FoldChange"]) > 1.5) & (deseq_results_by_time["padj"] < 0.05), "yes", "no")
                significant_genes.update(set(deseq_results_by_time.loc[deseq_results_by_time["Significant"] == "yes"].index))
            contrasts_significant_genes[contrast] = significant_genes
        return contrasts_significant_genes

class SampleMetadataManager:
    def __init__(self, sample_metadata_path: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata_path)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> pd.DataFrame:
        data = pd.read_csv(sample_metadata_path, index_col=0, dtype={"blacklist": "Int32"})
        data = data[data["blacklist"].isnull()]
        return data[data["technical_replicate_parent_uid"].isnull()]

class DegPatternManager:
    def __init__(self, outdir: Path) -> None:
        self.outdir = outdir

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, contrasts_significant_genes: dict[set[str]], cell_line: str,
            cell_line_col: str, treatment_col: str, covariate_cols: list[str]) -> None:
        design_factors = covariate_cols
        print(f"\nStarting on cell line: {cell_line}\n----------")
        gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line, cell_line_col)
        
        for contrast, significant_genes in contrasts_significant_genes.items():
            print(f"\nStarting on contrast: {contrast}")
            treatment, control = self.extract_conditions_from_contrast(contrast)
            gene_counts_by_virus, sample_metadata_by_virus = \
                self.filter_for_virus_contrasts(gene_counts_by_cell_line,sample_metadata_by_cell_line, treatment, control, treatment_col)
            
            dds = self.pydeseq2_normalization(gene_counts_by_virus, sample_metadata_by_virus, design_factors)
            
            self.perform_degpatterns_clustering(significant_genes, dds, sample_metadata_by_virus, self.outdir,
                                                cell_line, treatment, control, treatment_col)

    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                             cell_line: str, cell_line_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata[cell_line_col] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata
    
    @staticmethod
    def extract_conditions_from_contrast(contrast: str) -> tuple[str, str]:
        conditions = contrast.split("_vs_")
        return conditions[0], conditions[1]

    @staticmethod
    def filter_for_virus_contrasts(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                                   treatment: str, control: str, treatment_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        virus_list = [treatment, control]
        sample_metadata = sample_metadata[sample_metadata[treatment_col].isin(virus_list)]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata

    @staticmethod
    def pydeseq2_normalization(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, design_factors: list[str]) -> DeseqDataSet:
        print("Starting pyDEseq2 analysis")
        with catch_warnings():
            simplefilter("ignore")
            dds = DeseqDataSet(counts=gene_counts,
                            metadata=sample_metadata,
                            design_factors=design_factors)
        with suppress_stdout_stderr(): # NOTE: Primarily suppresses redundant messages in prod, but could suppress actual errors important for dev work
            dds.deseq2()
        return dds

    @classmethod
    def perform_degpatterns_clustering(cls, significant_genes: set[str], dds: DeseqDataSet, sample_metadata: pd.DataFrame,
                                       write_directory: Path, cell_line: str, treatment: str, control: str, treatment_col: str) -> None:
        normalized_counts = cls.extract_normalized_count_df_from_dds(dds)
        normalized_counts_of_interest = normalized_counts[list(significant_genes)].copy()
        # NOTE: A stupid way to get the control variable at the end of a sort so the labels are consistently colored by degpatterns:
        sample_metadata[treatment_col] = sample_metadata[treatment_col].apply(lambda x: f"~{control}" if x == control else x)

        gene_count_path, sample_metadata_path = cls.write_input_data_for_degpattern_clustering(normalized_counts_of_interest, sample_metadata, write_directory)
        degpattern_results_path = write_directory / f"{cell_line}_{treatment}_vs_{control}_gene_clusters.csv"
        degpattern_figure_path = write_directory / f"{cell_line}_{treatment}_vs_{control}_gene_clusters.pdf"
        cls.run_r_degpatterns_command(gene_count_path, sample_metadata_path, degpattern_results_path, degpattern_figure_path)

    @staticmethod
    def extract_normalized_count_df_from_dds(dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts = dds.layers["normed_counts"]
        gene_ids = dds._var.index.to_list()
        sample_ids = dds.obsm["design_matrix"].index.to_list()
        return pd.DataFrame(normalized_counts, index=sample_ids, columns=gene_ids)

    @staticmethod
    def write_input_data_for_degpattern_clustering(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, directory: Path) -> tuple[Path, Path]:
        gene_count_outpath = directory / "gene_counts_degpattern_input.csv"
        sample_metadata_outpath = directory / "sample_metadata_degpattern_input.csv"
        gene_counts.to_csv(gene_count_outpath, index=True)
        sample_metadata.to_csv(sample_metadata_outpath, index=True)
        return gene_count_outpath, sample_metadata_outpath

    @staticmethod
    def run_r_degpatterns_command(gene_count_path: Path, sample_metadata_path: Path, degpattern_results_path: Path, degpattern_figure_path: Path) -> None:
        print("Starting degPatterns in R")
        r_degpatterns_command = ["Rscript", "/src/app/deg_patterns/degpattern_clustering.R",
                         gene_count_path, sample_metadata_path, degpattern_results_path, degpattern_figure_path]

        p = subprocess.Popen(r_degpatterns_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline().strip()) != "":
            pass
        p.wait()

        if p.poll() != 0:
            print(f"Exit code: {p.poll()}")
            for line in p.stderr.readlines():
                print(line.strip())
            raise Exception("R degPatterns did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_counts", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-gene_contrast_path", type=str, required=True)
    parser.add_argument("-cell_line", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-cell_line_col", type=str, default="cell_line")
    parser.add_argument("-treatment_col", type=str, default="treatment")
    parser.add_argument("-covariate_cols", nargs="*", default=["treatment", "time_point"])

    args = parser.parse_args()
    
    cfcm = ContrastFoldChangeManager()
    contrasts_significant_genes = cfcm.run(Path(args.gene_contrast_path))
    
    gene_counts = pd.read_csv(args.gene_counts, index_col=0).sort_index()
    smm = SampleMetadataManager(Path(args.sample_metadata))

    dpm = DegPatternManager(Path(args.outdir))
    dpm.run(gene_counts, smm.sample_metadata, contrasts_significant_genes, args.cell_line,
            args.cell_line_col, args.treatment_col, args.covariate_cols)