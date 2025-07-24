from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.ds import DeseqStats # type: ignore
from pydeseq2.dds import DeseqDataSet # type: ignore
from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class SampleMetadataManager:
    def __init__(self, sample_metadata_path: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata_path)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> pd.DataFrame:
        data = pd.read_csv(sample_metadata_path, index_col=0, dtype={"blacklist": "Int32"})
        data = data[data["blacklist"].isnull()]
        return data[data["technical_replicate_parent_uid"].isnull()]

class CalculateContrasts:
    def __init__(self, cell_lines: list[str], virus_contrasts: list[str], outdir: Path) -> None:
        self.cell_lines = cell_lines
        self.virus_contrasts = [contrast.split("_vs_") for contrast in virus_contrasts]
        self.outdir = outdir

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
            cell_line_col: str, treatment_col: str, time_point_col: str,
            covariate_cols: list[str], time_points_of_interest: list[str]) -> None:

        for cell_line in self.cell_lines:
            print(f"\nStarting on cell line: {cell_line}")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = \
                self.filter_for_cell_line(gene_counts, sample_metadata, cell_line, cell_line_col)
            
            contrast_results = []
            contrast_labels = []
            
            for treatment, control in self.virus_contrasts:
                print(f"\nStarting on virus contrast: {treatment} vs {control}")
                gene_counts_by_virus, sample_metadata_by_virus = \
                    self.filter_for_virus_contrasts(gene_counts_by_cell_line, sample_metadata_by_cell_line,
                                                    treatment, control, treatment_col)
                
                for time_point in time_points_of_interest:
                    print(f"Starting on time: {time_point}")
                    gene_counts_by_time, sample_metadata_by_time = \
                        self.filter_for_time(gene_counts_by_virus, sample_metadata_by_virus, time_point, time_point_col)

                    dds = self.pydeseq2_normalization(gene_counts_by_time, sample_metadata_by_time, covariate_cols)
                    deseq_results = self.deseq_stats_calculation(dds, treatment, control, treatment_col)
                    self.extract_and_append_foldchange_padj(deseq_results, contrast_results, contrast_labels, treatment, control, time_point)

            collated_df = pd.concat(contrast_results, axis=1)
            collated_df.columns = contrast_labels
            collated_df = collated_df.round(6)
            print("\nSaving collated results")
            print(collated_df)
            
            out_filename = f"{cell_line}_virus_contrasts_per_time.csv"
            outpath = self.outdir / out_filename
            collated_df.to_csv(outpath)
    
    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                             cell_line: str, cell_line_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata[cell_line_col] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata

    @staticmethod
    def filter_for_virus_contrasts(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                                   treatment: str, control: str, treatment_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        virus_list = [treatment, control]
        sample_metadata = sample_metadata[sample_metadata[treatment_col].isin(virus_list)]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata

    @staticmethod
    def filter_for_time(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                        time: float, time_point_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata[time_point_col] == time]
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

    @staticmethod
    def deseq_stats_calculation(dds: DeseqDataSet, treatment: str, control: str, treatment_col: str) -> pd.DataFrame:
        results = DeseqStats(dds, contrast=[treatment_col, treatment, control.replace("_", "-")])
        with suppress_stdout_stderr():
            results.summary()
        return results.results_df

    @staticmethod
    def extract_and_append_foldchange_padj(deseq_results: pd.DataFrame, contrast_results: list[pd.DataFrame], contrast_labels: list[str],
                                           treatment: str, control: str, time_point: float) -> None:
        columns_of_interest = ["log2FoldChange", "padj"]
        for column in columns_of_interest:
            contrast_results.append(deseq_results[column])
            label = (f"{treatment}_vs_{control}", f"{time_point}", column)
            contrast_labels.append(label)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_counts", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-cell_lines", nargs="+", required=True)
    parser.add_argument("-virus_contrasts", nargs="+", required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-cell_line_col", type=str, default="cell_line")
    parser.add_argument("-treatment_col", type=str, default="treatment")
    parser.add_argument("-time_point_col", type=str, default="time_point")
    parser.add_argument("-covariate_cols", nargs="*", default=["treatment"])
    parser.add_argument("-time_points", nargs="*", default=[6.0, 12.0, 24.0, 48.0])

    args = parser.parse_args()

    gene_counts = pd.read_csv(args.gene_counts, index_col=0).sort_index()
    smm = SampleMetadataManager(Path(args.sample_metadata))
    
    cc = CalculateContrasts(args.cell_lines, args.virus_contrasts, Path(args.outdir))
    cc.run(gene_counts, smm.sample_metadata,
           args.cell_line_col, args.treatment_col, args.time_point_col,
           args.covariate_cols, args.time_points)