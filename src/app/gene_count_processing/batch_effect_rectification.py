from argparse import ArgumentParser
from inmoose.pycombat import pycombat_seq # type: ignore
import pandas as pd # type: ignore
from pathlib import Path

class SampleMetadataManager:
    def __init__(self, sample_metadata_path: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata_path)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> pd.DataFrame:
        data = pd.read_csv(sample_metadata_path, index_col=0, dtype={"blacklist": "Int32"})
        data = data[data["blacklist"].isnull()]
        return data[data["technical_replicate_parent_uid"].isnull()]

class BatchEffectRectifier:
    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, outpath: Path,
            cell_line_col: str, batch_col: str, covariate_cols: list[str]) -> None:
        cell_lines = self.extract_cell_lines(sample_metadata, cell_line_col)
        print(f"Metadata listed cell lines: {cell_lines}")

        gene_counts_by_cell_line_rectified = []

        for cell_line in cell_lines:
            print(f"\nStarting on cell line: {cell_line}\n----------")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line, cell_line_col)
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.rectify_batch_effect(gene_counts_by_cell_line, sample_metadata_by_cell_line,
                                                                                               batch_col, covariate_cols)
            
            gene_counts_by_cell_line_rectified.append(gene_counts_by_cell_line)
            
        gene_counts_rectified = pd.concat(gene_counts_by_cell_line_rectified, axis=0)
        self.write_table(gene_counts_rectified, outpath)

    @staticmethod
    def extract_cell_lines(sample_metadata: pd.DataFrame, cell_line_col: str) -> set[str]:
        cell_lines = set(sample_metadata[cell_line_col].to_list())
        cell_lines.discard("nan")
        return cell_lines
    
    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                             cell_line: str, cell_line_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata[cell_line_col] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata
    
    @classmethod
    def rectify_batch_effect(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                             batch_col: str, covariate_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
        print("Starting batch effect correction (this may take a while)")
        batches = sample_metadata[batch_col].to_list()
        covariates = sample_metadata[covariate_cols]
        gene_counts = pycombat_seq(gene_counts.T, batches, covar_mod=covariates).T.astype("int")
        return cls.correct_negative_numbers(gene_counts), sample_metadata

    @staticmethod
    def correct_negative_numbers(gene_counts: pd.DataFrame) -> pd.DataFrame:
        negative_count = (gene_counts < 0).sum().sum()
        if negative_count > 0:
            print(f"\nWARNING: {negative_count} negative values detected after batch effect correction")
            gene_counts[gene_counts < 0] = 0
            print("Negative values have been set to 0\n")
            return gene_counts
        return gene_counts

    @staticmethod
    def write_table(table: pd.DataFrame, outpath: Path) -> None:
        table.to_csv(outpath)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_counts", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-outpath", type=str, required=True)
    parser.add_argument("-cell_line_col", type=str, default="cell_line")
    parser.add_argument("-batch_col", type=str, default="library_prep_batch")
    parser.add_argument("-covariate_cols", nargs="*", default=["treatment", "time_point"])
    args = parser.parse_args()

    gene_counts = pd.read_csv(args.gene_counts, index_col=0)
    smm = SampleMetadataManager(args.sample_metadata)

    r = BatchEffectRectifier()
    r.run(gene_counts, smm.sample_metadata, Path(args.outpath),
          args.cell_line_col, args.batch_col, args.covariate_cols)