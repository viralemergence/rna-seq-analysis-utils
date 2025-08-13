from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import matplotlib.pyplot as plt # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.dds import DeseqDataSet # type: ignore
import seaborn as sns # type: ignore
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

class GeneCountCorrelations:
    def __init__(self, cell_lines: list[str], outdir: Path) -> None:
        self.cell_lines = cell_lines
        self.outdir = outdir

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
            cell_line_col: str, covariate_cols: list[str], custom_order: list[str]) -> None:
        design_factors = covariate_cols

        for cell_line in self.cell_lines:
            print(f"\nStarting on cell line: {cell_line}")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata,
                                                                                               cell_line, cell_line_col)
            dds = self.pydeseq2_normalization(gene_counts_by_cell_line, sample_metadata_by_cell_line, design_factors)
            
            normalized_counts_by_cell_line = self.extract_normalized_count_df_from_dds(dds)

            normalized_counts = self.transform_normalized_counts(normalized_counts_by_cell_line, dds, design_factors)
            normalized_counts = self.remove_low_gene_counts(normalized_counts)
            normalized_counts = self.reformat_dataframe(normalized_counts, custom_order)
            
            method = "spearman"
            correlation_matrix = self.calculate_correlation_matrix(normalized_counts, method)

            self.generate_heatmap(correlation_matrix, cell_line, method, self.outdir)
    
    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                             cell_line: str, cell_line_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata[cell_line_col] == cell_line]
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
    def extract_normalized_count_df_from_dds(dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts = dds.layers["normed_counts"]
        gene_ids = dds._var.index.to_list()
        sample_ids = dds.obsm["design_matrix"].index.to_list()
        return pd.DataFrame(normalized_counts, index=sample_ids, columns=gene_ids)

    @staticmethod
    def transform_normalized_counts(normalized_counts: pd.DataFrame, dds: DeseqDataSet, design_factors: list[str]) -> pd.DataFrame:
        # Depending on experimental design, this may need to be changed
        normalized_counts = normalized_counts.copy()
        if "time_point" in design_factors:
            normalized_counts["time_point"] = dds.obs["time-point"].astype(float).astype(int)
            # Note pydeseq2 does a dumb thing where it makes the underscore a hyphen, so this has to be corrected
        if "treatment" in design_factors:
            normalized_counts["treatment"] = dds.obs["treatment"]
        grouped_counts = normalized_counts.groupby(design_factors)
        return grouped_counts.mean().round(2)

    @staticmethod
    def remove_low_gene_counts(normalized_counts: pd.DataFrame) -> pd.DataFrame:
        per_gene_stats = normalized_counts.aggregate(["mean", "std"], axis=0).round(2)

        drop_genes = []
        for gene in per_gene_stats:
            if per_gene_stats.loc["mean", gene] == 0:
                drop_genes.append(gene)
        return normalized_counts.drop(drop_genes, axis=1)

    @staticmethod
    def reformat_dataframe(normalized_counts: pd.DataFrame, custom_order: list[str]) -> pd.DataFrame:
        normalized_counts = normalized_counts.sort_index(key=lambda x: x.map(dict(zip(custom_order, range(len(custom_order))))))
        return normalized_counts

    @staticmethod
    def calculate_correlation_matrix(normalized_counts: pd.DataFrame, method="spearman") -> pd.DataFrame:
        correlation_matrix = normalized_counts.T.corr(method=method)
        correlation_matrix = correlation_matrix.apply(lambda x: round(x, 4))
        max_value = correlation_matrix.max().max()
        min_value = correlation_matrix.min().min()
        print(f"{max_value}, {min_value}")
        return correlation_matrix

    @staticmethod
    def generate_heatmap(correlation_matrix: pd.DataFrame, cell_line: str, method: str, outdir: Path) -> None:
        plt.rcParams["svg.fonttype"] = "none"
        sns.set_theme(font="sans-serif", font_scale=0.6)
        fig, ax = plt.subplots(figsize=(2.5, 2.5), dpi=1200)
    
        heatmap = sns.heatmap(correlation_matrix, ax=ax, vmin=0.9, square=True, cbar_kws={"shrink": 0.75})
        cbar = heatmap.collections[0].colorbar
        cbar.set_label("Spearman Coefficient", labelpad=10, fontweight="bold")

        columns = [4, 8]
        line_width = 2
        for col in columns:
            heatmap.axvline(col, color="white", lw=line_width)
            heatmap.axhline(col, color="white", lw=line_width)
        heatmap.set_xlabel("")
        heatmap.set_ylabel("")
        
        heatmap.tick_params(left=False, bottom=False, pad=0)

        figure_filename = f"{cell_line}_{method}_heatmap.svg"
        figure_outpath = outdir / figure_filename
        fig.savefig(figure_outpath, bbox_inches="tight")
        plt.close(fig)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_counts", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-cell_lines", nargs="+", required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-cell_line_col", type=str, default="cell_line")
    parser.add_argument("-covariate_cols", nargs="*", default=["treatment", "time_point"])
    parser.add_argument("-custom_order", nargs="*", default=["No-Virus", "MR766", "PRVABC59"])

    args = parser.parse_args()

    gene_counts = pd.read_csv(args.gene_counts, index_col=0).sort_index()
    smm = SampleMetadataManager(args.sample_metadata)
    
    gcc = GeneCountCorrelations(args.cell_lines, Path(args.outdir))
    gcc.run(gene_counts, smm.sample_metadata,
            args.cell_line_col, args.covariate_cols, args.custom_order)