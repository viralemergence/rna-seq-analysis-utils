from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import matplotlib.pyplot as plt # type: ignore
from os import devnull, environ
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.dds import DeseqDataSet # type: ignore

environ["NUMBA_CACHE_DIR"] = "/tmp/" # Needed for scanpy to import properly
import scanpy as sc # type: ignore
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

class PrincipalComponentAnalysisManager:
    def __init__(self, pca_outdir: Path):
        self.pca_figure_dir = pca_outdir
        sc.settings.figdir = self.pca_figure_dir

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
            treatments: list[str], control: str,
            cell_lines: list[str], covariate_cols: list[str], cell_line_col: str,
            treatment_col: str, batch_col: str) -> None:
        design_factors = covariate_cols
        viruses = treatments
        pca_color_factors = [batch_col] + covariate_cols

        for cell_line in cell_lines:
            print(f"\nStarting on cell line: {cell_line}\n----------")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = \
                self.filter_for_cell_line(gene_counts, sample_metadata, cell_line, cell_line_col)
            
            for virus in viruses:
                print(f"\nStarting on virus: {virus}")
                gene_counts_by_virus, sample_metadata_by_virus = \
                    self.filter_for_virus(gene_counts_by_cell_line, sample_metadata_by_cell_line, virus, control, treatment_col)

                print("Starting pyDEseq2 analysis")
                with catch_warnings():
                    simplefilter("ignore")
                    dds = DeseqDataSet(counts=gene_counts_by_virus,
                                    metadata=sample_metadata_by_virus,
                                    design_factors=design_factors)
                with suppress_stdout_stderr(): # NOTE: Primarily suppresses redundant messages in prod, but could suppress actual errors important for dev work
                    dds.deseq2()
                
                dds.obs[batch_col] = dds.obs[batch_col].astype(int).astype(str)
                self.perform_principal_component_analysis(dds, cell_line, virus, pca_color_factors)

    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                             cell_line: str, cell_line_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata[cell_line_col] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata
    
    @staticmethod
    def filter_for_virus(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame,
                         virus: str, control: str, treatment_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        virus_list = [virus, control]
        sample_metadata = sample_metadata[sample_metadata[treatment_col].isin(virus_list)]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata

    @staticmethod
    def perform_principal_component_analysis(dds: DeseqDataSet, cell_line: str, virus: str, design_factors: list[str]) -> None:
        plt.rcParams["svg.fonttype"] = "none"
        sc.tl.pca(dds)
        parameters = {"size": 200, "annotate_var_explained": True}
        for design_factor in design_factors:
            color = design_factor
            if design_factor == "time_point":
                color = design_factor.replace("_", "-")
            sc.pl.pca(dds, color=color, save=f"_{cell_line}_{virus}_{design_factor}.svg", **parameters)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-gene_counts", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-cell_lines", nargs="+", required=True)
    parser.add_argument("-treatments", nargs="+", required=True)
    parser.add_argument("-control", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-covariate_cols", nargs="*", default=["treatment", "time_point"])
    parser.add_argument("-cell_line_col", type=str, default="cell_line")
    parser.add_argument("-treatment_col", type=str, default="treatment")
    parser.add_argument("-batch_col", type=str, default="library_prep_batch")
    args = parser.parse_args()

    gene_counts = pd.read_csv(args.gene_counts, index_col=0).sort_index()
    smm = SampleMetadataManager(args.sample_metadata)
    
    pcam = PrincipalComponentAnalysisManager(Path(args.outdir))
    pcam.run(gene_counts, smm.sample_metadata,
             args.treatments, args.control,
             args.cell_lines, args.covariate_cols, args.cell_line_col,
             args.treatment_col, args.batch_col)