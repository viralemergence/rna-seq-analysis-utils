from argparse import ArgumentParser
from collections import defaultdict
from csv import DictWriter, reader, writer
import pandas as pd # type: ignore
from pathlib import Path
from typing import Any, Iterator

class SampleMetadataManager:
    def __init__(self, sample_metadata_path: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata_path)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> list[dict[Any]]:
        data = pd.read_csv(sample_metadata_path, dtype={"uid": str, "technical_replicate_parent_uid": str})
        return data.to_dict(orient="records")

class StarGeneCountManager:
    def __init__(self, input_dir: Path, output_file: Path) -> None:
        self.input_files = self.get_input_files(input_dir)
        self.output_file = output_file

    @staticmethod
    def get_input_files(input_dir: Path) -> list[Path]:
        files = sorted([file for file in input_dir.iterdir() if file.is_file()])
        return files

    def run(self, sample_metadata: list[dict[Any]]) -> None:
        collated_gene_counts = self.get_collated_gene_counts(self.input_files)
        gene_name_union = self.get_gene_name_union(collated_gene_counts)
        self.fill_missing_zeros(collated_gene_counts, gene_name_union)

        gene_counts = self.convert_gene_counts_to_pd_dataframe(collated_gene_counts, gene_name_union)
        gene_counts = self.combine_technical_replicates(gene_counts, sample_metadata)
        gene_counts = self.remove_low_count_rows(gene_counts)

        gene_counts = self.transpose_gene_counts(gene_counts)
        gene_counts = self.sort_by_index(gene_counts)
        self.write_table(gene_counts, self.output_file)

    @classmethod
    def get_collated_gene_counts(cls, input_files: list[Path]) -> list[defaultdict[int]]:
        collated_gene_counts = list()
        for file in input_files:
            sample_info = defaultdict(int)
            sample_info["sample_uid"] = cls.extract_sample_uid(file)
            sample_info.update(cls.extract_gene_counts(file))
            collated_gene_counts.append(sample_info)
        return collated_gene_counts
        
    @staticmethod
    def extract_sample_uid(file: Path) -> str:
        suffix = "_ReadsPerGene.out.tab"
        return file.name.replace(suffix, "")
    
    @classmethod
    def extract_gene_counts(cls, file: Path) -> dict[int]:
        gene_counts = dict()
        with file.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            cls.skip_gene_count_header(reader_iterator)
            for line in reader_iterator:
                gene = line[0]
                strand_1_count = int(line[2])
                strand_2_count = int(line[3])
                total_count = strand_1_count + strand_2_count

                gene_counts[gene] = total_count
        return gene_counts

    @staticmethod
    def skip_gene_count_header(reader_iterator: Iterator) -> None:
        header_set = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}
        while len(header_set) > 0:
            header = next(reader_iterator)[0]
            header_set.remove(header)
            
    @staticmethod
    def get_gene_name_union(collated_gene_counts: list[defaultdict[int]]) -> list[str]:
        gene_names = set()
        for sample in collated_gene_counts:
            for gene, _ in sample.items():
                gene_names.add(gene)
        gene_names.remove("sample_uid")
        return sorted(list(gene_names))

    @staticmethod
    def fill_missing_zeros(collated_gene_counts: list[defaultdict[int]], gene_name_union: list[str]) -> None:
        for sample in collated_gene_counts:
            for gene in gene_name_union:
                sample[gene] = sample[gene]

    @classmethod
    def convert_gene_counts_to_pd_dataframe(cls, gene_counts: list[defaultdict[int]], gene_name_union: list[str]) -> pd.DataFrame:
        all_keys = ["sample_uid"] + gene_name_union

        data = [all_keys]
        for sample_data in gene_counts:
            counts = [sample_data[key] for key in all_keys]
            data.append(counts)

        gene_counts = cls.transpose_list_of_lists(data)
        gene_counts = pd.DataFrame(gene_counts[1:], columns=gene_counts[0])
        return cls.set_gene_id_as_index(gene_counts)

    @staticmethod
    def transpose_list_of_lists(data: list[list]) -> list[list]:
        return [list(i) for i in zip(*data)]

    @staticmethod
    def set_gene_id_as_index(gene_counts: pd.DataFrame) -> pd.DataFrame:
        gene_counts = gene_counts.rename(columns={"sample_uid": "gene_id"})
        return gene_counts.set_index("gene_id")

    @staticmethod
    def combine_technical_replicates(gene_counts: pd.DataFrame, sample_metadata: list[dict[Any]]) -> pd.DataFrame:
        for data in sample_metadata:
            uid = data["uid"]
            parent_uid = data["technical_replicate_parent_uid"]
            if pd.notna(parent_uid):
                gene_counts[parent_uid] = gene_counts[parent_uid] + gene_counts[uid]
                gene_counts = gene_counts.drop([uid], axis=1)
        return gene_counts

    @staticmethod
    def remove_low_count_rows(gene_counts: pd.DataFrame, threshold: int = 10) -> pd.DataFrame:
        return gene_counts[gene_counts.sum(axis = 1) > threshold]

    @staticmethod
    def transpose_gene_counts(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.T

    @staticmethod
    def sort_by_index(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.sort_index()

    @staticmethod
    def write_table(table: pd.DataFrame, outpath: Path) -> None:
        table.to_csv(outpath)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-input_dir", type=str, required=True)
    parser.add_argument("-collated_gene_counts", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sample_metadata))

    sgcm = StarGeneCountManager(Path(args.input_dir), Path(args.collated_gene_counts))
    sgcm.run(smm.sample_metadata)