from argparse import ArgumentParser
import pandas as pd # type: ignore
from pathlib import Path
from shutil import move, rmtree
import subprocess
from typing import Any

class SampleMetadataManager:
    def __init__(self, sample_metadata_path: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata_path)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> list[dict[Any]]:
        return pd.read_csv(sample_metadata_path).to_dict(orient="records")

class FilePathManager:
    def __init__(self, input_dir: Path, output_dir: Path, genomes_dir: Path, sample_metadata: list[dict[Any]], array_number: int) -> None:
        input_fastq_paths = self.extract_input_fastq_paths(input_dir, sample_metadata, array_number)

        self.output_count_dir = self.set_output_count_dir(input_fastq_paths[0], output_dir)
        self.input_fastq_paths = input_fastq_paths[1:]
        self.genome_index_path = self. extract_genome_index_path(genomes_dir, sample_metadata, array_number)

    @staticmethod
    def extract_input_fastq_paths(input_dir: Path, sample_metadata: list[dict[Any]], array_number: int) -> list[Any]:
        sample_data = sample_metadata[array_number]
        print(sample_data)

        uid = sample_data["uid"]
        r1_name = f"{uid}_r1_fastp.fastq.gz"
        r2_name = f"{uid}_r2_fastp.fastq.gz"

        r1_path = input_dir / r1_name
        r2_path = input_dir / r2_name
        return [uid, r1_path, r2_path]

    @staticmethod
    def set_output_count_dir(uid: int, output_dir: Path) -> list[Path]:
        output_count_dir = output_dir / (f"{uid}/")
        output_count_dir.mkdir(exist_ok=True)
        return output_count_dir

    @staticmethod
    def extract_genome_index_path(genomes_dir: Path, sample_metadata: list[dict[Any]], array_number: int) -> Path:
        sample_data = sample_metadata[array_number]
        genome = sample_data["genome"]
        genome_index_path = genomes_dir / (f"{genome}/")
        return genome_index_path

class StarManager:
    def __init__(self, input_fastq_paths: list[Path], genome_index_path: Path,
                 output_count_dir: Path, star_counts_dir: Path) -> None:
        self.input_fastq_paths = input_fastq_paths
        self.genome_index_path = genome_index_path
        self.output_count_dir = output_count_dir
        self.star_counts_dir = star_counts_dir

    def run_star_gene_count(self) -> None:
        input_r1 = self.input_fastq_paths[0]
        input_r2 = self.input_fastq_paths[1]

        star_command = ["STAR",
                         "--genomeDir", self.genome_index_path, "--runThreadN", "10",
                         "--readFilesCommand", "zcat",
                         "--readFilesIn", input_r1, input_r2,
                         "--outFileNamePrefix", f"{self.output_count_dir}/",
                         "--quantMode", "GeneCounts"
                         ]

        p = subprocess.Popen(star_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            print(p.stderr.readlines())
            raise Exception("STAR did not complete successfully")
        
    def move_gene_counts(self) -> None:
        gene_count_path = self.output_count_dir / "ReadsPerGene.out.tab"
        new_file_path = self.star_counts_dir / (f"{self.output_count_dir.name}_ReadsPerGene.out.tab")
        move(gene_count_path, new_file_path)

    def remove_star_files(self) -> None:
        rmtree(self.output_count_dir)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-input_dir", type=str, required=True)
    parser.add_argument("-genomes_dir", type=str, required=True)
    parser.add_argument("-star_counts_dir", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-array_number", type=int, required=True)
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sample_metadata))
    fpm = FilePathManager(Path(args.input_dir), Path(args.star_counts_dir), Path(args.genomes_dir),
                          smm.sample_metadata, args.array_number)

    sm = StarManager(fpm.input_fastq_paths, fpm.genome_index_path, fpm.output_count_dir, Path(args.star_counts_dir))
    sm.run_star_gene_count()
    sm.move_gene_counts()
    sm.remove_star_files()