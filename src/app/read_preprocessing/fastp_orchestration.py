from argparse import ArgumentParser
import pandas as pd # type: ignore
from pathlib import Path
import subprocess
from typing import Any

class SampleMetadataManager:
    def __init__(self, sample_metadata_path: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata_path)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> list[dict[Any]]:
        return pd.read_csv(sample_metadata_path).to_dict(orient="records")

class FastqInputManager:
    def __init__(self, input_dir: Path, output_dir: Path, sample_metadata: list[dict[Any]], array_number: int) -> None:
        input_fastq_paths = self.extract_input_fastq_paths(input_dir, sample_metadata, array_number)
        self.output_fastq_paths = self.set_output_fastq_paths(input_fastq_paths, output_dir)

        self.input_fastq_paths = input_fastq_paths[1:]

    @staticmethod
    def extract_input_fastq_paths(input_dir: Path, sample_metadata: list[dict[Any]], array_number: int) -> list[Any]:
        sample_data = sample_metadata[array_number]
        print(sample_data)

        uid = sample_data["uid"]
        r1_name = sample_data["raw_r1_file"]
        r2_name = r1_name.replace("_R1_", "_R2_")

        r1_path = input_dir / r1_name
        r2_path = input_dir / r2_name
        return [uid, r1_path, r2_path]

    @staticmethod
    def set_output_fastq_paths(input_fastq_paths: list[Any], output_dir: Path) -> list[Path]:
        uid = input_fastq_paths[0]
        r1_path = output_dir / f"{uid}_r1_fastp.fastq.gz"
        r2_path = output_dir / f"{uid}_r2_fastp.fastq.gz"
        return [r1_path, r2_path]

class FastpManager:
    def __init__(self) -> None:
        pass

    def run_fastp(self, input_fastq_paths: list[Path], output_fastq_paths: list[Path],
                  r1_adapter: str, r2_adapter: str) -> None:
        input_r1 = input_fastq_paths[0]
        input_r2 = input_fastq_paths[1]
        output_r1 = output_fastq_paths[0]
        output_r2 = output_fastq_paths[1]

        fastp_command = ["fastp", "-V",
                         "-i", input_r1, "-I", input_r2,
                         "-o", output_r1, "-O", output_r2,
                         "--adapter_sequence", r1_adapter,
                         "--adapter_sequence_r2", r2_adapter]

        p = subprocess.Popen(fastp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stderr.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            raise Exception("Fastp did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-input_dir", type=str, required=True)
    parser.add_argument("-output_dir", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-array_number", type=int, required=True)
    parser.add_argument("-r1_adapter", type=str, required=True)
    parser.add_argument("-r2_adapter", type=str, required=True)
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sample_metadata))
    fim = FastqInputManager(Path(args.input_dir), Path(args.output_dir), smm.sample_metadata, args.array_number)
    fpm = FastpManager()
    fpm.run_fastp(fim.input_fastq_paths, fim.output_fastq_paths, args.r1_adapter, args.r2_adapter)