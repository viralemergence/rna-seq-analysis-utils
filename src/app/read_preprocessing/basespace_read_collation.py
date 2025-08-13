from argparse import ArgumentParser
from pathlib import Path
from shutil import move

class FastqManager:
    def __init__(self, basespace_dir: Path, collated_dir: Path) -> None:
        self.basespace_dir = basespace_dir
        self.collated_dir = collated_dir

    def run(self) -> None:
        fastq_paths = self.get_fastq_paths(self.basespace_dir)
        self.move_fastqs(fastq_paths, self.collated_dir)

    @staticmethod
    def get_fastq_paths(basespace_dir: Path) -> set[Path]:
        fastq_paths = set()

        subdirectories = sorted([subdir for subdir in basespace_dir.iterdir() if subdir.is_dir()])
        for subdirectory in subdirectories:
            fastq_files = [fastq for fastq in subdirectory.iterdir() if fastq.is_file()]
            
            if len(fastq_files) > 2:
                raise Exception("More than 2 files detected in subdirectory of interest")
        
            for fastq in fastq_files:
                if fastq.stem.endswith("_R1_001.fastq"):
                    fastq_paths.add(fastq)
                elif fastq.stem.endswith("_R2_001.fastq"):
                    fastq_paths.add(fastq)
                else:
                    print(f"No fastq detected in directory: {subdirectory}")

        return fastq_paths
        
    @staticmethod
    def move_fastqs(fastq_paths: set[Path], collated_dir: Path) -> None:
        for fastq_path in fastq_paths:
            move(fastq_path, collated_dir)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-basespace_dir", type=str, required=True)
    parser.add_argument("-collated_dir", type=str, required=True)
    args = parser.parse_args()

    fm = FastqManager(Path(args.basespace_dir), Path(args.collated_dir))
    fm.run()