import subprocess
import os
import sys
from typing import Dict, List
import shutil

from conf import (
    BASE_DIR, HISAT2_PATH, SAMTOOLS_PATH,
    GENOME_FASTA, GFF_FILE, PREPROCESSING_DIR,
    ALIGNMENT_DIR, INDEX_BASE_NAME, THREADS, HISAT2_BUILD_PATH, HISAT2_DIR
)

def extract_splice_sites_and_exons(gff_file: str, splice_sites_file: str, exons_file: str) -> None:
    for script, output_file in [
        ("hisat2_extract_splice_sites.py", splice_sites_file),
        ("hisat2_extract_exons.py", exons_file)
    ]:
        script_path = os.path.join(HISAT2_DIR, script)
        print(f"Extracting from {gff_file} using {script}")
        with open(output_file, 'w') as outfile:
            subprocess.run([sys.executable, script_path, gff_file], stdout=outfile, check=True)

def build_hisat2_index(genome_fasta: str, index_base_name: str, splice_sites_file: str, exons_file: str) -> None:
    os.makedirs(os.path.dirname(index_base_name), exist_ok=True)

    # NAO passar --ss/--exon aqui: isso faz o hisat2-build construir o
    # indice grafico (HGFM), que para um genoma deste tamanho precisa de
    # muito mais RAM do que um portatil tem (e foi o que matou o processo
    # com SIGKILL). Em vez disso, indexa-se so o genoma (leve), e os
    # splice sites conhecidos sao passados ao HISAT2 na hora do
    # alinhamento via --known-splicesite-infile (ve run_hisat2_alignment).
    command = [
        HISAT2_BUILD_PATH, "-p", str(THREADS), genome_fasta, index_base_name
    ]

    print(f"Running HISAT2-build command: {' '.join(command)}")
    for file in [genome_fasta, splice_sites_file, exons_file]:
        print(f"{os.path.basename(file)}: {'Exists' if os.path.exists(file) else 'Not found'}")

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"HISAT2-build error (return code {e.returncode}):")
        print(e.stderr)
        raise

def parse_hisat2_stats(stderr_output: str) -> Dict[str, float]:
    stats = {}
    for line in stderr_output.split('\n'):
        if 'aligned 0 times' in line:
            stats['unaligned'] = int(line.split()[0])
        elif 'aligned exactly 1 time' in line:
            stats['uniquely_aligned'] = int(line.split()[0])
        elif 'aligned >1 times' in line:
            stats['multi_aligned'] = int(line.split()[0])
        elif 'overall alignment rate' in line:
            stats['overall_alignment_rate'] = float(line.split('%')[0])
    return stats

def run_hisat2_alignment(index_base_name: str, reads: List[str], output_bam: str, strandness: str = "unstranded", splice_sites_file: str = None, **kwargs) -> None:
    hisat2_exec = HISAT2_PATH
    os.makedirs(ALIGNMENT_DIR, exist_ok=True)

    output_sam = output_bam.replace('.bam', '.sam')
    hisat2_command = [
        hisat2_exec, "-x", index_base_name,
        "--new-summary", "--summary-file", f"{output_bam}.hisat2_summary.txt",
        "-p", str(THREADS), "-S", output_sam
    ]

    # Como o indice foi construido sem --ss/--exon (para poupar RAM na
    # construcao), os splice sites conhecidos entram aqui, no momento do
    # alinhamento, em vez de estarem embutidos no indice.
    if splice_sites_file:
        hisat2_command.extend(["--known-splicesite-infile", splice_sites_file])

    hisat2_command.extend(["-1", reads[0], "-2", reads[1]] if len(reads) == 2 else ["-U", reads[0]])

    if strandness == "forward":
        hisat2_command.extend(["--rna-strandness", "FR"])
    elif strandness == "reverse":
        hisat2_command.extend(["--rna-strandness", "RF"])

    hisat2_command.extend([f"--{key}", str(value)] for key, value in kwargs.items())
    print(f"Executing HISAT2 command: {' '.join(hisat2_command)}")

    samtools_view_command = [SAMTOOLS_PATH, "view", "-bS", "-h", output_sam, "-o", output_bam]
    samtools_sort_command = [SAMTOOLS_PATH, 'sort', '-m', '1G', '-@', '5', '-o', output_bam.replace('.bam', '_sorted.bam'), output_bam]
    print(f"Executing pipeline:")
    print(f"HISAT2 command: {' '.join(hisat2_command)}")
    print(f"Samtools view command: {' '.join(samtools_view_command)}")
    print(f"Samtools sort command: {' '.join(samtools_sort_command)}")

    try:
        hisat2_process = subprocess.Popen(hisat2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, hisat2_stderr = hisat2_process.communicate()
        if hisat2_process.returncode != 0:
            print(f"HISAT2 error (return code {hisat2_process.returncode}):")
            print(hisat2_stderr.decode())
            return

        samtools_view_process = subprocess.Popen(samtools_view_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, samtools_view_stderr = samtools_view_process.communicate()
        if samtools_view_process.returncode != 0:
            print(f"Samtools view error (return code {samtools_view_process.returncode}):")
            print(samtools_view_stderr.decode())
            return

        samtools_sort_process = subprocess.Popen(samtools_sort_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, samtools_sort_stderr = samtools_sort_process.communicate()
        if samtools_sort_process.returncode != 0:
            print(f"Samtools sort error (return code {samtools_sort_process.returncode}):")
            print(samtools_sort_stderr.decode())
            return

        print("HISAT2 Alignment Statistics:")
        alignment_stats = parse_hisat2_stats(hisat2_stderr.decode())
        for key, value in alignment_stats.items():
            print(f"{key}: {value}")

        save_alignment_stats(output_bam, alignment_stats, hisat2_stderr.decode())

        sorted_bam = output_bam.replace('.bam', '_sorted.bam')
        if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 0:
            os.rename(sorted_bam, output_bam)
            print(f"Finished alignment and sorting: {output_bam}")
            print(f"BAM file size: {os.path.getsize(output_bam)} bytes")
        else:
            print(f"Warning: Output file {sorted_bam} is missing or empty.")

        if os.path.exists(output_sam):
            os.remove(output_sam)

    except Exception as e:
        print(f"An error occurred: {str(e)}")

def is_genome_indexed(index_base_name: str) -> bool:
    return all(os.path.isfile(f"{index_base_name}.{i}.ht2") for i in range(1, 9))

def check_required_files() -> None:
    required_files = [
        HISAT2_BUILD_PATH,
        HISAT2_PATH,
        os.path.join(HISAT2_DIR, "hisat2_extract_splice_sites.py"),
        os.path.join(HISAT2_DIR, "hisat2_extract_exons.py")
    ]

    for file in required_files:
        if not os.path.isfile(file):
            print(f"Error: Required HISAT2 file '{os.path.basename(file)}' not found in {os.path.dirname(file)}")
            sys.exit(1)

def align_reads(preprocessing_dir: str, alignment_dir: str, index_base_name: str, splice_sites_file: str = None) -> None:
    print(f"Searching for FASTQ files in: {preprocessing_dir}")

    fastq_files = [f for f in os.listdir(preprocessing_dir) if f.startswith("trimmed_") and f.endswith("_R1.gz")]
    print(f"Found {len(fastq_files)} R1 FASTQ files")

    for filename in fastq_files:
        forward_reads = os.path.join(preprocessing_dir, filename)
        reverse_filename = filename.replace("_R1.gz", "_R2.gz")
        reverse_reads = os.path.join(preprocessing_dir, reverse_filename)

        print(f"Processing: {filename}")
        print(f"Forward reads: {forward_reads}")
        print(f"Reverse reads: {reverse_reads}")

        if not os.path.exists(reverse_reads):
            print(f"Warning: Paired file {reverse_reads} not found. Skipping {filename}")
            continue

        base_name = filename.replace("trimmed_", "").replace("_R1.gz", "")
        output_bam = os.path.join(alignment_dir, f"{base_name}_aligned_reads_hisat2.bam")

        if os.path.exists(output_bam):
            print(f"File {output_bam} already exists. Skipping.")
            continue

        print(f"Aligning {forward_reads} and {reverse_reads}")
        print(f"Output BAM file: {output_bam}")
        run_hisat2_alignment(INDEX_BASE_NAME, [forward_reads, reverse_reads], output_bam, strandness="unstranded", splice_sites_file=splice_sites_file)

    print("Finished alignment.")

def main():
    print(f"Current working directory: {os.getcwd()}")
    print(f"HISAT2_PATH: {HISAT2_PATH}")
    print(f"HISAT2_BUILD_PATH: {HISAT2_BUILD_PATH}")
    print(f"ALIGNMENT_DIR: {ALIGNMENT_DIR}")
    print(f"INDEX_BASE_NAME: {INDEX_BASE_NAME}")

    os.makedirs(ALIGNMENT_DIR, exist_ok=True)
    print(f"Alignment directory created/exists: {os.path.exists(ALIGNMENT_DIR)}")

    check_required_files()

    splice_sites_file = os.path.join(BASE_DIR, "splice_sites.txt")
    exons_file = os.path.join(BASE_DIR, "exons.txt")

    if not os.path.isfile(splice_sites_file) or not os.path.isfile(exons_file):
        extract_splice_sites_and_exons(GFF_FILE, splice_sites_file, exons_file)

    if not is_genome_indexed(INDEX_BASE_NAME):
        build_hisat2_index(GENOME_FASTA, INDEX_BASE_NAME, splice_sites_file, exons_file)
    else:
        print("The genome is already indexed")

    print("Starting alignment...")
    align_reads(PREPROCESSING_DIR, ALIGNMENT_DIR, INDEX_BASE_NAME, splice_sites_file=splice_sites_file)

if __name__ == "__main__":
    main()
