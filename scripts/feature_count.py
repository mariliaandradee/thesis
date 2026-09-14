import os
import subprocess
import sys
import argparse
from typing import List, Dict

from conf import (
    BASE_DIR, FEATURECOUNTS_PATH, GFF_FILE,
    ALIGNMENT_DIR, FEATURECOUNTS_DIR, THREADS
)


def clean_gtf_file(input_gtf: str, output_gtf: str, log_file: str = None):
    """Clean the GTF's gene_id attribute and, along the way, keep a record
    of every line that needed fixing (line number, feature type, chromosome/
    contig, and the attributes field before and after) so the correction can
    be shown/reproduced later, e.g. for the advisor or the Supplementary
    Material, instead of having to re-scan the raw GTF separately.
    """
    fixed_examples: List[Dict[str, str]] = []
    line_number = 0
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if not line.startswith("#"):
                line_number += 1
                fields = line.strip().split('\t')
                attributes = fields[-1].split(';')
                new_attributes = []
                attributes_before = fields[-1]
                fixed_this_line = False
                for attribute in attributes:
                    if attribute.strip().startswith("gene_id"):
                        gene_id = attribute.strip().split(' ')[1].strip('"')
                        if not gene_id:
                            attribute = 'gene_id "unknown_gene"'
                            fixed_this_line = True
                    new_attributes.append(attribute)
                fields[-1] = '; '.join(new_attributes)

                if fixed_this_line:
                    fixed_examples.append({
                        "line_number": line_number,
                        "chromosome": fields[0] if len(fields) > 0 else "",
                        "feature_type": fields[2] if len(fields) > 2 else "",
                        "start": fields[3] if len(fields) > 3 else "",
                        "end": fields[4] if len(fields) > 4 else "",
                        "attributes_before": attributes_before.strip(),
                        "attributes_after": fields[-1].strip(),
                    })

                outfile.write('\t'.join(fields) + '\n')
            else:
                outfile.write(line)

    print(f"GTF cleaning: {len(fixed_examples)} exon entries had an empty "
          f"gene_id and were reassigned to 'unknown_gene' (out of {line_number} "
          f"non-comment lines).")

    if log_file:
        with open(log_file, 'w') as logf:
            logf.write("line_number,chromosome,feature_type,start,end,attributes_before,attributes_after\n")
            for ex in fixed_examples:
                logf.write(
                    f'{ex["line_number"]},{ex["chromosome"]},{ex["feature_type"]},'
                    f'{ex["start"]},{ex["end"]},"{ex["attributes_before"]}","{ex["attributes_after"]}"\n'
                )
        print(f"Full list of corrected lines saved to: {log_file}")

    return fixed_examples

def run_featurecounts(input_bams: List[str], output_file: str, annotation_file: str, threads: int):
    command = [
        FEATURECOUNTS_PATH,
        "-T", str(threads),
        "-a", annotation_file,
        "-o", output_file,
        "-t", "exon",
        "-g", "gene_id",
        "-p",
        "-P",
        "-B",
        "-C",
        *input_bams
    ]

    print(f"Running featureCounts with command: {' '.join(command)}")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("featureCounts output:")
        print(result.stdout)
        if result.stderr:
            print("featureCounts errors/warnings:")
            print(result.stderr)

        # Verificação do arquivo de resumo
        summary_file = output_file + ".summary"
        if os.path.exists(summary_file):
            print(f"Summary file created: {summary_file}")
            print("First few lines of the summary file:")
            with open(summary_file, 'r') as f:
                for _ in range(5):
                    print(f.readline().strip())
        else:
            print(f"Warning: Summary file not found at {summary_file}")

    except subprocess.CalledProcessError as e:
        print(f"Error running featureCounts: {e}")
        print("Error output:")
        print(e.stderr)
        sys.exit(1)


def parse_featurecounts_summary(summary_file: str) -> Dict[str, Dict[str, int]]:
    print(f"Parsing summary file: {summary_file}")
    summary = {}
    try:
        with open(summary_file, 'r') as f:
            all_lines = f.readlines()

            if len(all_lines) < 2:
                print("Error: Summary file has less than 2 lines")
                return summary

            headers = all_lines[0].strip().split('\t')
            if len(headers) < 2:
                print("Error: Header line does not have enough columns")
                return summary

            bam_files = headers[1:]

            for line in all_lines[1:]:
                parts = line.strip().split('\t')
                status = parts[0]

                if len(parts) - 1 != len(bam_files):
                    print(f"Warning: Skipping line due to insufficient columns: {line.strip()}")
                    continue

                for i, bam_file in enumerate(bam_files):
                    try:
                        count = int(parts[i+1])
                    except ValueError:
                        print(f"Warning: Could not convert '{parts[i+1]}' to int for {bam_file} in status {status}")

                    if bam_file not in summary:
                        summary[bam_file] = {}
                    summary[bam_file][status] = count

    except FileNotFoundError:
        print(f"Summary file not found: {summary_file}")
    except Exception as e:
        print(f"Error parsing summary file: {e}")

    return summary

def print_summary_stats(summary: Dict[str, Dict[str, int]]):
    if not summary:
        print("No summary statistics available.")
        return
    for bam_file, stats in summary.items():
        print(f"\nSummary for {bam_file}:")
        total_reads = sum(stats.values())
        for status, count in stats.items():
            percentage = (count / total_reads) * 100 if total_reads > 0 else 0
            print(f"{status}: {count} ({percentage:.2f}%)")

def check_featurecounts_output(output_file: str):
    try:
        with open(output_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                print(f"First few lines of featureCounts output:")
                for line in lines[:5]:
                    print(line.strip())
            else:
                print("featureCounts output file is empty or contains only header.")
    except FileNotFoundError:
        print(f"featureCounts output file not found: {output_file}")

def main(threads: int):
    os.makedirs(FEATURECOUNTS_DIR, exist_ok=True)

    if not os.access(FEATURECOUNTS_DIR, os.W_OK):
        print(f"Error: FEATURECOUNTS_DIR '{FEATURECOUNTS_DIR}' is not writable.")
        sys.exit(1)

    bam_files = [os.path.join(ALIGNMENT_DIR, f) for f in os.listdir(ALIGNMENT_DIR) if f.endswith('.bam')]

    if not bam_files:
        print(f"No BAM files found in {ALIGNMENT_DIR}")
        sys.exit(1)

    print(f"Found {len(bam_files)} BAM files:")
    for bam_file in bam_files:
        print(f"  - {bam_file}")

    output_file = os.path.join(FEATURECOUNTS_DIR, "featurecounts_output.txt")
    summary_file = output_file + ".summary"

    print(f"Output file path: {output_file}")
    print(f"Summary file path: {summary_file}")

    cleaned_gtf = os.path.join(FEATURECOUNTS_DIR, "cleaned_annotations.gtf")
    gene_id_fix_log = os.path.join(FEATURECOUNTS_DIR, "gene_id_fixes_log.csv")
    clean_gtf_file(GFF_FILE, cleaned_gtf, log_file=gene_id_fix_log)

    run_featurecounts(bam_files, output_file, cleaned_gtf, threads)

    check_featurecounts_output(output_file)

    if os.path.exists(summary_file):
        print(f"Summary file found: {summary_file}")
        print(f"Summary file size: {os.path.getsize(summary_file)} bytes")
        summary = parse_featurecounts_summary(summary_file)
        print_summary_stats(summary)
    else:
        print(f"Error: Summary file not found at {summary_file}")
        print("Check if featureCounts created the summary file.")
        sys.exit(1)

    print("Feature counting process completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run featureCounts on aligned BAM files.")
    parser.add_argument("-t", "--threads", type=int, default=THREADS,
                        help=f"Number of threads to use (default: {THREADS})")
    args = parser.parse_args()

    main(args.threads)
