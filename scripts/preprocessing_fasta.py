import sys
import os
import subprocess
from conf import BASE_DIR, PREPROCESSING_DIR, THREADS, TRIMMOMATIC_PATH, TRUQUE3_PATH

# "java" tem de estar instalado e no PATH do teu Mac (confirma com
# "java -version"; se nao tiveres, "brew install openjdk" e depois segue
# as instrucoes do brew para o pores no PATH)
java_path = "java"

def trimming_plots(forward_file, reverse_file, output_forward, output_reverse):
    unpaired_forward = os.path.join(os.path.dirname(output_forward), f"unpaired_{os.path.basename(forward_file)}")
    unpaired_reverse = os.path.join(os.path.dirname(output_reverse), f"unpaired_{os.path.basename(reverse_file)}")

    trimmomatic_jar = TRIMMOMATIC_PATH

    try:
        trimmomatic_command = [
            java_path, "-jar", trimmomatic_jar, "PE",
            "-threads", str(THREADS),
            "-phred33",
            forward_file, reverse_file,
            output_forward, unpaired_forward,
            output_reverse, unpaired_reverse,
            f"ILLUMINACLIP:{TRUQUE3_PATH}:2:30:10:2:True", "LEADING:3", "TRAILING:3",
            "SLIDINGWINDOW:4:15", "MINLEN:36"
        ]

        print("Running Trimmomatic with command:", " ".join(trimmomatic_command))
        subprocess.run(trimmomatic_command, check=True)
        print(f"Trimming completed for {forward_file} and {reverse_file}")
    except subprocess.CalledProcessError as e:
        print("Error with Trimmomatic:", e)
        sys.exit(1)

def process_files_in_directories(directories):
    output_dir = PREPROCESSING_DIR
    os.makedirs(output_dir, exist_ok=True)

    for directory in directories:
        if not os.path.isdir(directory):
            print(f"Directory not found: {directory}")
            continue

        print(f"Processing files in directory: {directory}")
        files = os.listdir(directory)
        if not files:
            print(f"No files found in the directory: {directory}")
            continue

        for filename in files:
            if filename.endswith("_R1_001.fastq.gz"):
                forward_file = os.path.join(directory, filename)
                reverse_file = os.path.join(directory, filename.replace("R1", "R2"))

                if os.path.exists(reverse_file):
                    output_forward = os.path.join(output_dir, f"trimmed_{os.path.basename(forward_file)}")
                    output_reverse = os.path.join(output_dir, f"trimmed_{os.path.basename(reverse_file)}")

                    if os.path.exists(output_forward) and os.path.exists(output_reverse):
                        print(f"Files already processed: {output_forward}, {output_reverse}")
                        continue

                    print(f"Processing pair: {forward_file}, {reverse_file}")

                    try:
                        trimming_plots(forward_file, reverse_file, output_forward, output_reverse)

                        if os.path.exists(output_forward) and os.path.exists(output_reverse):
                            print(f"Files created successfully: {output_forward}, {output_reverse}")
                        else:
                            print("Error: Output files were not created.")
                            raise FileNotFoundError("Output files were not created.")
                    except Exception as e:
                        print(f"An error occurred while processing {forward_file} and {reverse_file}: {e}")
                else:
                    print(f"Reverse file not found for {forward_file}")

if __name__ == "__main__":
    # ATENCAO: preenche com as pastas locais que tem os FASTQ em bruto
    # (antes do Trimmomatic), no teu Mac. Os 2 caminhos do cluster que
    # estavam aqui antes nao existem neste computador.
    # Nota: como o codigo original juntava BASE_DIR a um caminho ja
    # absoluto, o os.path.join() acima ignorava sempre o BASE_DIR (nao
    # tinha efeito nenhum) - por isso os caminhos aqui devem ja vir
    # completos, sem depender do join.
    sample_directories = [
        # os.path.join(BASE_DIR, "nome_da_pasta_com_os_fastq_em_bruto"),
    ]

    if not sample_directories:
        print("Nenhuma pasta de amostras em bruto configurada em sample_directories — "
              "edita este ficheiro e adiciona os caminhos locais antes de correr.")
        sys.exit(1)

    print(f"Starting processing in directories: {', '.join(sample_directories)}")
    process_files_in_directories(sample_directories)
