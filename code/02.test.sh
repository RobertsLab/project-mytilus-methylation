# Set directories and files
reads_dir="data/raw-wgbs/"
genome_folder="output/01-bismark-init/"
output_dir="output/02-bismark-klone-array/"
checkpoint_file="output/02-bismark-klone-array/completed_samples.log"

# Print variables to check their values
echo "reads_dir: $reads_dir"
echo "genome_folder: $genome_folder"
echo "output_dir: $output_dir"
echo "checkpoint_file: $checkpoint_file"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}*_1.fastq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_1.fastq.gz")

# Print the list of files and selected file and sample name
echo "files: ${files[@]}"
echo "file: $file"
echo "sample_name: $sample_name"

# Check if the sample has already been processed
if grep -q "^${sample_name}$" ${checkpoint_file}; then
    echo "Sample ${sample_name} already processed. Skipping..."
    exit 0
fi

# Define log files for stdout and stderr
stdout_log="${output_dir}${sample_name}_stdout.log"
stderr_log="${output_dir}${sample_name}_stderr.log"

# Print log file paths
echo "stdout_log: $stdout_log"
echo "stderr_log: $stderr_log"