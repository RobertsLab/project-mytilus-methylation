#!/bin/bash
# Set variables
reads_dir="../data/raw-wgbs/"
#bismark_dir="/path/to/bismark/"
#bowtie2_dir="/path/to/bowtie2/"
genome_folder="../output/01-bismark-init/"
output_dir="../output/01-bismark-init/"
checkpoint_file="../output/01-bismark-init/completed_samples.log"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Run Bismark for each sample
for file in ${reads_dir}*_1.fastq.gz; do
    sample_name=$(basename "$file" _1.fastq.gz)

    # Check if the sample has already been processed
    if grep -q "^${sample_name}$" ${checkpoint_file}; then
        echo "Sample ${sample_name} already processed. Skipping..."
        continue
    fi

    # Define log files for stdout and stderr
    stdout_log="${output_dir}${sample_name}_stdout.log"
    stderr_log="${output_dir}${sample_name}_stderr.log"

    # Run Bismark and redirect stdout and stderr
    bismark \
    -genome ${genome_folder} \
    -p 8 \
    -score_min L,0,-0.6 \
    --non_directional \
    -1 ${reads_dir}${sample_name}_1.fastq.gz \
    -2 ${reads_dir}${sample_name}_2.fastq.gz \
    -o ${output_dir} \
    > ${stdout_log} 2> ${stderr_log}

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        # Append the sample name to the checkpoint file
        echo ${sample_name} >> ${checkpoint_file}
        echo "Sample ${sample_name} processed successfully."
    else
        echo "Sample ${sample_name} failed. Check ${stderr_log} for details."
        break
    fi
done