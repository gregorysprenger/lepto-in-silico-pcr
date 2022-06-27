#!/usr/bin/bash -l

# Assign Job-Name instead of defauly which is the name of job-script
#$ -N download_genomes
# Submit job to the default queue "short"
#$ -q all.q
# Start the script in the current working directory
#$ -cwd
# Email notifications
#$ -M $USER
# Email notifications if the job aborts, at beginning of job and upon exit
#$ -m abe


# Create conda env and dependency files
# -------------------------------------
# conda create -n genomedl
# conda activate genomedl
# mamba install ncbi-genome-download


conda activate genomedl


# Pass all taxids from above to ncbi-genome-downloader
ngd --genera "Leptospira" --formats fasta --human-readable -o from_NCBI --verbose bacteria --parallel 5 --retries 5 -m genome_list.tsv 1> ngd.stdout.log 2> ngd.stderr.log


# Rename all files and place in human_readable_long
cd from_NCBI
mkdir human_readable_long
for f in human_readable/refseq/bacteria/*/*/*/*.fna.gz; do
  accn=$(basename $f | cut -d _ -f 1,2);
  long_name=$(echo $f | awk -v var="${accn}" -F'/' '{print $(NF-3)"_"$(NF-2)"_"$(NF-1)"_("var").fna.gz"}');
  ln -sv ../"${f}" human_readable_long/"${long_name}";
done


# gunzip all files
cd human_readable_long
gunzip -f *


conda deactivate
