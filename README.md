# Lepto in silico PCR

## Requirements

Create conda environment and install packages:

```
conda create -n genomedl
conda activate genomedl
mamba install ncbi-genome-download ispcr
```

## Usage

Bash files were ran on UGE using qsub command:

1. Download all 769 refseq files, make filenames readable, and ungzip.
```
qsub fetch_genomes.sh
```

2. Create new folder in same directory as the from_NCBI directory. 

```
mkdir primers
```

3. Create file primers.tab with primer sequences and place into primers directory. Make sure spaces are tabs.
```
echo "LeptoPrimer    GAGTAACACGTGGGTAATCTTCCT        TTTACCCCACCAACTAGCTAATC" > primers.tab
```

4. Find primers from genomes.

```
qsub find_primers.sh
```

5. Concatenate forward and reverse primers by using concat_files function in python script.


6. Run mafft on concatenated files.

```
qsub mafft.sh
```

7. Analyze genome and save to file:

```
python3 analyze_lepto_genome.py
```

8. Check lepto_in_silico.xlsx for finished analysis.


## Notes

The analyze_lepto_genome.py script was written where speed was key. Time was taken afterwards to add comments, show some data, and fix a few errors in the analyze_lepto_genomes.ipynb. 

The word 'qsub' is used to send scripts to the HPC (UGE). This can be substituted with 'bash' to run locally.