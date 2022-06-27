#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import fnmatch
import re
from itertools import groupby
from Bio import SeqIO


def concat_files(path1, path2, split_char, ext):
    # path1 = current file path
    # path2 = where the new file will go
    # split_char = what to split by, ex. '_'
    # ext = file extension, ex = ".fasta"
    
    files = fnmatch.filter(os.listdir(path1), "*" + ext)
    def split_name(fname): return fname.split(split_char)[1]
    for name, folder in groupby(sorted(files, key=split_name), key=split_name):
        new_name = str(name) + "-merged" + ext
        with open(os.path.join(path2, new_name), 'w') as outfile:
            for fname in folder:
                with open(os.path.join(path1, fname), 'r') as infile:
                    outfile.write(infile.read())


def alignment(path_name, which_format):
    # path_name = path to files
    # which_format = merged or individual
    names = []
    alignments = []
    sequence_num = []
    accessions = []
    lengths = []
    name = ''
    records = []
    for file in fnmatch.filter(os.listdir(path_name), "*" + '.txt'):
        alignment = 0
        matches = 0
        mismatches = 0
        if which_format == 'serovar':
            fasta_name_list = file.split('_')[:-2]
            fasta_name_list = [name + "_" for name in fasta_name_list]
            fasta_name = ''.join(map(str, fasta_name_list))
            records = list(SeqIO.parse(os.path.join(path_name, fasta_name + 'concat.fna'), "fasta"))
        elif which_format == 'merged':
            name = file.split('-')[0]
            fasta_name = file.split('_')[0]
            records = list(SeqIO.parse(os.path.join(path_name,fasta_name), "fasta"))
        elif which_format == 'individual':
            fasta_name_list = file.split('_')[:-1]
            
            acc = re.search(r'.+\((GCF_.+)\)_.+', file)
            accessions.append(acc.group(1))
            
            fasta_name_list = [name + "_" for name in fasta_name_list]
            fasta_name = ''.join(map(str, fasta_name_list))
            records = list(SeqIO.parse(os.path.join(path_name, fasta_name + 'concat.fna'), "fasta"))
        elif which_format == 'bacterio':
            fasta_name = file.split('_')[0]
            # for sequences
            records = list(SeqIO.parse(os.path.join(path_name, fasta_name + '_concat.fa'), "fasta"))
            #for name
            get_name = list(SeqIO.parse(os.path.join(path_name, fasta_name + '.fna'), 'fasta'))
            for record in get_name:
                split_id = record.id
                split_id = split_id.split('_')
                name = split_id[1]
                accessions.append(split_id[3])
        else:
            break

        sequence_num.append(len(records))

        if which_format != 'serovar':
            tmp_lengths = []
            for record in records:
                tmp_lengths.append(len(record.seq))

            lengths.append(tmp_lengths)

        with open(os.path.join(path_name, file)) as mafft_file:
            next(mafft_file)
            for line in mafft_file:
                if line.startswith("NZ_"):
                    continue
                else:
                    matches += line.count("*")
                    mismatches += line.count(".")
        if matches != 0:
            alignment = (matches/(matches+mismatches))*100
        else:
            alignment = 0
        names.append(name)
        alignments.append(alignment)
    return names, alignments, sequence_num, accessions, lengths


path = "/scicomp/home-pure/top8/lepto/lepto_genomes/primers"
merge_path = "/scicomp/home-pure/top8/lepto/lepto_genomes/primers/merged_primers"
large_serovar_path = "/scicomp/home-pure/top8/lepto/lepto_genomes/primers/large_serovar_concat/"
#ext = ".fasta"

#concat_files(path, merge_path, '_', '.fasta')

merge_names, merge_alignments, merge_sequence_num, merge_acc, merge_lengths = alignment(merge_path, 'merged')
ind_names, ind_alignments, ind_sequence_num, ind_acc, ind_lengths = alignment(path, 'individual')
serovar_names, serovar_alignments, serovar_sequence_num, serovar_acc, serovar_lengths = alignment(large_serovar_path, 'serovar')

#merged data
merged_data = {'species':merge_names, '# of seq':merge_sequence_num, 'nuc. identity':merge_alignments}
merged = pd.DataFrame(merged_data)

# large serovar data
large_serovar_data = {'# of seq':serovar_sequence_num, 'nuc. identity':serovar_alignments}
large_serovar = pd.DataFrame(large_serovar_data)
print("large servoar folder\n", large_serovar.head(), "\n\n")

#percent length per AE010300 (153 bp)
percent_lengths = []
for i in ind_lengths:
    tmp_perc_lengths = []
    for j in i:
        tmp_perc_lengths.append(round((j/153)*100, 2))
        
    percent_lengths.append(tmp_perc_lengths)

#ind data
ind_data = {'accession':ind_acc, '# of seq':ind_sequence_num, 'nuc. identity':ind_alignments, 'length per AE010300':percent_lengths}
ind = pd.DataFrame(ind_data)


# genome_list.tsv
df = pd.read_csv('/scicomp/home-pure/top8/lepto/lepto_genomes/genome_list.tsv', sep="\t")
df = df[['assembly_accession', 'organism_name', 'infraspecific_name']]
df.rename({'assembly_accession':'accession', 'organism_name':'species_org', 'infraspecific_name':'strain'}, axis=1, inplace=True)
df['species_org'] = df['species_org'].str.replace("Leptospira", "")
df['strain'] = df['strain'].str.replace('strain=', '')
df['serovar'] = df['species_org'].str.contains('serovar')
df['serovar'] = df['serovar'].map({True: 'serovar', False: ''})
move_serovar = df.pop('serovar')
df.insert(2, 'serovar', move_serovar)
df['species'] = df['species_org'].str.split(' ').str[1]
move_species = df.pop('species')
df.insert(1, 'species', move_species)
remove_species_org = df.pop('species_org')


# merge genome_list and each acc data
df2 = pd.merge(df, ind)

# remove brackets from percent length per AE10300 column
for c in df2.columns:
    df2[c]= df2.apply(lambda x: str(x[c]).replace("[",'').replace(']',''), axis=1)
df2['length per AE010300'] = df2['length per AE010300'].replace('', 0)


df2.to_csv('genome_by_accession.tsv', sep="\t", index=False)
merged.to_csv('genome_by_species.tsv', sep="\t", index=False)

zeros = 0
less_than_5 = 0 
greater_than_10 = 0
for i in merge_sequence_num:
    if i == 0:
        zeros += 1
    elif 0 < i >= 5:
        less_than_5 += 1
    else:
        greater_than_10 += 1

print('zeros: ', zeros, ' --- <5: ', less_than_5, ' --- >10: ', greater_than_10)

# Obtain accessions where row == 'serovar'
serovar_df = df2.copy()
rows_to_drop = serovar_df[serovar_df['serovar'] != 'serovar'].index
serovar_df.drop(rows_to_drop, inplace=True)
serovar_acc = serovar_df['accession']
serovar_acc.to_csv('serovar_accession.tsv', sep="\t", index=False)

###################
#while IFS= read -r file; do cat *\("$file"\)_concat.fna >> serovar_concat.fna; done < "../serovar_accession.tsv"
#############################

approved_species = pd.read_csv('/scicomp/home-pure/top8/lepto/lepto_genomes/approved_species.csv', sep=",")
approved_species = approved_species[['Name', 'Clade']]
approved_species.rename({'Clade':'clade'})
approved_species.to_csv('check_names.tsv', sep='\t')
approved_species['species'] = approved_species['Name'].astype(str).str.split().str[1]
approved_species.drop('Name', axis=1, inplace=True)


bacterio_path = "/scicomp/home-pure/top8/lepto/lepto_genomes/primers/bacterio_serovars"
bacterio_names, bacterio_alignments, bacterio_sequence_num, bacterio_acc, bacterio_lengths = alignment(bacterio_path, 'bacterio')

#percent length per AE010300 (153 bp)
bacterio_perc_lengths = []
for i in bacterio_lengths:
    tmp_perc_lengths = []
    for j in i:
        tmp_perc_lengths.append(round((j/153)*100, 2))

    bacterio_perc_lengths.append(tmp_perc_lengths)
    

bacterio_data = {"species":bacterio_names, "# of seq":bacterio_sequence_num, "nuc. identity":bacterio_alignments, "length per AE10300":bacterio_perc_lengths}
bacterio_df = pd.DataFrame(bacterio_data)
bacterio_merge = pd.merge(bacterio_df, approved_species, how="outer")
bacterio_merge = bacterio_merge.groupby("species").first().reset_index()

for c in bacterio_merge.columns:
    bacterio_merge[c]= bacterio_merge.apply(lambda x: str(x[c]).replace("[",'').replace(']',''), axis=1)

bacterio_merge['length per AE10300'] = bacterio_merge['length per AE10300'].str.replace(' ', '0')

bacterio_merge.to_csv('approved_serovars.tsv', sep="\t", index=False)
