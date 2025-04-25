import ast
from Bio import SeqIO
from genomic_benchmarks.loc2seq import download_dataset
import os
import random
from Bio.Seq import Seq
import pysam
import csv
import numpy as np
import json
import pandas as pd

csv_file = 'source_dataset/dataset_0/promoter/notata_promoter.csv'
df = pd.read_csv(csv_file)
cnt = 0


def modify_sequence(idx, chrom, start, end, seq, vcf_file):
    tabix_file = pysam.TabixFile(vcf_file)
    try:
        records = list(tabix_file.fetch(chrom, start, end))
    except ValueError:
        return seq, 'N'

    if not records:
        return seq, 'N'

    person_variant = {}
    variants = []
    modified_seq = ''
    global cnt
    for record in records:
        fields = record.strip().split('\t')
        for index, field in enumerate(fields):
            if field == '1|0' or field == '1|1':
                if (index, 0) not in person_variant:
                    person_variant[(index, 0)] = [
                        [int(fields[1]), fields[3], fields[4].split(',')]]
                else:
                    person_variant[(index, 0)].append(
                        [int(fields[1]), fields[3], fields[4].split(',')])
            if field == '1|1' or field == '0|1':
                if (index, 1) not in person_variant:
                    person_variant[(index, 1)] = [
                        [int(fields[1]), fields[3], fields[4].split(',')]]
                else:
                    person_variant[(index, 1)].append(
                        [int(fields[1]), fields[3], fields[4].split(',')])
    flag = random.random()
    for (index, type), tmp_variants in person_variant.items():
        if flag < 0.5 and len(tmp_variants) == 1:
            variants = tmp_variants
            break
        if flag >= 0.5 and len(tmp_variants) > 1:
            variants = tmp_variants
            break
    seq_list = list(seq)
    if len(variants) != 0 and all(len(variant[2]) == 1 for variant in variants):
        for variant in variants:
            pos = variant[0]
            original_base = variant[1]
            mutated_base = variant[2][0]

            if seq_list[pos - start] == original_base:
                seq_list[pos - start] = mutated_base
            else:
                cnt += 1
                return seq, 'N'
    else:
        return seq, 'N'

    return ''.join(seq_list), json.dumps(variants)


chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

for index, row in df.iterrows():
    if row['chr'] in chr:
        modified_seq, variants = modify_sequence(
            idx=index,
            chrom=str(row['chr']),
            start=int(row['start']),
            end=int(row['end']),
            seq=row['seq'],
            vcf_file='source_dataset/variants/ALL.{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'.format(
                row['chr'])
        )
    else:
        modified_seq, variants = row['seq'], 'N'
    df.at[index, 'seq'] = modified_seq
    df.at[index, 'variants'] = variants

print(cnt)
df.to_csv(
    'source_dataset/dataset_0/promoter/notata_promoter_with_variants.csv', index=False)
