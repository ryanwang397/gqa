from datasets import load_dataset
import ast
from scipy.stats import norm
import re
from collections import Counter, defaultdict
import glob
import ast
import random
import numpy as np
import pandas as pd

# ------------------------------------------------------------------------

# import pandas as pd
# import random
# from Bio import SeqIO

# chromosome_mapping = {
#     'chr1': 'NC_000001.11', 'chr2': 'NC_000002.12', 'chr3': 'NC_000003.12', 'chr4': 'NC_000004.12', 'chr5': 'NC_000005.10',
#     'chr6': 'NC_000006.12', 'chr7': 'NC_000007.14', 'chr8': 'NC_000008.11', 'chr9': 'NC_000009.12', 'chr10': 'NC_000010.11',
#     'chr11': 'NC_000011.10', 'chr12': 'NC_000012.12', 'chr13': 'NC_000013.11', 'chr14': 'NC_000014.9', 'chr15': 'NC_000015.10',
#     'chr16': 'NC_000016.10', 'chr17': 'NC_000017.11', 'chr18': 'NC_000018.10', 'chr19': 'NC_000019.10', 'chr20': 'NC_000020.11',
#     'chr21': 'NC_000021.9', 'chr22': 'NC_000022.11', 'chrX': 'NC_000023.11', 'chrY': 'NC_000024.10'
# }


# def extract_and_mutate_sequence(chr_seq, chrom, pos, alt):
#     # 1-based to 0-based conversion for position
#     pos = pos - 1
#     pre_length = random.randint(200, 250)
#     post_length = random.randint(200, 250)
#     start = max(0, pos - pre_length)
#     end = pos + post_length + 1  # Adjust end index to include the last base
#     sequence = chr_seq.seq[start:end].upper()

#     # Extract reference allele at position (0-based index in substring)
#     ref_allele = sequence[pre_length]

#     # Generate mutated sequence
#     # mutated_sequence = sequence[:pre_length] + \
#     #     alt + sequence[pre_length+1:]

#     # Convert indices back to 1-based for output
#     return str(sequence), [pos + 1, ref_allele, alt], start + 1, end


# def process_snps(csv_file, seq_index):
#     df = pd.read_csv(csv_file)
#     mutation_info_list = []
#     seq_list = []
#     start_list = []
#     end_list = []

#     for index, row in df.iterrows():
#         if index == 1000:
#             break
#         chrom = 'chr' + str(row['CHR_ID'])
#         pos = int(row['CHR_POS'])
#         risk_allele_info = row['STRONGEST SNP-RISK ALLELE'].split('-')
#         rsid = risk_allele_info[0]
#         alt_allele = risk_allele_info[1]

#         # Extract sequence and apply mutation
#         mutated_seq, mutation_info, start_pos, end_pos = extract_and_mutate_sequence(
#             seq_index[chromosome_mapping[chrom]], chrom, pos, alt_allele)
#         mutation_info_list.append(mutation_info)
#         seq_list.append(mutated_seq)
#         start_list.append(start_pos)
#         end_list.append(end_pos)

#     # Add new columns to DataFrame
#     df = df.head(1000)
#     print(df.shape)
#     df['seq'] = seq_list
#     df['start'] = start_list
#     df['end'] = end_list
#     df['variants'] = mutation_info_list

#     return df


# # Using the example
# csv_path = 'source_dataset/disease/filtered_cancer.csv'
# fasta_file = 'GRCh38_latest_genomic.fna'
# seq_index = SeqIO.index(fasta_file, "fasta")
# processed_df = process_snps(csv_path, seq_index)

# # Save the processed DataFrame to a new CSV file
# processed_df.to_csv(
#     'source_dataset/dataset_0/disease/cancer_temp.csv', index=False)


# def adjust_sequence(seq, start, end, chr, target_length, ref_seqs):
#     current_length = len(seq)
#     diff = target_length - current_length
#     if diff < 0:
#         cut_each_side = abs(diff) // 2
#         new_seq = seq[cut_each_side:-cut_each_side]
#         new_start = start + cut_each_side
#         new_end = end - cut_each_side
#     else:
#         extend_each_side = diff // 2
#         chr_seq = seq_index[chromosome_mapping[chr]]
#         new_start = max(0, start - extend_each_side)
#         new_end = end + extend_each_side
#         new_seq = chr_seq.seq[new_start:new_end]

# ------------------------------------------------------------------------
# import ast

# csv_file = 'source_dataset/dataset_0/disease/cancer_temp.csv'
# df = pd.read_csv(csv_file)


# def modify_sequence(idx, chrom, start, end, seq, vcf_file, posi, ori, vt):
#     tabix_file = pysam.TabixFile(vcf_file)
#     try:
#         records = list(tabix_file.fetch(chrom, start, end))
#     except ValueError:
#         return seq, 'N'

#     if not records:
#         return seq, 'N'

#     person_variant = {}
#     variants = []
#     modified_seq = ''
#     global cnt
#     for record in records:
#         fields = record.strip().split('\t')
#         for index, field in enumerate(fields):
#             if field == '1|0' or field == '1|1':
#                 if (index, 0) not in person_variant:
#                     person_variant[(index, 0)] = [
#                         [int(fields[1]), fields[3], fields[4].split(',')]]
#                 else:
#                     person_variant[(index, 0)].append(
#                         [int(fields[1]), fields[3], fields[4].split(',')])
#             if field == '1|1' or field == '0|1':
#                 if (index, 1) not in person_variant:
#                     person_variant[(index, 1)] = [
#                         [int(fields[1]), fields[3], fields[4].split(',')]]
#                 else:
#                     person_variant[(index, 1)].append(
#                         [int(fields[1]), fields[3], fields[4].split(',')])
#     # flag = random.random()
#     flag = 0
#     for (index, type), tmp_variants in person_variant.items():
#         if flag == 1:
#             break
#         for sublist in tmp_variants:
#             if sublist[0] == posi and sublist[1] == ori and sublist[2][0] == vt:
#                 variants = tmp_variants
#                 flag = 1
#                 break
#     seq_list = list(seq)
#     if len(variants) != 0 and all(len(variant[2]) == 1 for variant in variants):
#         for variant in variants:
#             pos = variant[0]
#             original_base = variant[1]
#             mutated_base = variant[2][0]

#             if seq_list[pos - start] == original_base:
#                 seq_list[pos - start] = mutated_base
#             else:
#                 return seq, 'N'
#     else:
#         return seq, 'N'

#     return ''.join(seq_list), json.dumps(variants)


# chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
#        'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

# for index, row in df.iterrows():
#     v_list = ast.literal_eval(row['variants'])
#     if 'chr' + str(row['CHR_ID']) in chr:
#         modified_seq, variants = modify_sequence(
#             idx=index,
#             chrom='chr' + str(row['CHR_ID']),
#             start=int(row['start']),
#             end=int(row['end']),
#             seq=row['seq'],
#             vcf_file='source_dataset/variants/ALL.{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'.format(
#                 'chr' + str(row['CHR_ID'])),
#             posi=int(v_list[0]),
#             ori=v_list[1],
#             vt=v_list[2]
#         )
#     else:
#         modified_seq, variants = row['seq'], 'N'
#     df.at[index, 'seq'] = modified_seq
#     df.at[index, 'variants'] = variants

# # 保存修改后的CSV
# df.to_csv(
#     'source_dataset/dataset_0/disease/cancer_with_variants.csv', index=False)

# #---------------------------------------------------------------------

df = pd.read_csv('disease/cancer_with_variants.csv')

df['cancer_related'] = df['variants'].apply(
    lambda x: 'yes' if x != 'N' else 'no')

df['cancer_kind'] = df.apply(lambda row: 'N' if row['variants'] == 'N' else row['MAPPED_TRAIT'].split(
    ',')[0] if pd.notnull(row['MAPPED_TRAIT']) else '', axis=1)


def get_variant_for_chr_pos(variants_str, chr_pos):
    try:
        variants = ast.literal_eval(variants_str)
        for variant in variants:
            if variant[0] == chr_pos:
                return f"{variant[1]}->{variant[2][0]}"
        return 'N'
    except:
        return 'N'


df['cancer_related_variant'] = df.apply(
    lambda row: get_variant_for_chr_pos(row['variants'], row['CHR_POS']), axis=1)

df.to_csv('disease/cancer_with_variants.csv', index=False)
