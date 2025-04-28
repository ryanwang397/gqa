# from genomic_benchmarks.data_check import list_datasets
# from genomic_benchmarks.data_check import info
# from genomic_benchmarks.loc2seq import download_dataset

# info("human_ocr_ensembl", version=0)

# download_dataset("human_ensembl_regulatory", version=0)

from Bio import SeqIO
import csv
import pysam
from Bio.Seq import Seq
import numpy as np
import random
import os
from genomic_benchmarks.loc2seq import download_dataset
import pandas as pd


# def max_difference_from_csv(file_path, column1, column2):
#     # Load the Excel file
#     df = pd.read_csv(file_path)

#     # Calculate the absolute difference between the two columns
#     differences = (df[column1] - df[column2]).abs()

#     # Find the maximum difference
#     max_diff = differences.max()
#     min_diff = differences.min()

#     return max_diff, min_diff


# max_diff, min_diff = max_difference_from_csv(
#     './source_dataset/train/promoter.csv', 'end', 'start')
# print(max_diff, min_diff)
# file_path = 'VQA_RAD.xlsx'
# df = pd.read_excel(file_path, engine='openpyxl')

# unique_counts = {}
# for column in df.columns:
#     unique_counts[column] = df[column].nunique()

# print("unique counts:", unique_counts)

# element_to_count = 'NULL'
# element_counts = 0
# element_counts = df['Q_REPHASE'].value_counts().get(element_to_count, 0)

# print(f"'{element_to_count}' appears:", element_counts)

# ----------------------------------------------------------------------------

# import pandas as pd

# # df = pd.read_parquet('source_dataset/train/0000.parquet')

# # df.to_csv('splice_site.csv', index=False)

# import requests
# import sys
# import re

# chromosome_mapping = {
#     'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5,
#     'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10,
#     'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15,
#     'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20,
#     'chr21': 21, 'chr22': 22, 'X': 'X', 'Y': 'Y'
# }


# def fetch_endpoint(server, request, content_type):
#     """
#     Fetch an endpoint from the server, allow overriding of default content-type
#     """
#     r = requests.get(server+request, headers={"Accept": content_type})

#     if not r.ok:
#         r.raise_for_status()
#         sys.exit()

#     if content_type == 'application/json':
#         return r.json()
#     else:
#         return r.text


# def fetch_blat_result(user_seq, db, output_format='json'):
#     """
#     Fetch BLAT search result from UCSC server.
#     """
#     # UCSC BLAT URL
#     url = "https://genome.ucsc.edu/cgi-bin/hgBlat"

#     params = {
#         'userSeq': user_seq,
#         'type': 'DNA',
#         'db': db,
#         'output': output_format
#     }

#     r = requests.get(url, params=params)

#     if not r.ok:
#         r.raise_for_status()
#         sys.exit()

#     return r.json()


# def preprocess(result):
#     fields = result['fields']
#     seq_info = next(
#         (item for item in result['blat'] if re.match(r'^chr\d+$', item[13])), None)
#     info_dic = {key: value for key, value in zip(fields, seq_info)}
#     # print(info_dic)
#     block_sizes = [int(size) for size in info_dic['blockSizes'].split(',')]
#     q_starts = [int(start) for start in info_dic['qStarts'].split(',')]
#     t_starts = [int(start) for start in info_dic['tStarts'].split(',')]
#     t_starts[0] += 1
#     t_start = info_dic['tStart'] + 1
#     t_starts = [x - t_start for x in t_starts]
#     q_end = info_dic['qEnd'] - 1
#     t_end = info_dic['tEnd'] - t_start
#     t_name = chromosome_mapping[info_dic['tName']]
#     strand = 1 if info_dic['strand'] == '+' else -1
#     return t_name, block_sizes, q_starts, t_starts, q_end, t_end, strand, t_start


# def process_row(row):
#     q_seq = row['seq']
#     # BLAT search query sequence
#     db = 'hg38'
#     result = fetch_blat_result(q_seq, db)
#     # print(result)
#     chromosome, block_sizes, q_starts, t_starts, q_end, t_end, strand, t_start = preprocess(
#         result)

#     # ensembl search reference sequence
#     server = "https://rest.ensembl.org"
#     req = f"/sequence/region/human/{chromosome}:{t_start}..{t_start + t_end}:{strand}?"
#     content_type = "text/plain"
#     t_seq = fetch_endpoint(server, req, content_type)
#     print(t_seq)

#     row['seq'] = t_seq
#     row['chr'] = chromosome
#     row['start'] = t_start
#     row['end'] = t_start + t_end

#     return t_seq


# df = pd.read_csv('source_dataset/dataset_0/enhancer/strong_enhancer.csv')

# df['chr'] = pd.NA
# df['start'] = pd.NA
# df['end'] = pd.NA

# df = df.apply(process_row, axis=1)

# print(df)

# df.to_csv(
#     'source_dataset/dataset_0/enhancer/modified_strong_enhancer.csv', index=False)
# # -----------------------------------------------------------------------------

# ./gfClient -minScore=20 -minIdentity=0 127.0.0.1 1234 . input.fa out.psl

# df = pd.read_csv('source_dataset/dataset_0/enhancer/weak_enhancer.csv')

# with open('input.fa', 'w') as fasta_file:
#     for index, row in df.iterrows():
#         fasta_file.write(f">{index}\n")
#         fasta_file.write(f"{row['seq']}\n")

# # -----------------------------------------------------------------------------

# import csv
# import re


# def read_psl(filename, output_filename):
#     pattern = re.compile(r'chr(\d{1,2}|X|Y)$')
#     with open(filename, 'r') as file, open(output_filename, 'w', newline='') as csvfile:
#         csvwriter = csv.writer(csvfile)
#         csvwriter.writerow(['name', 'chr', 'start', 'end'])
#         for line in file:
#             if line.startswith('psLayout') or line.startswith('match') or line.startswith('-'):
#                 continue
#             fields = line.strip().split('\t')
#             if len(fields) > 17 and (fields[11] == '0' and fields[12] == '200') and fields[0] == '200' and pattern.match(fields[13]):
#                 strand = fields[8]
#                 name = fields[9]
#                 chrom = fields[13]
#                 start = fields[15]
#                 end = fields[16]
#                 csvwriter.writerow([name, chrom, start, end])


# read_psl('source_dataset/dataset_0/enhancer/weak_enhancer.psl',
#          'source_dataset/dataset_0/enhancer/weak_enhancer_processed.csv')

# # -----------------------------------------------------------------------------


# df = pd.read_csv(
#     'source_dataset/dataset_0/enhancer/weak_enhancer_processed.csv')

# is_duplicate_with_next = df['name'] == df['name'].shift(-1)
# is_duplicate_with_prev = df['name'] == df['name'].shift(1)
# df_filtered = df[~(is_duplicate_with_next | is_duplicate_with_prev)]
# # df_filtered = df[df['is_duplicated'] == ]

# # df['Sign'] = 0
# # df.loc[df['Name'] == False, 'Name'] = 1

# # df.drop('is_duplicated', axis=1, inplace=True)

# # print(df)

# df_filtered.to_csv(
#     'source_dataset/dataset_0/enhancer/weak_processed_repeat.csv', index=False)

# -------------------------------------------------------------

# df = pd.read_csv('source_dataset/train/splice_site.csv')
# df_processed = pd.read_csv('source_dataset/train/processed_repeat.csv')
# df['chr'] = ''
# df['start'] = -1
# df['end'] = -1
# df['strand'] = ''
# i = 0
# for index, row in df.iterrows():
#     if row['name'] == df_processed.iloc[i]['Name']:
#         df.loc[index, 'chr'] = df_processed.iloc[i]['Chromosome']
#         df.loc[index, 'start'] = df_processed.iloc[i]['Start']
#         df.loc[index, 'end'] = df_processed.iloc[i]['End']
#         df.loc[index, 'strand'] = df_processed.iloc[i]['Strand']
#         i += 1

# df.to_csv('source_dataset/train/splice_site_processed.csv', index=False)

# ------------------------------------------------------
import pandas as pd

file1 = pd.read_csv(
    'source_dataset/dataset_0/enhancer/weak_processed_repeat.csv')

file2 = pd.read_csv('source_dataset/dataset_0/enhancer/weak_enhancer.csv')

seq_values = []

for index, row in file1.iterrows():
    seq_value = file2.iloc[row['name']]['seq']
    seq_values.append(seq_value)

file1['seq'] = seq_values

file1 = file1.drop(columns=['name'])

file1.to_csv('weak_enhancer_final.csv', index=False)

# ------------------------------------------------------

# from pyensembl import EnsemblRelease

# data = EnsemblRelease(104)

# data.download()
# data.index()

# chromosome_name = "1"
# start = 100000
# end = 100100

# # sequence = data.sequence(region=chromosome_name, start=start, end=end)
# sequence = data.

# from genomic_benchmarks.data_check import list_datasets, info
# print(list_datasets())
# print(info("human_ensembl_regulatory", version=0))
# download_dataset("human_ensembl_regulatory", version=0)


# fasta_file = 'GRCh38_latest_genomic.fna'

# record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# start = 100000 - 1
# end = 100100
# chr1_seq = record_dict[list(record_dict.keys())[0]].seq

# for seq_record in SeqIO.parse("GRCh38_latest_genomic.fna", "fasta"):
#     if seq_record.id.startswith("NC_"):
#         print(seq_record.id)
# NC_000001.11
# NC_000002.12
# NC_000003.12
# NC_000004.12
# NC_000005.10
# NC_000006.12
# NC_000007.14
# NC_000008.11
# NC_000009.12
# NC_000010.11
# NC_000011.10
# NC_000012.12
# NC_000013.11
# NC_000014.9
# NC_000015.10
# NC_000016.10
# NC_000017.11
# NC_000018.10
# NC_000019.10
# NC_000020.11
# NC_000021.9
# NC_000022.11
# NC_000023.11
# NC_000024.10
# NC_012920.1

# segment = chr1_seq[start:end]

# print(segment)


# seq_index = SeqIO.index(fasta_file, "fasta")

# chromosome_id = "NC_000011.10"

# seq_record = seq_index[chromosome_id]

# start, end = 86708800, 86709200

# subseq = seq_record.seq[start:end]
# print(subseq)

# seq_index.close()

# # -----------------------------------------------------------------------------


# csv_path = 'source_dataset/train/ocr.csv'
# df = pd.read_csv(csv_path)

# filtered_df = df[(df['end'] - df['start'] >= 200) &
#                  (df['end'] - df['start'] <= 500)]

# sample_df = filtered_df.sample(n=min(10000, len(filtered_df)), random_state=1)
# diff = df['end'] - df['start']
# print(diff.value_counts())

# txt_folder_path = 'genomic_benchmarks/human_ensembl_regulatory/train/ocr'


# def read_seq_from_txt(id_value):
#     try:
#         with open(os.path.join(txt_folder_path, f'{id_value}.txt'), 'r') as file:
#             return file.read().strip()
#     except FileNotFoundError:
#         return ''


# sample_df['seq'] = sample_df['id'].apply(read_seq_from_txt)

# sample_df = sample_df.drop('strand', axis=1)

# new_csv_path = 'used_data/ocr.csv'
# sample_df.to_csv(new_csv_path, index=False)

# # -----------------------------------------------------------------------------


# NC_000001.11
# NC_000002.12
# NC_000003.12
# NC_000004.12
# NC_000005.10
# NC_000006.12
# NC_000007.14
# NC_000008.11
# NC_000009.12
# NC_000010.11
# NC_000011.10
# NC_000012.12
# NC_000013.11
# NC_000014.9
# NC_000015.10
# NC_000016.10
# NC_000017.11
# NC_000018.10
# NC_000019.10
# NC_000020.11
# NC_000021.9
# NC_000022.11
# NC_000023.11
# NC_000024.10

# chromosome_mapping = {
#     'chr1': 'NC_000001.11', 'chr2': 'NC_000002.12', 'chr3': 'NC_000003.12', 'chr4': 'NC_000004.12', 'chr5': 'NC_000005.10',
#     'chr6': 'NC_000006.12', 'chr7': 'NC_000007.14', 'chr8': 'NC_000008.11', 'chr9': 'NC_000009.12', 'chr10': 'NC_000010.11',
#     'chr11': 'NC_000011.10', 'chr12': 'NC_000012.12', 'chr13': 'NC_000013.11', 'chr14': 'NC_000014.9', 'chr15': 'NC_000015.10',
#     'chr16': 'NC_000016.10', 'chr17': 'NC_000017.11', 'chr18': 'NC_000018.10', 'chr19': 'NC_000019.10', 'chr20': 'NC_000020.11',
#     'chr21': 'NC_000021.9', 'chr22': 'NC_000022.11', 'chrX': 'NC_000023.11', 'chrY': 'NC_000024.10'
# }

# df = pd.read_csv('source_dataset/train/splice_site_processed.csv')

# fasta_file = 'GRCh38_latest_genomic.fna'
# seq_index = SeqIO.index(fasta_file, "fasta")


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
#     return new_start, new_end, str(new_seq)


# for index, row in df.iterrows():
#     if index == 10000:
#         break
#     if df.at[index, 'start'] == -1:
#         continue
#     new_start, new_end, new_seq = adjust_sequence(
#         row['sequence'], row['start'], row['end'], row['chr'], random.randint(200, 500), seq_index)
#     df.at[index, 'start'] = new_start
#     df.at[index, 'end'] = new_end
#     df.at[index, 'sequence'] = new_seq
#     print(index)

# df = df.drop('strand', axis=1)

# df.to_csv('used_data/test.csv', index=False)

# # -----------------------------------------------------------------------------

# questions = ['Is this an enhancer sequence?', 'Could this sequence function as an enhancer?', 'Is it possible that this DNA sequence acts as an enhancer?',
#              'Does this sequence serve the role of an enhancer in the genome?',
#              'Can we identify this sequence as an enhancer?',
#              'Is this sequence characterized as an enhancer?',
#              'Might this DNA sequence be classified as an enhancer?',
#              'Does this particular sequence have the properties of an enhancer?',
#              'Would this sequence be considered an enhancer?',
#              'Does this sequence exhibit the typical features of an enhancer?',
#              ]

# questions = ['Is this a promoter sequence?',
#              'Could this sequence be acting as a promoter?',
#              'Does this sequence have the characteristics of a promoter?',
#              'Can this sequence be identified as a promoter in genetic transcription?',
#              'Might this DNA sequence function in the capacity of a promoter?',
#              'Is it possible for this DNA sequence to serve as a promoter?',
#              'Do the features of this sequence align with those of known promoters?',
#              'Is this sequence functioning in the role of a promoter?',
#              'Is this sequence responsible for starting the transcription process?',
#              'Is this sequence located upstream of a gene?']

questions = ['Is this a DNA sequence belonging to open chromatin regions?',
             'Does this sequence reside within regions of open chromatin?',
             'Does this DNA fragment correspond to an open chromatin area?',
             'Is this sequence a component of open chromatin regions?',
             'Is it possible that this sequence is situated within an open chromatin region?',
             'Does this sequence represent a section of the open chromatin region?',
             'Does this particular sequence align with open chromatin regions?',
             'Would this sequence be considered part of the open chromatin regions?',
             'Does the location of this sequence correspond to open chromatin areas?',
             'Is this an ocr sequence?']


# df = pd.read_csv('used_data/ocr.csv')

# question_column = np.repeat(questions, 1000)

# answer_column = np.repeat("yes", 10000)

# df['question'] = question_column
# df['answer'] = answer_column

# # df.to_csv('used_data/dataset.csv', index=False)
# with open('used_data/dataset.csv', 'a') as f:
#     df.to_csv(f, header=False, index=False)


# ocr_df = pd.read_csv('used_data/promoter.csv').iloc[:5000]
# promoter_df = pd.read_csv('used_data/enhancer.csv').iloc[:5000]

# combined_df = pd.concat([ocr_df, promoter_df], ignore_index=True)

# shuffled_df = combined_df.sample(frac=1).reset_index(drop=True)

# shuffled_df['question'] = np.repeat(questions, 1000)
# shuffled_df['answer'] = "no"

# with open('used_data/dataset.csv', 'a') as f:
#     shuffled_df.to_csv(f, header=False, index=False)


# df = pd.read_csv('used_data/test.csv')

# df.rename(columns={'sequence': 'seq'}, inplace=True)

# df = df.head(10000)

# df = df.dropna(subset=['seq'])

# df = df[df['start'] != -1]

# df['seq'] = df['seq'].str.upper()

# df.to_csv('used_data/splice_site.csv', index=False)


# df = pd.read_csv('used_data/splice_site.csv')

# df_acceptor = df[df['label'] == 0]
# df_donor = df[df['label'] == 1]

# df_acceptor['site'] = 'acceptor'
# df_donor['site'] = 'donor'

# df_acceptor.to_csv('used_data/acceptor.csv', index=False)
# df_donor.to_csv('used_data/donor.csv', index=False)

# -------------------------------------------------------------------------------


# csv_file = 'used_data/dataset.csv'
# df = pd.read_csv(csv_file)

# def modify_sequence(chrom, start, end, seq, vcf_file):
#     tabix_file = pysam.TabixFile(vcf_file)
#     try:
#         records = list(tabix_file.fetch(chrom, start + 1, end))
#     except ValueError:
#         return seq, 'N'

#     if not records:
#         return seq, 'N'

# #     random_record = random.choice(records)
# #     fields = random_record.strip().split('\t')
# #     var_pos, ref, alts = int(fields[1]), fields[3], fields[4].split(',')

# #     if len(ref) == 1 and all(len(alt) == 1 for alt in alts):
# #         seq_list = list(seq)
# #         random_alt = random.choice(alts)
# #         seq_list[var_pos - start - 1] = random_alt
# #         return ''.join(seq_list), ref + '->' + random_alt
# #     else:
# #         return seq, 'N'
#     # print(records[0])
# #     fields = records[0].strip().split('\t')
# #     fields1 = records[1].strip().split('\t')
# #     print(fields)


# for index, row in df.iterrows():
#     modified_seq, variant = modify_sequence(
#         chrom=str(row['region']),
#         start=int(row['start']),
#         end=int(row['end']),
#         seq=row['seq'],
#         vcf_file='source_dataset/variants/ALL.{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'.format(
#             row['region'])
#     )
#     df.at[index, 'seq'] = modified_seq
#     df.at[index, 'variant'] = variant

# df.to_csv('modified_sequences.csv', index=False)

# # -------------------------------------------------------------------------------


# csv_file_path = 'add_variants/modified_sequences.csv'

# df = pd.read_csv(csv_file_path)

# variant_column = df['answer1']

# unique_variants = variant_column.unique()
# number_of_unique_variants = len(unique_variants)

# print("Total number of unique variants in 'variant' column:",
#       number_of_unique_variants)

# print("Unique variants:", unique_variants)

# variant_counts = variant_column.value_counts()

# print("Counts of each unique variant:", variant_counts)
# # -------------------------------------------------------------------------------

# csv_file_path = 'add_variants/modified_sequences.csv'

# df = pd.read_csv(csv_file_path)

# df['question1'] = "what's the variant in this sequence?"

# df = df.rename(columns={'variant': 'answer1'})

# df.to_csv(csv_file_path, index=False)

# ----------------------------------------------------------------------------------
# csv_file = 'used_data/dataset.csv'
# df = pd.read_csv(csv_file)

# def search_sequence(chrom, start, end, seq, vcf_file):
#     tabix_file = pysam.TabixFile(vcf_file)
#     try:
#         records = list(tabix_file.fetch(chrom, start + 1, end))
#     except ValueError:
#         return seq, 'N'

#     if not records:
#         return seq, 'N'

#     person_variant = {}
#     for record in records:
#         fields = record.strip().split('\t')
#         for index, field in enumerate(fields):
#             if field == '1|0' or field == '1|1' or field == '0|1':
#                 person_variant[index] = person_variant.get(index, 0) + 1
#     cnt = len(person_variant)
#     num_per_person = sum(person_variant.values()) / cnt if cnt > 0 else 0
#     return cnt, num_per_person


# for index, row in df.iterrows():
#     variant_person_num, variant_num_per_person = search_sequence(
#         chrom=str(row['region']),
#         start=int(row['start']),
#         end=int(row['end']),
#         seq=row['seq'],
#         vcf_file='source_dataset/variants/ALL.{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'.format(
#             row['region'])
#     )
#     df.at[index, 'variant_person_num'] = variant_person_num
#     df.at[index, 'variant_num_per_person'] = variant_num_per_person

# df.to_csv('cnt_sequences.csv', index=False)

# ----------------------------------------------------------------------------------


# def read_chromosome_positions_from_tsv(file_path):
#     chromosome_positions = []
#     with open(file_path, 'r', newline='') as file:
#         reader = csv.DictReader(file, delimiter='\t')
#         for row in reader:
#             chrom_values = row['CHR_ID'].split(';')
#             pos_values = row['CHR_POS'].split(';')
#             for chrom, pos in zip(chrom_values, pos_values):
#                 if chrom and pos:  # Check if both chrom and pos are not empty
#                     try:
#                         pos = int(pos)
#                     except ValueError:
#                         # If pos cannot be converted to int, skip this iteration
#                         continue
#                     chromosome_positions.append((chrom, int(pos)))
#     return chromosome_positions


# def read_chromosome_intervals_from_csv(file_path):
#     chromosome_intervals = {}
#     with open(file_path, 'r') as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             chrom = row['region']
#             if chrom.startswith('chr'):
#                 chrom = chrom[3:]
#             start = int(row['start'])
#             end = int(row['end'])
#             if chrom not in chromosome_intervals:
#                 chromosome_intervals[chrom] = []
#             chromosome_intervals[chrom].append((start, end))
#     return chromosome_intervals


# def count_positions_in_intervals(chromosome_intervals, chromosome_positions):
#     count = 0
#     for chrom, pos in chromosome_positions:
#         if chrom in chromosome_intervals:
#             for start, end in chromosome_intervals[chrom]:
#                 if start <= pos <= end:
#                     count += 1
#                     break  # If position is found in one interval, break the inner loop
#     return count


# # Example data file paths
# positions_file = 'disease/gwas-association-downloaded_2024-05-04-EFO_0000216-withChildTraits.tsv'
# intervals_file = 'cnt_sequences.csv'

# # Read data from CSV files
# chromosome_positions = read_chromosome_positions_from_tsv(positions_file)
# chromosome_intervals = read_chromosome_intervals_from_csv(intervals_file)

# # Count positions in intervals
# result = count_positions_in_intervals(
#     chromosome_intervals, chromosome_positions)
# print("Number of positions in intervals:", result)

# ----------------------------------------------------------------------------------
# tsv_file_path = 'source_dataset/disease/gwas-association-downloaded_2024-05-04-EFO_0000216-withChildTraits.tsv'

# # def fetch_variants(chr_id, pos, vcf_path):
# #     tabix_file = pysam.TabixFile(vcf_path)
# #     try:
# #         records = list(tabix_file.fetch(chr_id, pos - 1, pos))
# #         for record in records:
# #             fields = record.strip().split('\t')
# #             print(fields)
# #         return len(records) > 0
# #     except ValueError:
# #         print('wrong')

# def fetch_variants(chr_id, pos, vcf_file_path):
#     try:
#         with pysam.VariantFile(vcf_file_path) as vcf:
#             chr = 'chr' + chr_id
#             records = list(vcf.fetch(chr, pos - 1, pos))
#             # for record in records:
#             #     if record.pos == pos:
#             #         print(f"Variant found at {chr_id}:{pos}")
#             #         print(
#             #             f"ID: {record.id}, REF: {record.ref}, ALT: {record.alts}")
#             #     else:
#             #         print(
#             #             f"No exact match found at {chr_id}:{pos}, but nearby variant exists.")
#             return len(records) > 0
#     except FileNotFoundError:
#         print(f"The file {vcf_file_path} does not exist.")
#     except ValueError:
#         print(
#             f"Error accessing variants for {chr_id}:{pos} in {vcf_file_path}.")


# cnt = 0
# with open(tsv_file_path, 'r') as tsvfile:
#     reader = csv.DictReader(tsvfile, delimiter='\t')
#     for i, row in enumerate(reader):
#         # if i == 1:
#         #     break
#         chr_id = row['CHR_ID']
#         if chr_id == '':
#             continue
#         try:
#             pos = int(row['CHR_POS'])
#         except ValueError:
#             print(
#                 f"Skipping invalid position: {row['CHR_POS']} for chromosome {chr_id}")
#             continue
#         vcf_file_path = f'source_dataset/variants/ALL.chr{chr_id}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'

#         if fetch_variants(chr_id, pos, vcf_file_path):
#             print(
#                 f"Variant found at chr{chr_id}:{pos} in file {vcf_file_path}")
#             cnt += 1

# print(cnt)

# ----------------------------------------------------------------------------------

# def fetch_variants(chr_id, pos, vcf_file_path):
#     try:
#         with pysam.VariantFile(vcf_file_path) as vcf:
#             chr = 'chr' + chr_id
#             records = list(vcf.fetch(chr, pos - 1, pos))
#             # for record in records:
#             #     if record.pos == pos:
#             #         print(f"Variant found at {chr_id}:{pos}")
#             #         print(
#             #             f"ID: {record.id}, REF: {record.ref}, ALT: {record.alts}")
#             #     else:
#             #         print(
#             #             f"No exact match found at {chr_id}:{pos}, but nearby variant exists.")
#             return len(records) > 0
#     except FileNotFoundError:
#         print(f"The file {vcf_file_path} does not exist.")
#     except ValueError:
#         print(
#             f"Error accessing variants for {chr_id}:{pos} in {vcf_file_path}.")


# cnt = 0
# rows_to_save = []

# csv_file_path = 'cancer.csv'
# with open(csv_file_path, 'r') as csvfile:
#     reader = csv.DictReader(csvfile, delimiter=',')
#     for row in reader:
#         chr_id = row['CHR_ID']
#         if chr_id == '':
#             continue
#         try:
#             pos = int(row['CHR_POS'])
#         except ValueError:
#             print(
#                 f"Skipping invalid position: {row['CHR_POS']} for chromosome {chr_id}")
#             continue

#         vcf_file_path = f'source_dataset/variants/ALL.chr{chr_id}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'

#         if fetch_variants(chr_id, pos, vcf_file_path):
#             print(
#                 f"Variant found at chr{chr_id}:{pos} in file {vcf_file_path}")
#             cnt += 1
#             rows_to_save.append(row)

# output_file_path = 'filtered_cancer.csv'
# with open(output_file_path, 'w', newline='') as f:
#     writer = csv.DictWriter(f, fieldnames=reader.fieldnames)
#     writer.writeheader()
#     writer.writerows(rows_to_save)

# print(f"Total variants found: {cnt}")
