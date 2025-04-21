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
bs_ds = 'regulatory/acceptor_with_variants.csv'
ot_ds = ['regulatory/ocr_with_variants.csv', 'regulatory/donor_with_variants.csv',
         'regulatory/promoter_with_variants.csv', 'regulatory/enhancer_with_variants.csv',
         'region/3utr_with_variants.csv', 'region/5utr_with_variants.csv',
         'region/coding_sequences.csv', 'region/introns_with_variants.csv']
# qs = ['Is this an enhancer sequence?',
#       'Could this sequence function as an enhancer?',
#       'Is it possible that this DNA sequence acts as an enhancer?',
#       'Does this sequence serve the role of an enhancer in the genome?',
#       'Can we identify this sequence as an enhancer?',
#       'Is this sequence characterized as an enhancer?',
#       'Might this DNA sequence be classified as an enhancer?',
#       'Does this particular sequence have the properties of an enhancer?',
#       'Would this sequence be considered an enhancer?',
#       'Does this sequence exhibit the typical features of an enhancer?',
#       'Is there a chance this sequence is an enhancer?']
# qs = ['Is this a promoter sequence?',
#       'Could this sequence be acting as a promoter?',
#       'Does this sequence have the characteristics of a promoter?',
#       'Can this sequence be identified as a promoter in genetic transcription?',
#       'Might this DNA sequence function in the capacity of a promoter?',
#       'Is it possible for this DNA sequence to serve as a promoter?',
#       'Do the features of this sequence align with those of known promoters?',
#       'Is this sequence functioning in the role of a promoter?',
#       'Is this sequence responsible for starting the transcription process?',
#       'Is this sequence located upstream of a gene?',
#       'Is this sequence upstream or downstream of a particular gene?'
#       ]

# qs = ['Is this a DNA sequence belonging to open chromatin regions?',
#       'Does this sequence reside within regions of open chromatin?',
#       'Does this DNA fragment correspond to an open chromatin area?',
#       'Is this sequence a component of open chromatin regions?',
#       'Is it possible that this sequence is situated within an open chromatin region?',
#       'Does this sequence represent a section of the open chromatin region?',
#       'Does this particular sequence align with open chromatin regions?',
#       'Would this sequence be considered part of the open chromatin regions?',
#       'Does the location of this sequence correspond to open chromatin areas?',
#       'Is this an ocr sequence?',
#       ]

# qs = ['Is this a splice site sequence?',
#       'Does this sequence indicate a splice site?',
#       'Can this sequence be identified as a splice site?',
#       'Might this sequence represent a splice site?',
#       'Does this represent a sequence found at splice sites?',
#       'Is this sequence associated with splice sites?',
#       'Can this sequence be recognized as a splice site sequence?',
#       'Is this sequence reflective of a splice site?',
#       'Is there evidence to suggest this is a splice site sequence?',
#       'Is this sequence characteristic of a splice site?',
#       ]


# bs_df = pd.read_csv(bs_ds)
# bs_df = bs_df.iloc[:2000][['chr', 'start', 'end', 'seq']]
# bs_df['question'] = np.random.choice(qs, size=len(bs_df))
# bs_df['ans'] = 'yes'

# frames = [bs_df]
# for file in ot_ds:
#     df = pd.read_csv(file)
#     df = df.sample(n=250)[['chr', 'start', 'end', 'seq']]
#     df['question'] = np.random.choice(qs, size=len(df))
#     df['ans'] = 'no'
#     frames.append(df)

# result_df = pd.concat(frames)
# existing_df = pd.read_csv('dataset.csv')
# combined_df = pd.concat([existing_df, result_df])
# combined_df.to_csv('dataset.csv', index=False)

# ————————————————————————————————————————————————

def read_csv_partial(file, read_front_80=True):
    total_rows = sum(1 for _ in open(file)) - 1

    if read_front_80:
        nrows_to_read = int(total_rows * 0.8)
        df = pd.read_csv(file, nrows=nrows_to_read)
    else:
        rows_to_skip = int(total_rows * 0.8)
        df = pd.read_csv(file, skiprows=rows_to_skip + 1)

    return df


with open('index.txt', 'r') as file:
    data = file.read().split('*')

frames = []
stats_frames = []

for i, group in enumerate(data):
    questions = []
    answers = {}
    if 'q:' in group and 'a:' in group:
        q_part = group.split('q:')[1].split('a:')[0].strip()
        a_part = group.split('a:')[1].strip()

        questions = [q.strip() for q in q_part.split('\n') if q.strip()]

        for a in a_part.split('\n'):
            if ':' in a:
                key, value = a.split(':')
                answers[key.strip()] = value.strip()

        all_sampled_dfs = []

        for file, answer in answers.items():

            sampled_df = read_csv_partial(file, read_front_80=True)

            sampled_df['question'] = np.random.choice(
                questions, size=len(sampled_df))

            if re.match(r'\[.*\]', answer):
                column_name = re.findall(r'\[(.*)\]', answer)[0]

                if column_name in sampled_df.columns:
                    # sampled_df['ans'] = df[column_name].values[:len(
                    #     sampled_df)]

                    sampled_df['ans'] = sampled_df[column_name]
                else:
                    print(
                        f"Warning: Column '{column_name}' not found in file {file}. Skipping file.")

            else:
                sampled_df['ans'] = answer

            all_sampled_dfs.append(
                sampled_df[['chr', 'start', 'end', 'seq', 'question', 'ans']])

        all_sampled_df = pd.concat(all_sampled_dfs, ignore_index=True)

        value_counts = all_sampled_df['ans'].value_counts()

        if i == 24 or i == 31:
            target_values = ['N', 'C->T', 'G->A', 'A->G', 'T->C']

            min_count = value_counts[target_values].min()

            sampled_dfs = []
            for value in target_values:
                sampled_dfs.append(all_sampled_df[all_sampled_df['ans'] == value].sample(
                    n=min_count, random_state=42))

            others_df = all_sampled_df[~all_sampled_df['ans'].isin(
                target_values)]
            others_sampled = others_df.sample(n=min_count, random_state=42)
            others_sampled['ans'] = 'others'

            all_sampled_df = pd.concat(
                sampled_dfs + [others_sampled], ignore_index=True)

        elif i == 25:
            target_values = [0, 1, 2]

            min_count = value_counts[target_values].min()

            sampled_dfs = []
            for value in target_values:
                sampled_dfs.append(all_sampled_df[all_sampled_df['ans'] == value].sample(
                    n=min_count, random_state=42))

            all_sampled_df = pd.concat(sampled_dfs, ignore_index=True)

        elif i == 30:
            top_4_values = value_counts.nlargest(4).index.tolist()

            min_count = value_counts[top_4_values].min()

            sampled_dfs = []
            for value in top_4_values:
                sampled_dfs.append(all_sampled_df[all_sampled_df['ans'] == value].sample(
                    n=min_count, random_state=42))

            all_sampled_df = pd.concat(sampled_dfs, ignore_index=True)

        else:
            min_count = value_counts.min()

            sampled_dfs = [all_sampled_df[all_sampled_df['ans'] == value].sample(
                n=min_count, random_state=42) for value in value_counts.index]

            all_sampled_df = pd.concat(sampled_dfs, ignore_index=True)

        final_counts = all_sampled_df['ans'].value_counts()
        total = len(all_sampled_df)

        stats_df = pd.DataFrame({
            'group_index': i,
            'ans_category': final_counts.index,
            'count': final_counts.values,
            'proportion': final_counts.values / total
        })

        stats_frames.append(stats_df)

        frames.append(all_sampled_df)

final_df = pd.concat(frames)
final_df.to_csv('train_sampled.csv', index=False)

final_stats_df = pd.concat(stats_frames, ignore_index=True)
final_stats_df.to_csv('train_sampled_stats.csv', index=False)


# ————————————————————————————————————————————————

# # create labels
# # variant_yes_or_no, more_than_one_variant, variant, variant_number

# def process_csv(file_path):
#     df = pd.read_csv(file_path)
#     df['variant_yes_or_no'] = df['variants'].apply(
#         lambda x: 'yes' if x != 'N' else 'no')

#     df['more_than_one_variant'] = df['variants'].apply(
#         lambda x: 'yes' if x != 'N' and len(ast.literal_eval(x)) > 1 else 'no')

#     def select_random_variant(variants):
#         if variants == 'N':
#             return 'N'
#         else:
#             variants_list = ast.literal_eval(variants)
#             random_variant = random.choice(variants_list)
#             ref = random_variant[1]
#             alt = random_variant[2][0]
#             return f'{ref}->{alt}'

#     df['variant'] = df['variants'].apply(select_random_variant)

#     df['variant_number'] = df['variants'].apply(
#         lambda x: len(ast.literal_eval(x)) if x != 'N' else 0)

#     df.to_csv(file_path, index=False)
#     print(f"Processed file: {file_path}")


# csv_files = ['regulatory/ocr_with_variants.csv', 'regulatory/donor_with_variants.csv', 'regulatory/acceptor_with_variants.csv',
#              'regulatory/promoter_with_variants.csv', 'regulatory/enhancer_with_variants.csv',
#              'region/3utr_with_variants.csv', 'region/5utr_with_variants.csv',
#              'region/introns_with_variants.csv']

# for csv_file in csv_files:
#     process_csv(csv_file)

# ————————————————————————————————————————————————


# df = pd.read_csv('disease/cancer_with_variants.csv')

# df.rename(columns={'CHR_ID': 'chr'}, inplace=True)

# df['chr'] = 'chr' + df['chr'].astype(str)

# df.to_csv('disease/cancer_with_variants.csv', index=False)

# ————————————————————————————————————————————————

# df = pd.read_csv('disease/cancer_with_variants.csv')

# df['cancer_related'] = df['variants'].apply(
#     lambda x: 'yes' if x != 'N' else 'no')

# df['cancer_kind'] = df.apply(lambda row: 'N' if row['variants'] == 'N' else row['MAPPED_TRAIT'].split(
#     ',')[0] if pd.notnull(row['MAPPED_TRAIT']) else '', axis=1)


# def get_variant_for_chr_pos(variants_str, chr_pos):
#     try:
#         variants = ast.literal_eval(variants_str)
#         for variant in variants:
#             if variant[0] == chr_pos:
#                 return f"{variant[1]}->{variant[2][0]}"
#         return 'N'
#     except:
#         return 'N'


# df['cancer_related_variant'] = df.apply(
#     lambda row: get_variant_for_chr_pos(row['variants'], row['CHR_POS']), axis=1)

# df.to_csv('disease/cancer_with_variants.csv', index=False)

# ————————————————————————————————————————————————


# df = pd.read_csv('species/worm_seq.csv')

# df_filtered = df[df['label'] != 0]

# df_filtered['chr'] = '*'
# df_filtered['start'] = '*'
# df_filtered['end'] = '*'

# df_filtered.to_csv('species/worm_seq.csv', index=False)
