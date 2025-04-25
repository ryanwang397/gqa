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
