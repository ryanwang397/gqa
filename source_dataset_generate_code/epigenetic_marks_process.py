# import pandas as pd

# def process_fna_to_csv(input_file):
#     # Dictionary to hold sequences categorized by their identifier (0, 1, 2, etc.)
#     file_data = {"0": [], "1": [], "2": []}

#     # Read the entire content of the fna file
#     with open(input_file, 'r') as f:
#         content = f.read()

#     # Split the content into blocks based on '>'
#     blocks = content.strip().split('>')

#     for block in blocks:
#         if block.strip():  # Skip any empty blocks
#             lines = block.strip().splitlines()

#             # First line is the header like ***|0 or chr10:70478753-70479053|1
#             header = lines[0].strip()
#             # Combine remaining lines as sequence
#             seq = ''.join(line.strip() for line in lines[1:])

#             # Check the header and extract key (0, 1, 2)
#             if not header.startswith("***"):
#                 # For cases like ***|0 or ***|2
#                 parts = header.split("|")
#                 current_chr, key = "*", parts[-1]
#                 current_start, current_end = "*", "*"
#             else:
#                 # For cases like chr10:70478753-70479053|1
#                 parts = header.split(":")
#                 current_chr = parts[0]
#                 start_end = parts[1].split("-")
#                 current_start = start_end[0]
#                 current_end = start_end[1].split(
#                     "|")[0]  # In case of chr with |key
#                 key = header.split("|")[-1]  # Get the last part as key

#             # Append the processed data to the corresponding list in file_data
#             if key not in file_data:
#                 file_data[key] = []
#             file_data[key].append(
#                 [current_chr, current_start, current_end, seq])

#     # Convert file data to DataFrames and save to CSV files
#     for key, data in file_data.items():
#         if data:
#             df = pd.DataFrame(data, columns=["chr", "start", "end", "seq"])
#             df.to_csv(f"output_{key}.csv", index=False)

#     print(
#         f"CSV files have been successfully generated for keys: {', '.join(file_data.keys())}")

# process_fna_to_csv('more_tasks/epigenetic_marks/yeast/H4ac.fna')

# #-----------------------------------------------------------------------------

# import os
# import pandas as pd
# import random

# input_folder = 'epigenetic_marks/yeast_ini'
# output_folder = 'epigenetic_marks/yeast'

# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)

# for filename in os.listdir(input_folder):
#     if filename.endswith(".csv"):
#         file_path = os.path.join(input_folder, filename)
#         df = pd.read_csv(file_path)

#         new_data = []
#         for index, row in df.iterrows():
#             start = row['start']
#             end = row['end']
#             seq = row['seq']
#             chr = row['chr']

#             seq_len = len(seq)
#             mid = seq_len // 2

#             sub_seq_len = random.randint(200, 500)

#             half_len = sub_seq_len // 2
#             new_start_idx = max(0, mid - half_len)
#             new_end_idx = min(seq_len, mid + half_len)

#             new_seq = seq[new_start_idx:new_end_idx]

#             # new_start = start + new_start_idx
#             # new_end = start + new_end_idx - 1

#             new_data.append({
#                 'chr': chr,
#                 'start': start,
#                 'end': end,
#                 'seq': new_seq
#             })

#         new_df = pd.DataFrame(new_data)

#         output_file_path = os.path.join(output_folder, filename)
#         new_df.to_csv(output_file_path, index=False)

# print("Processing complete. Modified files are saved in:", output_folder)
