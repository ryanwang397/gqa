import os
import csv

# 输入文件夹路径
input_folder = "RefSeq/plasmid/noncoding"
output_csv = "RefSeq_plasmid_noncoding.csv"

# 初始化存储结果的列表
output_rows = []

# 遍历文件夹中的每个文件
for file_name in os.listdir(input_folder):
    # 只处理 .fasta 文件
    if file_name.endswith(".fasta"):
        file_path = os.path.join(input_folder, file_name)
        with open(file_path, "r") as file:
            # 逐行读取文件内容
            for line in file:
                # 跳过注释行或标题行
                if line.startswith(">"):
                    continue
                # 去除行尾换行符并保存
                sequence = line.strip()
                output_rows.append(["*", "*", "*", sequence.upper()])

# 将结果写入CSV文件
with open(output_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    # 写入表头
    csvwriter.writerow(["chr", "start", "end", "seq"])
    # 写入数据行
    csvwriter.writerows(output_rows)

print(f"Data successfully written to {output_csv}")
