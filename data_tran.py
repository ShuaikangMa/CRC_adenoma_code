import os
import pandas as pd

def extract_genus_abundance(file_path):
    genus_abundance = {}
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            name = parts[5].strip()
            read_count = int(parts[1])
            if rank == 'G':
                genus_abundance[name] = read_count
            if rank == 'G1':
                 genus_abundance[name] = genus_abundance.get(name, 0) + read_count
                 genus_abundance['unclassified'] = genus_abundance.get('unclassified', 0) + read_count
                 if 'unclassified' in name:
                     genus_name = name.split('unclassified ')[-1].strip()
                     genus_abundance[genus_name] = genus_abundance.get(genus_name, 0) + read_count
                 else:
                     genus_abundance['unclassified'] = genus_abundance.get('unclassified', 0) + read_count
             if rank == 'G2':
                 genus_abundance[name] = genus_abundance.get(name, 0) + read_count

    return pd.Series(genus_abundance)

def merge_all_samples(folder):
  
    all_data = []
    sample_names = []

    for file in os.listdir(folder):
        if file.endswith('.bracken.kreport'):
            file_path = os.path.join(folder, file)
            sample_name = file.split('.')[0]
            sample_names.append(sample_name)
            genus_series = extract_genus_abundance(file_path)
            all_data.append(genus_series)

    df = pd.concat(all_data, axis=1)
    df.columns = sample_names
    df = df.fillna(0).astype(int)
    return df

folder_path = './down'
output_csv = './merged_g_abundance.csv'

merged_df = merge_all_samples(folder_path)
merged_df.to_csv(output_csv)

print("end:", output_csv)
