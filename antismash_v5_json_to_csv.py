#!/usr/bin/env python3
import argparse
import json
import os
import pandas as pd
parser = argparse.ArgumentParser(description="""Extracts regions from antismash version 5 json output and saves to a csv file.
Output will be similar to this:
+----------+-------------+-----------+--------+--------+-------------+
| contig   | contig_edge | product   | start  | end    | region_name |
+----------+-------------+-----------+--------+--------+-------------+
| contig_3 | False       | NRPS-like | 713302 | 756656 | 3.1         |
+----------+-------------+-----------+--------+--------+-------------+
| contig_3 | False       | terpene   | 857436 | 878671 | 3.2         |
+----------+-------------+-----------+--------+--------+-------------+
| contig_4 | False       | T1PKS     | 72613  | 112751 | 4.1         |
+----------+-------------+-----------+--------+--------+-------------+
""", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('antismash_v5_json_file')
args = parser.parse_args()
as5 = json.load(open(args.antismash_v5_json_file))
regions = []
for contig in as5['records']:
    for f in contig['features']:
        if 'region_number' in f['qualifiers']:
            regions.append(dict(contig=contig['name'], location=f['location'], **f['qualifiers']))
regions_df = pd.DataFrame(regions)
contig_to_num = {record['name']:num for num, record in enumerate(as5['records'], 1)}
regions_df['product'] = regions_df['product'].apply(', '.join)
regions_df['contig_edge'] = regions_df['contig_edge'].apply(', '.join)
regions_df['region_number'] = regions_df['region_number'].apply(', '.join)
regions_df[['start','end']] = regions_df['location'].str.replace('[\[\]]', '').str.split(':', expand=True)
def contig_and_region_num_to_name(row):
    contig_num = contig_to_num[row['contig']]
    return str(contig_num) + '.' + row['region_number']
regions_df['region_name'] = regions_df.apply(contig_and_region_num_to_name, axis=1)
regions_df = regions_df.drop(['candidate_cluster_numbers', 'probabilities', 'subregion_numbers', 'tool', 'rules', 'region_number', 'location'], axis=1)
base = os.path.basename(args.antismash_v5_json_file).replace('json', '')
regions_df.to_csv(base+'csv', index=False)