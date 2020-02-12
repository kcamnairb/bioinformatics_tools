#!/usr/bin/env python3
import argparse
import json
import os
import pandas as pd
from glob import glob
parser = argparse.ArgumentParser(description="""Extracts regions from antismash version 5 json output and saves to a csv file, 
and to 4 column bed file. Csv output will be similar to this:
+----------+-------------+------------------+--------+--------+-------------+----------------------------+------------+-----------------------------+
|  contig  | contig_edge |     product      | start  |  end   | region_name | Most similar known cluster | Similarity |     known_cluster_type      |
+----------+-------------+------------------+--------+--------+-------------+----------------------------+------------+-----------------------------+
| contig_3 | FALSE       | NRPS-like        | 713302 | 756656 |         3.1 |                            |            |                             |
| contig_4 | FALSE       | T1PKS            |  72613 | 112751 |         4.1 |                            |            |                             |
| contig_4 | FALSE       | NRPS-like, T1PKS | 467152 | 545384 |         4.2 | asparasone A               | 75%        | Polyketide:Iterative type I |
| contig_4 | FALSE       | NRPS             | 580949 | 628593 |         4.3 |                            |            |                             |
| contig_6 | FALSE       | NRPS             | 109529 | 153546 |         6.1 | astellolide A              | 100%       | Terpene                     |
+----------+-------------+------------------+--------+--------+-------------+----------------------------+------------+-----------------------------+
""", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('antismash_v5_directory')
args = parser.parse_args()
json_file = glob(args.antismash_v5_directory + '/*.json')[0]
html_index = os.path.join(args.antismash_v5_directory, 'index.html')
as5 = json.load(open(json_file))
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
regions_df['start'] = regions_df['start'].str.replace('<|>','')
regions_df['end'] = regions_df['end'].str.replace('<|>','')
def contig_and_region_num_to_name(row):
    contig_num = contig_to_num[row['contig']]
    return str(contig_num) + '.' + row['region_number']
regions_df['region_name'] = regions_df.apply(contig_and_region_num_to_name, axis=1)
regions_df = regions_df.drop(['candidate_cluster_numbers', 'probabilities', 'subregion_numbers', 'tool', 'rules', 'region_number', 'location'], axis=1)
base = os.path.basename(args.antismash_v5_directory).replace('json', '')
def get_known_clusterblast_results(html_index):
    dfs = []
    for df in pd.read_html(html_index)[:-1]:
        if df.shape[1] == 7:
            df.columns = ['region_name','Type','From','To','Most similar known cluster','known_cluster_type','Similarity']
            dfs.append(df)
        else:
            dfs.append(df)
    df = pd.concat(dfs)
    df.region_name = df.region_name.str.replace('Region&nbsp','')
    df = df[['region_name','Most similar known cluster','Similarity','known_cluster_type']]         
    return df
known_clusters = get_known_clusterblast_results(html_index)
regions_df = regions_df.merge(known_clusters, on='region_name', how='left')
regions_df.to_csv(os.path.join(args.antismash_v5_directory,base+'.csv'), index=False)
regions_df[['contig','start', 'end', 'region_name']].to_csv(os.path.join(args.antismash_v5_directory,base+'.bed'), index=False, header=False, sep='\t')
