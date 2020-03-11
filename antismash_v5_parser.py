#!/usr/bin/env python3
import argparse
import json
import os
import re
import pandas as pd
from glob import glob
parser = argparse.ArgumentParser(description="""Extracts regions from antismash version 5 output and saves to a csv file, 
and to 4 column bed file. Csv output will be similar to this:

+----------+-------------+------------------+--------+--------+-------------+----------------------------+------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  contig  | contig_edge |     product      | start  |  end   | region_name | Most similar known cluster | Similarity |     known_cluster_type      |                                                                      genes_in_cluster                                                                       |
+----------+-------------+------------------+--------+--------+-------------+----------------------------+------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| contig_3 | FALSE       | NRPS-like        | 713302 | 756656 |         3.1 |                            |            |                             | kd17_993_g;kd17_994_g;kd17_995_g;kd17_996_g;kd17_997_g;kd17_998_g;kd17_999_g;kd17_1000_g;kd17_1001_g                                                        |
| contig_4 | FALSE       | T1PKS            |  72613 | 112751 |         4.1 |                            |            |                             | kd17_1075_g;kd17_1076_g;kd17_1077_g;kd17_1078_g;kd17_1079_g;kd17_1080_g;kd17_1081_g;kd17_1082_g;kd17_1083_g;kd17_1084_g;kd17_1085_g;kd17_1086_g;kd17_1087_g |
| contig_4 | FALSE       | NRPS-like, T1PKS | 467152 | 545384 |         4.2 | asparasone A               | 75%        | Polyketide:Iterative type I | kd17_1222_g;kd17_1223_g;kd17_1224_g;kd17_1225_g;kd17_1226_g;kd17_1227_g;kd17_1228_g;kd17_1229_g;kd17_1230_g                                                 |
| contig_6 | FALSE       | NRPS             | 109529 | 153546 |         6.1 | astellolide A              | 100%       | Terpene                     | kd17_1610_g;kd17_1611_g;kd17_1612_g;kd17_1613_g;kd17_1614_g;kd17_1615_g;kd17_1616_g;kd17_1617_g                                                             |
| contig_7 | FALSE       | fungal-RiPP      | 290143 | 331036 |         7.1 | ustiloxin B                | 78%        | RiPP:Thiopeptide            | kd17_1896_g;kd17_1897_g;kd17_1898_g;kd17_1899_g;kd17_1900_g;kd17_1901_g                                                                                     |
+----------+-------------+------------------+--------+--------+-------------+----------------------------+------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
""", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('antismash_v5_directory')
args = parser.parse_args()
json_file = glob(args.antismash_v5_directory + '/*.json')[0]
html_index = os.path.join(args.antismash_v5_directory, 'index.html')
knownclusterblasts = glob(os.path.join(args.antismash_v5_directory, 'knownclusterblast') + '/*.txt')
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
def get_genes_in_cluster(clusterblast):
    region_name = re.search('(\d+?)_c(\d+?).txt', os.path.basename(clusterblast))
    region_name = region_name.group(1) + '.' + region_name.group(2)
    genes_in_cluster = []
    for idx, line in enumerate(open(clusterblast)):
        if idx > 2:
            if line.startswith('\n'): 
                break
            genes_in_cluster.append(line.split('\t')[0])
    return (region_name, ';'.join(genes_in_cluster))
genes_in_cluster = pd.DataFrame.from_records([get_genes_in_cluster(clusterblast) for clusterblast in knownclusterblasts], columns=['region_name', 'genes_in_cluster'])
regions_df = regions_df.merge(known_clusters, on='region_name', how='left')
regions_df = regions_df.merge(genes_in_cluster, on='region_name', how='left')
regions_df.to_csv(os.path.join(args.antismash_v5_directory,base+'.csv'), index=False)
regions_df[['contig','start', 'end', 'region_name']].to_csv(os.path.join(args.antismash_v5_directory,base+'.bed'), index=False, header=False, sep='\t')
