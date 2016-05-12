import argparse
parser = argparse.ArgumentParser(description="Takes proteinortho output (or similarly structured list of orthologs) and creates input for making venn diagram from sites like http://bioinformatics.psb.ugent.be/webtools/Venn/")
parser.add_argument("proteinortho")
args = parser.parse_args()
orthologs_list = [line.strip().split('\t')[3:] for line in open(args.proteinortho,'r') if line.count('*') < len(line.strip().split('\t')[3:])]
header = orthologs_list.pop(0)
print(header)

venn_input = []
for idx, line in enumerate(orthologs_list):
    temp = []
    for gene in line:
        if gene == '*':
            temp.append('')
        else:
            temp.append(idx)
    venn_input.append(temp)
with open('venn_input.txt','w') as outfile:
    for line in venn_input:
        outfile.write(str(line)+'\n')

for idx, species in enumerate(header):
    with open(species,'w') as outfile:
        for line in venn_input:
            outfile.write(str(line[idx])+'\n')
sets = [set(genes) for genes in zip(*venn_input)]
print('Number of unique genes:')
for idx, species in enumerate(header):
    unique_set = sets[idx]
    for idx2, species2 in enumerate(header):
        if idx != idx2:
            unique_set = unique_set - sets[idx2]
    print(': '.join([species, str(len(unique_set))]))
