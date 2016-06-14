import argparse
parser = argparse.ArgumentParser(description='''Takes SRA numbers (beggining with "SRR") and ouputs the url to download them from ENA.
Conversion done according to http://www.ebi.ac.uk/ena/browse/read-download''')
parser.add_argument('sras', nargs='+')
args = parser.parse_args()
sra = args.sras
def sra2ena_url(sra):   
   dir1 = sra[:6]
   base_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/' + dir1 + '/'
   dir2 = ''
   digits = [x for x in sra if x.isdigit()]
   if len(digits) == 6:
      url = base_url + sra
   elif len(digits) == 7:
      dir2 = '00' + digits[-1]
      url = base_url + dir2 + '/' + sra
   elif len(digits) == 8:
      dir2 = '0' + digits[-2:]
      url = base_url + dir2 + '/' + sra
   elif len(digits) == 9:
      dir2 = digits[-3:]
      url = base_url + dir2 + '/' + sra
   return url
for id in sra:
    print(sra2ena_url(id))
