import pandas as pd
import csv
from sys import maxsize
import warnings

outfmt6_files = snakemake.input['outfmt6']
metadata_file = snakemake.input['metadata']
output_file = str(snakemake.output)

params = dict(snakemake.params)
cols = params['cols']
min_id = params['id'] if 'id' in params else 0
min_score = params['score'] if 'score' in params else 0
min_length = params['length'] if 'length' in params else 0
max_evalue = params['evalue'] if 'evalue' in params else 1
clades = params['clades'] if 'clades' in params else []

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    df = pd.read_excel(metadata_file)
metadata = df[[ 'label', 'clade', 'subclade', 'alias' ]].set_index('label').T.to_dict()

def read_tsv_line(file, cols):
    for line in csv.DictReader(file, delimiter = '\t', fieldnames = cols):
        line['pident'] = float(line['pident'])
        line['bitscore'] = float(line['bitscore'])
        line['length'] = int(line['length'])
        line['evalue'] = float(line['evalue'])
        yield line

hits = {}
csv.field_size_limit(maxsize)
for i, outfmt6_file in enumerate(outfmt6_files):
    with open(outfmt6_file) as file:
        for line in read_tsv_line(file, cols):
            sseqid = line['sseqid']
            qseqid = line['qseqid']
            if line['bitscore'] > min_score and line['pident'] > min_id and line['length'] > min_length and line['evalue'] < max_evalue:
                if sseqid not in hits or hits[sseqid]['bitscore'] < line['bitscore']:
                    hits[sseqid] = line
                    for key, value in metadata[qseqid].items():
                        hits[sseqid][key] = value

df = pd.DataFrame.from_dict(hits, orient = "index")
if clades:
    df = df.query('clade in @clades or subclade in @clades')
df.to_csv(output_file, index = False)
