from bioformats import read_hmmsearch
import pandas as pd

input_file = str(snakemake.input)
output_file = str(snakemake.output)

rows = []
with open(input_file) as file:
    for record in read_hmmsearch(file):
        row = {}
        row['stitle'] = record['description']
        row['sseqid'] = record['seq_name']
        row['clade'] = row['subclade'] = record['profile']['name']
        row['evalue'] = record['full']['evalue']
        rows.append(row)

pd.DataFrame(rows).to_csv(output_file)
