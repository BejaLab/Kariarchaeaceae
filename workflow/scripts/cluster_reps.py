from Bio import SeqIO
from pandas import read_csv

ref_files = snakemake.input['refs']
fasta_file = snakemake.input['fasta']
clstr_file = snakemake.input['clstr']

output_file = str(snakemake.output)

queries = {}
for ref_file in ref_files:
    for record in SeqIO.parse(ref_file, 'fasta'):
        queries[record.id] = True

selected = {}
with open(clstr_file) as file:
    best_match = None
    contains_queries = False
    for line in file:
        if line.startswith('>'):
            if not contains_queries:
                selected[best_match] = True
            best_match = None
            contains_queries = False
        else:
            num, length, seqid_str, *rest, identity = line.split()
            seqid = seqid_str[1:-3]
            if seqid in queries:
                contains_queries = True
            elif identity == '*':
                best_match = seqid
    if best_match is not None and not contains_queries:
        selected[best_match] = True

with open(output_file, 'w') as file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in selected or record.id in queries:
            SeqIO.write(record, file, 'fasta')
