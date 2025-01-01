from Bio import SeqIO
from pandas import read_csv

matches_file = snakemake.input['matches']
blast_file = snakemake.input['blast']
clstr_file = snakemake.input['clstr']

output_file = str(snakemake.output)

col_names = snakemake.params['cols']
max_ident = snakemake.params['max_ident']

blast = (
    read_csv(blast_file, sep = '\t', names = col_names)
        .sort_values(by = ['bitscore'], ascending = False)
        .query('pident <= @max_ident')
        .drop_duplicates(subset = 'sseqid')
)
queries = set(blast['qseqid'].to_list())
selected = { query: True for query in queries }
ranks = { sseqid: i for i, sseqid in enumerate(blast['sseqid'].to_list()) }

with open(clstr_file) as file:
    best_match = None
    contains_queries = False
    for line in file:
        if line.startswith('>'):
            if best_match is not None and not contains_queries:
                selected[best_match] = True
            best_rank = len(ranks)
            best_match = None
            contains_queries = False
        else:
            num, length, sseqid_str, *identity = line.split()
            sseqid = sseqid_str[1:-3]
            if sseqid in queries:
                contains_queries = True
            elif sseqid in ranks:
                if ranks[sseqid] < best_rank:
                    best_rank = ranks[sseqid]
                    best_match = sseqid
    if best_match is not None and not contains_queries:
        selected[best_match] = True

with open(output_file, 'w') as file:
    for record in SeqIO.parse(matches_file, 'fasta'):
        if record.id in selected:
            SeqIO.write(record, file, 'fasta')
