from Bio import SeqIO

fasta_file = snakemake.input['fasta']
trimal_file = snakemake.input['trimal']
output_file = str(snakemake.output)

with open(trimal_file) as file:
    header, *pos = file.read().split()

start = int(pos[1].strip(',')) - 1
stop = int(pos[-1].strip(',')) + 1

with open(output_file, 'w') as out:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        record.seq = record.seq[start:stop]
        SeqIO.write(record, out, 'fasta')
