from Bio import SeqIO
from collections import defaultdict

input_file = str(snakemake.input)
output_file = str(snakemake.output)
residues = snakemake.params['residues']

ref_name = residues.pop("ref")

records = SeqIO.to_dict(SeqIO.parse(input_file, 'fasta'))

ref = records[ref_name]

aln_residues = defaultdict(list)

ref_pos = 0
for aln_pos, res in enumerate(ref.seq):
    ref_pos += res != '-'
    for feature, feature_pos in residues.items():
        if ref_pos in feature_pos:
            aln_residues[feature].append(aln_pos)

with open(output_file, 'w') as file:
    header = [ "label" ]
    header += aln_residues.keys()
    file.write("\t".join(header) + "\n")
    for record in records.values():
        file.write(record.id)
        for feature, aln_positions in aln_residues.items():
            seq = ''.join([ record.seq[aln_pos] for aln_pos in aln_positions ])
            file.write("\t" + seq)
        file.write("\n")
