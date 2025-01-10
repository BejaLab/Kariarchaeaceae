from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

maf_file = str(snakemake.input)
fasta_file = snakemake.output['fasta']
bed_file = snakemake.output['bed']

genomes = snakemake.params['genomes']

def get_genome(record_id, genomes):
    found = []
    for genome in genomes:
        if record_id.startswith(genome + '.'):
            found.append(genome)
    assert found, f"No genome found for record {record_id}"
    assert len(found) == 1, f"More than one genome found for record {record_id}"
    return found[0]

alignment = defaultdict(str)
lens = defaultdict(int)
for records in AlignIO.parse(maf_file, "maf"):
    ref = records[0]
    ref_len = len(ref.seq)
    blank = Seq('-' * ref_len)
    aln = {}
    lens[ref.id] += ref_len
    for record in records:
        genome = get_genome(record.id, genomes)
        aln[genome] = record.seq
    for genome in genomes:
        present = genome in aln
        alignment[genome] += aln[genome] if present else blank

with open(bed_file, 'w') as bed:
    start = 0
    for ref_id, ref_len in lens.items():
        bed.write(f"{ref_id}\t{start}\t{start + ref_len}\n")
        start += ref_len

with open(fasta_file, 'w') as fasta:
    for genome, aln in alignment.items():
        record = SeqRecord(aln, id = genome, description = '')
        SeqIO.write(record, fasta, 'fasta')

