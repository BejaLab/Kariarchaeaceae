import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import tempfile
from io import StringIO

fasta_files = snakemake.input['fasta']
blast_files = snakemake.input['blast']
query_file  = snakemake.input['query']
msa_file    = snakemake.input['msa']
output_fasta_file = snakemake.output['fasta']
output_part_file = snakemake.output['part']

genomes = snakemake.params['genomes']
evalue_threshold = snakemake.params['evalue']
min_genes = snakemake.params['min_genes']
query_genome = snakemake.params['query_genome']

aligned = {}
record = None
for record in SeqIO.parse(msa_file, 'fasta'):
    aligned[record.id] = record.seq
msa_end = len(record.seq)
aligned[query_genome] = '-' * msa_end
marker_coords = {}

data = {}
for record in SeqIO.parse(query_file, 'fasta'):
    data[record.id] = { query_genome: record }

for i, genome in enumerate(genomes):
    genes = {}
    with open(blast_files[i]) as file:
        for line in file:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
            if qseqid not in genes and float(evalue) <= evalue_threshold:
                genes[qseqid] = sseqid
    if len(genes) > 0:
        markers = dict([ reversed(g) for g in genes.items() ])
        for record in SeqIO.parse(fasta_files[i], 'fasta'):
            if record.id in markers:
                marker = markers[record.id]
                data[marker][genome] = record

for marker, marker_data in data.items():
    if len(marker_data) > min_genes:
        with tempfile.NamedTemporaryFile(mode = "w") as fasta_file:
            missing = list(aligned.keys())
            for genome, record in marker_data.items():
                record.id = genome
                SeqIO.write(record, fasta_file, 'fasta')
                missing.remove(genome)
            fasta_file.flush()
            result = subprocess.run([ 'mafft', '--auto', fasta_file.name ], capture_output = True, text = True, check = True)
            record = None
            for record in SeqIO.parse(StringIO(result.stdout), "fasta"):
                aligned[record.id] += record.seq
            aln_len = len(record.seq)
            for genome in missing:
                aligned[genome] += '-' * aln_len
            marker_coords[marker] = msa_end + 1, msa_end + aln_len
            msa_end += aln_len

with open(output_fasta_file, 'w') as file:
    for genome, seq in aligned.items():
        record = SeqRecord(Seq(seq), id = genome, description = '')
        SeqIO.write(record, file, 'fasta')
with open(output_part_file, 'w') as file:
    for marker, (start, end) in marker_coords.items():
        file.write(f"AA, {marker} = {start}-{end}\n")
