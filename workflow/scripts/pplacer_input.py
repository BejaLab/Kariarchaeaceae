import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import tempfile
from io import StringIO
from collections import defaultdict

fasta_files  = snakemake.input['fasta']
blast_files  = snakemake.input['blast']
query_file   = snakemake.input['query']
msa_dir      = snakemake.input['msas']
exclude_file = snakemake.input['exclude']
output_fasta_file = snakemake.output['fasta']
genes_file = snakemake.output['genes']

genomes = snakemake.params['genomes']
evalue_threshold = snakemake.params['evalue']
min_genes = snakemake.params['min_genes']
query_scaffold = snakemake.params['scaffold']

genome_set = set(genomes + [ query_scaffold ])

# Exclude ORF(s) matching rhodopsins
exclude = []
with open(exclude_file) as file:
    for line in file:
        query, match, *rest = line.split()
        exclude.append(match)

alignment = { genome: '' for genome in genomes }
alignment[query_scaffold] = ''

coordinates = {}

data = {}
for record in SeqIO.parse(query_file, 'fasta'):
    data[record.id] = { query_scaffold: record }

with open(genes_file, 'w') as genes_file:
    for msa_file in os.listdir(msa_dir):
        if msa_file.endswith('.aln'):
            record = None
            msa_path = os.path.join(msa_dir, msa_file)
            genomes_with_gene = set()
            for record in SeqIO.parse(msa_path, 'fasta'):
                alignment[record.id] += record.seq
                genes_file.write(f"{record.id}\t{msa_file}\treference\n")
                genomes_with_gene.add(record.id)
            gaps = '-' * len(record.seq)
            for genome in genome_set - genomes_with_gene:
                alignment[genome] += gaps

    coordinates = {}
    for i, genome in enumerate(genomes):
        genes = {}
        with open(blast_files[i]) as file:
            for line in file:
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
                if qseqid not in genes and float(evalue) <= evalue_threshold and qseqid not in exclude and sseqid not in exclude:
                    genes[qseqid] = sseqid
                    start, stop = int(qstart) - 1, int(qend)
                    if qseqid in coordinates:
                        coordinates[qseqid] = min(start, coordinates[qseqid][0]), max(stop, coordinates[qseqid][1])
                    else:
                        coordinates[qseqid] = int(qstart), int(qend)
        if len(genes) > 0:
            markers = dict([ reversed(g) for g in genes.items() ])
            for record in SeqIO.parse(fasta_files[i], 'fasta'):
                if record.id in markers:
                    marker = markers[record.id]
                    data[marker][genome] = record

    for marker, marker_data in data.items():
        if len(marker_data) > min_genes:
            with tempfile.NamedTemporaryFile(mode = "w") as fasta_file:
                genomes_with_gene = set()
                for genome, record in marker_data.items():
                    if record.id in coordinates:
                        start, stop = coordinates[record.id]
                        record.seq = record.seq[start:stop].strip('X')
                    record.id = genome
                    SeqIO.write(record, fasta_file, 'fasta')
                    genes_file.write(f"{genome}\t{marker}\tquery\n")
                    genomes_with_gene.add(genome)
                fasta_file.flush()
                result = subprocess.run([ 'mafft', '--auto', fasta_file.name ], capture_output = True, text = True, check = True)
                record = None
                for record in SeqIO.parse(StringIO(result.stdout), "fasta"):
                    alignment[record.id] += record.seq
                gap = '-' * len(record.seq)
                for genome in genome_set - genomes_with_gene:
                    alignment[genome] += gap

with open(output_fasta_file, 'w') as file:
    for genome, seq in alignment.items():
        record = SeqRecord(Seq(seq), id = genome, description = '')
        SeqIO.write(record, file, 'fasta')
