from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pandas import read_csv

tblastn_file = snakemake.input['tblastn']
fasta_file = snakemake.input['fasta']
output_file = str(snakemake.output)

params = dict(snakemake.params)

cols = params['cols']
pident_max = params['pident_max'] if 'pident_max' in params else 100
split_seqid = params['split_seqid'] if 'split_seqid' in params else False
evalue_max = params['evalue_max'] if 'evalue_max' in params else 1

tblastn_results = read_csv(tblastn_file, names = cols, sep = '\t').sort_values([ 'bitscore' ], ascending = [ False ])

pad = 999

def find_orf_plus(seq, start, end):
    frame = start % 3
    orf_start = start - pad
    if orf_start < 0:
        orf_start = frame
    i_from = orf_start
    stops = 0
    for i in range(i_from, len(seq), 3):
        codon = seq[i : i + 3]
        orf_end = i + 2
        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
            if i <= start:
                orf_start = i + 3
                stops = 0
            elif i >= end:
                orf_end = i - 1
                break
            else:
                stops += 1
    return orf_start, orf_end, stops

def find_orf_minus(seq, start, end):
    last_i = len(seq) - 1
    frame = (last_i - start) % 3
    orf_start = start + pad
    if orf_start > last_i:
        orf_start = last_i - frame
    i_from = orf_start - 2
    stops = 0
    for i in range(i_from, -1, -3):
        codon = seq[i : i + 3]
        orf_end = i
        if codon == 'TTA' or codon == 'CTA' or codon == 'TCA':
            if i >= start:
                orf_start = i - 1
                stops = 0
            elif i <= end:
                orf_end = i + 3
                break
            else:
                stops += 1
    return orf_start, orf_end, stops

def extract_sequence(full_record, start, end):
    seq = full_record.seq
    if start < end:
        orf_start, orf_end, stops = find_orf_plus(seq, start, end)
        orf_seq = seq[orf_start:orf_end + 1]
        strand = '+'
    else:
        orf_start, orf_end, stops = find_orf_minus(seq, start, end)
        orf_seq = seq[orf_end:orf_start + 1].reverse_complement()
        strand = '-'

    record = SeqRecord(orf_seq, id = f"{full_record.id}_{orf_start}_{orf_end}", description = '')
    return record, stops

fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))  # Use your fasta file name

ids = {}
with open(output_file, "w") as output_handle:
    for index, row in tblastn_results.iterrows():
        sseqid, sstart, send = row['sseqid'], row['sstart'], row['send']
        subject_id = sseqid.split('|')[1] if split_seqid else sseqid
        record, stops = extract_sequence(fasta_sequences[subject_id], row['sstart'] - 1, row['send'] - 1)
        if record.id not in ids and not stops and row['pident'] <= pident_max and row['evalue'] <= evalue_max:
            SeqIO.write(record, output_handle, "fasta")
        ids[record.id] = True
