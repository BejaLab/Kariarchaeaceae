from Bio import SeqIO

def split_fasta(input_file, output_file, gap_size):
    gap = 'N' * gap_size
    with open(output_file, 'w') as f:
        for record in SeqIO.parse(input_file, 'fasta'):
            seq = str(record.seq).upper()
            splits = seq.split(gap)  # split by stretches of 'N'
            if len(splits) > 1:
                for i, split_seq in enumerate(splits, start = 1):
                    if split_seq: # Skip empty sequences
                        new_record = record[:]
                        new_record.id = f"{record.id}.{i}"
                        new_record.seq = split_seq
                        SeqIO.write(new_record, f, 'fasta')
            else:
                SeqIO.write(record, f, 'fasta')

input_file = str(snakemake.input)   # input FASTA file
output_file = str(snakemake.output) # output file name
gap_size = snakemake.params['gap']

split_fasta(input_file, output_file, gap_size)
