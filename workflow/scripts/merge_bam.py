import pysam, re
from pysam import AlignmentFile

output_file = str(snakemake.output)

trim = snakemake.params['trim']
quals = snakemake.params['quals']
do_sort = snakemake.params['sort']
platform = snakemake.params['platform']
sample = snakemake.params['sample']

header = None
read_groups = []
all_reads = []
for bam_file in snakemake.input:
    with AlignmentFile(bam_file, "rb") as input:
        if reads := list(input):
            all_reads += reads
            read_groups += input.header['RG']
            if not header:
                header = input.header.to_dict()
assert header, "No mapped reads"

if trim:
    for reference in header["SQ"]:
        reference["SN"] = re.sub(trim, '', reference["SN"])
if platform or sample:
    for rg in read_groups:
        if platform:
            rg['PL'] = platform
        if sample:
            rg['SM'] = sample
header['RG'] = read_groups

with AlignmentFile(output_file, "wb", header = header) as output:
    for read in all_reads:
        if quals and read.qual is None and read.query_sequence is not None:
            read.qual = quals * len(read.query_sequence)
        output.write(read)

if do_sort:
    pysam.sort("-o", output_file, output_file)
