from pandas import read_csv
from io import StringIO
import pysam
from tempfile import NamedTemporaryFile
import os

min_cov = snakemake.params['min_cov']

bam_files = []
for bam_file in snakemake.input:
    data_str = pysam.coverage(bam_file)
    data = read_csv(StringIO(data_str), sep = "\t")
    covbases = sum(data["covbases"])
    if covbases >= min_cov:
        bam_files.append(bam_file)

assert bam_files, "No bam files passed the coverage threshold"

chunk_size = 1000

try:
    chunks = []
    for i in range(0, len(bam_files), chunk_size):
        chunk = bam_files[i:i + chunk_size]
        f = NamedTemporaryFile(delete = False, suffix = '.bam')
        pysam.merge("-f", "-o", f.name, *chunk)
        chunks.append(f.name)
    pysam.merge("-o", str(snakemake.output), *chunks)
except Exception as e:
    print("An error occurred: ", e)
finally:
    for chunk_file in chunks:
        try:
            os.unlink(chunk_file)
        except Exception as e:
            print(f"Unable to delete temporary file {chunk_file}: {e}")
