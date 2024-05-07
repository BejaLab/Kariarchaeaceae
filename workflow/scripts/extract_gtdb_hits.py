import csv
from sys import maxsize
import pandas as pd

input_file = str(snakemake.input)
output_file = str(snakemake.output)

params = dict(snakemake.params)

cols = params['cols']
min_id = params['pident']
max_evalue = params['evalue']
top_hits = params['top_hits']
ingroup_ratio = params['ingroup_ratio']
ingroup = params['ingroup']

def read_tsv_line(file, cols):
    for line in csv.DictReader(file, delimiter = '\t', fieldnames = cols):
        line['pident'] = float(line['pident'])
        line['bitscore'] = float(line['bitscore'])
        line['length'] = int(line['length'])
        line['evalue'] = float(line['evalue'])
        yield line

counts = {}
hits = {}
best_score = {}
csv.field_size_limit(maxsize)

with open(input_file) as file:
    for line in read_tsv_line(file, cols):
        id = line['qseqid']
        if line['pident'] >= min_id and line['evalue'] <= max_evalue:
            if id not in counts:
                counts[id] = { True: 0, False: 0 }
            if id not in best_score:
                best_score[id] = line['bitscore']
            if counts[id][True] + counts[id][False] < top_hits:
                is_ingroup = ingroup in line['stitle']
                counts[id][is_ingroup] += 1
                if is_ingroup and id not in hits:
                    hits[id] = line

hits = { id: hits[id] for id, n in counts.items() if n[True] / (n[True] + n[False]) >= ingroup_ratio }
df = pd.DataFrame.from_dict(hits, orient = "index")
df.to_csv(output_file, index = False)
