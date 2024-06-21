import re

input_file = str(snakemake.input)

nodes_file = snakemake.output['nodes']
names_file = snakemake.output['names']
mapping_file = snakemake.output['mapping']
merged_file = snakemake.output['merged']
delnodes_file = snakemake.output['delnodes']

# Initialize dictionaries and counters
ids = {"root": 1}
rank = {"d": "superkingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
taxCnt = 1

# Process each line in the input file
with open(input_file) as file:
    with open(nodes_file, 'w') as nodes_dmp, open(names_file, 'w') as names_dmp, open(mapping_file, 'w') as mapping:
        nodes_dmp.write("1\t|\t1\t|\tno rank\t|\t-\t|\n")
        names_dmp.write("1\t|\troot\t|\t-\t|\tscientific name\t|\n")
        for line in file:
            if line.startswith(">"):
                parts = line.split()
                str_val = " ".join(parts[4:])
                a = str_val.split(";")
                prevTaxon = 1
                for item in a:
                    if item in ids:
                        prevTaxon = ids[item]
                    else:
                        taxCnt += 1
                        b = item.split("_")
                        nodes_dmp.write(f"{taxCnt}\t|\t{prevTaxon}\t|\t{rank[b[0]]}\t|\t-\t|\n")
                        names_dmp.write(f"{taxCnt}\t|\t{b[2]}\t|\t-\t|\tscientific name\t|\n")
                        ids[item] = taxCnt
                        prevTaxon = ids[item]

                seq_id = parts[0].replace(">", "")
                mapping.write(f"{seq_id}\t{ids[a[-1]]}\n")

open(merged_file, 'w').close()
open(delnodes_file, 'w').close()
