from sys import argv
import os
import os.path
import gzip

taxonomy_file = snakemake.input['taxonomy']
assembly_dir = snakemake.input['assemblies']
output_file = str(snakemake.output)
taxa = snakemake.params['taxa']

tax = {}
with open(taxonomy_file) as file:
    for line in file:
        gtdb_acc, taxonomy = line.rstrip().split('\t')
        included = any([ taxon in taxonomy for taxon in taxa ])
        if included:
            acc = gtdb_acc[3:]
            tax[acc] = taxonomy

with open(output_file, 'w') as output:
    for dirpath, dirnames, filenames in os.walk(assembly_dir):
        for filename in filenames:
            acc_found = None
            for acc in tax:
                if acc in filename:
                    acc_found = acc
                    break
            if acc_found:
                taxonomy = tax.pop(acc_found)
                full_path = os.path.join(dirpath, filename)
                with gzip.open(full_path, 'rb') as file:
                    for line in file:
                        data = line.decode('utf-8')
                        if data.startswith('>'):
                            data = f"{data.split()[0]} {acc_found} {taxonomy}\n"
                        output.write(data)
