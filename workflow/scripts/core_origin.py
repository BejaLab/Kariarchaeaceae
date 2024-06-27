from pandas import read_csv, DataFrame

info_file = snakemake.input['info'] # 'analysis/pangenome/superpang/assembly.info.tsv'
origins_file = snakemake.input['origins'] # 'analysis/pangenome/superpang/NBP2origins.tsv'
output_file = str(snakemake.output)

info = read_csv(info_file, delimiter = '\t')
origins = read_csv(origins_file, delimiter = '\t')

nodes = []
for i, row in info.iterrows():
    if 'core' in row['regions']:
        nodes += row['path'].strip('[]').split(',')

contigs = []
genomes = []
for i, row in origins.iterrows():
    node = row['Node'].split('_')[1]
    if node in nodes:
        contigs += row['SourceContigs'].split(',')
        genomes += row['SourceGenomes'].split(',')

df = DataFrame({ 'genome': genomes, 'contig': contigs })
df.drop_duplicates(inplace = True)
df.to_csv(output_file, sep = '\t', header = False, index = False)
