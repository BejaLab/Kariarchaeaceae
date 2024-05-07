
from csv import DictReader
import pandas as pd

clades_file = snakemake.input['clades']
gtdb_file = snakemake.input['gtdb']
profile_file = snakemake.input['profile']

profiles_out_file = snakemake.output['profiles']
samples_out_file = snakemake.output['samples']
matches_out_file = snakemake.output['matches']

pident_not_in_gtdb = snakemake.params['pident_not_in_gtdb']

gtdb_df = pd.read_csv(gtdb_file)
ingroup_ids = list(gtdb_df['qseqid'])
clades_df = pd.read_csv(clades_file).query('sseqid in @ingroup_ids or pident >= @pident_not_in_gtdb')
clades_df['stitle'] = (clades_df['sseqid'] + ' ' + clades_df['stitle']).str.split(' // ')
clades_df = clades_df.explode('stitle')
clades_df['gene'] = clades_df['stitle'].str.extract(r'^(\S+)')

genes = clades_df.set_index('gene').T.to_dict()

all_chunks = []
for chunk in pd.read_csv(profile_file, delimiter = '\t', chunksize = 10**6, low_memory = False):
    chunk = chunk.query('`Locus Tag` in @genes')
    if len(chunk.index) > 0:
        all_chunks.append(chunk)

data = (
    pd.concat(all_chunks)
        .merge(clades_df, left_on = "Locus Tag", right_on = "gene")
        .set_index('Scaffold Name')
        .groupby(level = 0).first()
        .copy()
)

data['Length'] = data['Scaffold Length (bp)'].fillna(data['DNA Sequence Length (bp)'])
data['Coverage'] = data['Scaffold Read Depth'].fillna(1) # NB: if no depth is provided we consider the contig as coming from a single read (pair)
data['value'] = data['Coverage'] / data['Length']

profiles = data[[ 'Genome ID', 'value']]
sample_cols = [ 'Sequencing Depth', 'Sequencing Strategy', 'Geographic Location', 'Isolation', 'Habitat', 'Ecosystem', 'Ecosystem Category', 'Ecosystem Subtype', 'Ecosystem Type', 'Specific Ecosystem', 'Biotic Relationships', 'Depth In Meters', 'Latitude', 'Longitude', 'Salinity' ]
samples = data.set_index('Genome ID').groupby(level = 0).first()[sample_cols]
match_cols = list(clades_df.keys()) + [ 'Length', 'Scaffold Length (bp)', 'DNA Sequence Length (bp)', 'Coverage', 'Scaffold Read Depth' ]
matches = data[match_cols]

profiles.to_csv(profiles_out_file)
samples.to_csv(samples_out_file)
matches.to_csv(matches_out_file)
