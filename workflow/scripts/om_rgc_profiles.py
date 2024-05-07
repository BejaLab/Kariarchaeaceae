from csv import DictReader
import pandas as pd
import warnings

clades_file = snakemake.input['clades']
gtdb_file = snakemake.input['gtdb']
profile_metaG_file = snakemake.input['profile_metaG']
profile_metaT_file = snakemake.input['profile_metaT']
metadata_file = snakemake.input['metadata']

profiles_out_file = snakemake.output['profiles']
samples_out_file = snakemake.output['samples']
matches_out_file = snakemake.output['matches']

pident_not_in_gtdb = snakemake.params['pident_not_in_gtdb']

gtdb_df = pd.read_csv(gtdb_file)
ingroup_ids = list(gtdb_df['qseqid'])
clades_df = pd.read_csv(clades_file).query('sseqid in @ingroup_ids or pident >= @pident_not_in_gtdb')
clades_df['om_rgc'] = clades_df['stitle'].str.extract(r'(OM-RGC\S+)', expand = True)
matches = clades_df.query('~om_rgc.isnull()', engine='python').set_index('om_rgc')
records = matches.T.to_dict()

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    metadata = pd.read_excel(metadata_file, sheet_name = 'Table_W1')
samples = metadata
samples['sample'] = samples['PANGAEA sample id'] + '_' + samples['MetaG/MetaT']
samples['Sequencing Strategy'] = samples['MetaG/MetaT'].replace({ 'MetaG': 'Metagenome', 'MetaT': 'Metatranscriptome' })

def read_profiles(file_name, records, metaG_metaT):
    all_chunks = []
    for chunk in pd.read_csv(file_name, delimiter = '\t', chunksize = 10**6):
        chunk = chunk.query('OMRGC_ID in @records')
        if len(chunk.index) > 0:
            all_chunks.append(chunk)
    df = pd.melt(pd.concat(all_chunks), id_vars = "OMRGC_ID", var_name = 'sample', value_name = 'value')
    df['sample'] = df['sample'] + '_' + metaG_metaT
    return(df)

profiles = pd.concat([
    read_profiles(profile_metaG_file, records, 'MetaG'),
    read_profiles(profile_metaT_file, records, 'MetaT')
]).query('value > 0')

profiles.to_csv(profiles_out_file)
samples.to_csv(samples_out_file)
matches.to_csv(matches_out_file)
