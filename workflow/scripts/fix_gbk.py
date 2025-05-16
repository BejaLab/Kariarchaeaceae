from Bio import SeqIO
import re
from datetime import datetime

organism = snakemake.params['organism']
locus_prefix = snakemake.params['locus_prefix']

with open(str(snakemake.output), 'w') as file:
    for record in SeqIO.parse(str(snakemake.input), 'genbank'):
        record.id = re.sub('_length=\\d+', '', record.id)
        record.description = f"{organism} core pangenome scaffold {record.id}"

        date = record.annotations['structured_comment']['Genome-Annotation-Data']['Annotation Date']
        date_obj = datetime.strptime(date, '%m/%d/%Y %H:%M:%S')
        record.annotations['date'] = date_obj.strftime('%d-%b-%Y').upper()

        record.annotations['source'] = record.annotations['organism'] = record.features[0].qualifiers['organism'] = organism

        for feature in record.features:
            if 'organism' in feature.qualifiers:
                feature.qualifiers['organism'] = organism
            if 'locus_tag' in feature.qualifiers:
                feature.qualifiers['locus_tag'] = re.sub('tmp_', locus_prefix + '_', feature.qualifiers['locus_tag'][0])

        del record.annotations['comment']
        del record.annotations['references']
        record.dbxrefs = []

        SeqIO.write(record, file, 'genbank')
