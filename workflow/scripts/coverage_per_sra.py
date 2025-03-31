metadata = snakemake.params['metadata']
output_file = str(snakemake.output)

for (index, row), input_file in zip(metadata.iterrows(), snakemake.input):
    with open(input_file) as file:
        for line in file:
            chromosome, depth, bases, size, fraction = line.split()
            if chromosome == 'genome' and depth == '0':
                metadata.loc[index, 'cov_fraction'] = 1 - float(fraction)

metadata.to_csv(output_file)
