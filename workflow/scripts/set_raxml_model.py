
input_file = str(snakemake.input)
output_file = str(snakemake.output)

prefix = 'best-scoring AA model:'
matrix = None

with open(input_file) as file:
    for line in file:
        line = line.strip()
        if prefix in line:
            model = line.split(prefix)[1]
            matrix = model.split()[0]
            #if 'with empirical base frequencies' in line:
            #    matrix += 'F'

assert matrix, 'Model not found'

with open(input_file) as file, open(output_file, 'w') as out:
    for line in file:
        if line.startswith('Substitution Matrix:'):
            line = f'Substitution Matrix: {matrix}\n'
        out.write(line)
