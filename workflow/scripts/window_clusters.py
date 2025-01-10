from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Align import MultipleSeqAlignment
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

input_file = str(snakemake.input)
window = snakemake.params['window']
max_gaps = snakemake.params['max_gaps']
dist_threshold = snakemake.params['dist_threshold']
output_file = str(snakemake.output)

# Load the alignment file
alignment = AlignIO.read(input_file, "fasta")

def calculate_dist(seq1, seq2):
    return sum(s1 != s2 for s1, s2 in zip(seq1, seq2) if '-' not in (s1, s2))/len(seq1)
def calc_matrix(segment):
    names = [ record.id for record in segment ]
    matrix = []
    for i, rec1 in enumerate(segment):
        matrix.append([])
        for j, rec2 in enumerate(segment):
            if i >= j:
                matrix[i].append(calculate_dist(rec1.seq, rec2.seq))
    return DistanceMatrix(names, matrix)

# For each segment, calculate the distance and store in an array or similar data structure
clusters_assignment = {}
# Determine how many segments there are
segments = len(alignment[0].seq) // window

with open(output_file, 'w') as file:
    for i in range(segments):
        start = i*window
        # Extract the segment
        segment = alignment[:, start : start + window]
        # Make a new alignment including only sequences with less than X% gaps
        valid_align = [ record for record in segment if record.seq.count('-') / window < max_gaps ]
        if len(valid_align) > 1:
            # Create a new MultipleSeqAlignment object for calculation
            valid_segment = MultipleSeqAlignment(valid_align)
            # Calculate the distance matrix
            dm = calc_matrix(valid_align)
            # Store the matrix
            dv = squareform(dm)
            # Apply UPGMA hierarchical clustering
            link = linkage(dv, "average")
            # Define clusters depending on the distance threshold
            clusters = fcluster(link, t = dist_threshold, criterion = 'distance')
            for seq, cluster in zip(valid_segment, clusters):
                file.write(f"{start}\t{seq.id}\t{cluster}\n")
