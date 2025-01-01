from dendropy import Tree

from_file = str(snakemake.input[0])
to_file = str(snakemake.input[1])
output_file = str(snakemake.output)

from_tree = Tree.get_from_path(from_file, schema = "newick", preserve_underscores = True)
to_tree = Tree.get_from_path(to_file, schema = "newick", preserve_underscores = True, taxon_namespace = from_tree.taxon_namespace)

shared = len(list(from_tree.leaf_node_iter()))

supports = {}

def get_bits(bp, shared, inverse = False):
    symbol0 = '.' if inverse else '*'
    symbol1 = '*' if inverse else '.'
    return bp.leafset_as_bitstring(symbol0, symbol1, reverse = True)[0:shared]

for bp, edge in from_tree.bipartition_edge_map.items():
    if edge.head_node.label is not None:
        left  = get_bits(bp, shared, inverse = False)
        right = get_bits(bp, shared, inverse = True)
        supports[left]  = edge.head_node.label
        supports[right] = edge.head_node.label

for bp, edge in to_tree.bipartition_edge_map.items():
    left  = get_bits(bp, shared, inverse = False)
    right = get_bits(bp, shared, inverse = True)
    if left in supports:
        edge.head_node.label = supports[left]
    elif right in supports:
        edge.head_node.label = supports[right]

to_tree.write_to_path(output_file, schema = "newick", unquoted_underscores = True)
