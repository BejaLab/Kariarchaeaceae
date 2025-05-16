# Combine all rhodopsin sets in one fasta file
rule rhodopsins_cat:
    input:
        expand("rhodopsins/{rhodopsin_set}.faa", rhodopsin_set = rhodopsin_sets)
    output:
        "analysis/rhodopsins/proteins.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit grep -vsrp XXX -o {output} {input}"

# Extract TM3 motif residues for rhodopsins in mafft alignment
rule residues:
    input:
        "analysis/rhodopsins/mafft.faa"
    output:
        "analysis/rhodopsins/residues.txt"
    params:
        residues = config['residues']
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/residues.py"

# Align rhodopsins with mafft
rule rhodopsins_mafft:
    input:
        "analysis/rhodopsins/proteins.faa"
    output:
        "analysis/rhodopsins/mafft.faa"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "cat {input} | mafft --auto --reorder - > {output}"

# Trim the mafft alignment
rule rhodopsins_trimal:
    input:
        "analysis/rhodopsins/mafft.faa"
    output:
        "analysis/rhodopsins/trimal.faa"
    params:
        gt = 0.9
    conda:
        "envs/phylophlan.yaml"
    shell:
        "trimal -gt {params.gt} -in {input} -out {output}"

# Do SplitsTree network
# NB: SplitsTree not available from conda
rule splitstree_run:
    input:
        fasta = "analysis/rhodopsins/trimal.faa",
        nexus = "metadata/splitstree.nex"
    output:
        "output/Rhodopsins_NeighborNet_network.nex"
    shell:
        "xvfb-run -a SplitsTreeCMD -x 'IMPORT FILE={input.fasta} DATATYPE=PROTEIN; EXECUTE FILE={input.nexus}; UPDATE; SAVE FILE={output}; QUIT;'"

# Plot the SplitsTree network
rule plot_rhodopsins:
    input:
        residues = "analysis/rhodopsins/residues.txt",
        network = "output/Rhodopsins_NeighborNet_network.nex",
        metadata = "metadata/rhodopsins.xlsx"
    output:
        "output/Rhodopsins_NeighborNet_network.svg"
    conda:
        "envs/r.yaml"
    script:
        "scripts/plot_protein_network.R"

# Cluster rhodopsins with cd-hit
rule cluster_rhodopsins:
    input:
        "analysis/rhodopsins/proteins.faa"
    output:
        "analysis/rhodopsins/proteins.cdhit"
    params:
        c = 0.9
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -c {params.c} -d 0 -i {input} -o {output}"

# Search a nucleotide database with clustered rhodopsins
rule tblastn_fna:
    input:
        query = "analysis/rhodopsins/proteins.cdhit",
        ndb = "databases/fna/{database}.ndb"
    output:
        "analysis/tblastn/{database}.txt"
    params:
        db = "databases/fna/{database}",
        evalue = 1e-5
    conda:
        "envs/blast.yaml"
    threads:
        40
    shell:
        "tblastn -query {input.query} -db {params.db} -outfmt 6 -out {output} -max_target_seqs 1000000000 -num_threads {threads} -evalue {params.evalue}"

# Search a protein database with clustered rhodopsins
rule blastp_database:
    input:
        query = "rhodopsins/{family}.faa",
        pdb = "databases/faa/{database}.pdb"
    output:
        "analysis/blastp/{database}-{family}.txt"
    params:
        db = "databases/faa/{database}",
        evalue = 1e-10,
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qseq', 'sseq', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' ]
    conda:
        "envs/blast.yaml"
    threads:
        20
    shell:
        "blastp -query {input.query} -db {params.db} -outfmt '6 {params.cols}' -out {output} -max_target_seqs 1000000000 -num_threads {threads} -evalue {params.evalue}"

# Search genomes using particular rhodopsin queries 
rule blastp_genomes:
    input:
        query = "rhodopsins/{family}.faa",
        pdb = "analysis/genomes/all_genomes.faa.pdb"
    output:
        "analysis/blastp/all_genomes-{family}.txt"
    params:
        db = "analysis/genomes/all_genomes.faa",
        evalue = 1e-10,
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qseq', 'sseq', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' ]
    conda:
        "envs/blast.yaml"
    threads:
        20
    shell:
        "blastp -query {input.query} -db {params.db} -outfmt '6 {params.cols}' -out {output} -max_target_seqs 1000000000 -num_threads {threads} -evalue {params.evalue}"

# Extract sequences for family-assigned hits
rule blastp_faa:
    input:
        "analysis/blastp/{database}_merged.csv"
    output:
        "analysis/blastp/{database}_merged.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "csvcut -c sseqid,stitle,sseq {input} | csvformat -T | sed 1d | sed 's/\\t/ /' | seqkit tab2fx | seqkit seq -g -o {output}"

# mmseqs easy-taxonomy for taxonomic classification of rhodopsins
rule easy_taxonomy:
    input:
        db_tax = "analysis/rhodopsin_cds/GTDB/coding_taxonomy",
        db = "analysis/rhodopsin_cds/GTDB/coding",
        query = "analysis/blastp/{database}_merged.faa"
    output:
        "analysis/blastp/{database}_merged_lca.tsv"
    params:
        lca_mode = 4,
        prefix = "analysis/blastp/{database}_merged"
    conda:
        "envs/kits.yaml"
    threads:
        32
    shell:
        "mmseqs easy-taxonomy {input.query} {input.db} {params.prefix} tmp --lca-mode {params.lca_mode} --tax-lineage 1 --threads {threads}"

# Assign family to rhodopsin hits in particular database
# based on blastp searches and filtering criteria
rule blastp_merge_hits:
    input:
        outfmt6 = expand("analysis/blastp/{{database}}-{family}.txt", family = rhodopsin_sets),
        metadata = "metadata/rhodopsins.xlsx"
    output:
        "analysis/blastp/{database}_merged.csv"
    params:
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qseq', 'sseq', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' ],
        id = 50,
        score = 90,
        clades = [ 'HeimdallR', 'ACA', 'ACB' ]
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/extract_hits.py"

# Filter the rhodopsin family assignment based on
# taxonomy annotations (retaining only archaeal sequences)
rule blastp_filter_taxa:
    input:
        csv = "analysis/blastp/{database}_merged.csv",
        lca = "analysis/blastp/{database}_merged_lca.tsv"
    output:
        "analysis/blastp/{database}_clades.csv"
    params:
        taxon = "d_Archaea"
    conda:
        "envs/kits.yaml"
    shell:
        "grep {params.taxon} {input.lca} | cut -f1 | csvgrep -c sseqid -f- {input.csv} > {output}"

# Extract family-assigned and taxonomy-filterd rhodopsin sequences
rule blastp_clade_faa:
    input:
        "analysis/blastp/{database}_clades.csv"
    output:
        "analysis/blastp/{database}-{family}.faa"
    params:
        clade = lambda w: { 'HeimdallR': 'HeimdallR', 'ACB': 'ACB', 'ACA': 'PR' }[w.family],
        subclade = lambda w: { 'HeimdallR': '', 'ACB': '', 'ACA': 'ACA' }[w.family]
    conda:
        "envs/kits.yaml"
    shell:
        "csvgrep -c clade -m '{params.clade}' {input} | csvgrep -c subclade -m '{params.subclade}' | csvcut -c sseqid,stitle,sseq | csvformat -T | sed 1d | sed 's/\\t/ /' | seqkit tab2fx | seqkit seq -g -o {output}"

# Copy the fasta file containing pre-filtered rhodopsin CDS's
# extracted from GTDB representative genomes
rule copy_gtdb_coding:
    input:
        "databases/rhodopsins/GTDB_coding.fna"
    output:
        "analysis/rhodopsin_cds/GTDB/coding.fna"
    shell:
        "cp {input} {output}"

# Translate the GTDB rhodopsin CDS's
rule gtdb_coding_translate:
    input:
        "analysis/rhodopsin_cds/GTDB/coding.fna"
    output:
        "analysis/rhodopsin_cds/GTDB/coding.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -g {input} | seqkit translate -o {output}"

# Generate files needed for the database of taxonomically-annotated
# rhodopsin sequences
rule gtdb_coding_tax:
    input:
        "analysis/rhodopsin_cds/GTDB/coding.fna"
    output:
        nodes = "analysis/rhodopsin_cds/GTDB/nodes.dmp",
        names = "analysis/rhodopsin_cds/GTDB/names.dmp",
        mapping = "analysis/rhodopsin_cds/GTDB/mapping",
        merged = "analysis/rhodopsin_cds/GTDB/merged.dmp",
        delnodes = "analysis/rhodopsin_cds/GTDB/delnodes.dmp"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/gtdb_tax.py"

# Create the sequence database for the GTDB rhodopsins
rule gtdb_codon_createdb:
    input:
        "analysis/rhodopsin_cds/GTDB/coding.faa"
    output:
        "analysis/rhodopsin_cds/GTDB/coding"
    conda:
        "envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb {input} {output}"

# Create the taxonomic database based on the rhodopsin CDS's
# extracted from GTDB representative genomes
rule gtdb_codon_createtaxdb:
    input:
        db = "analysis/rhodopsin_cds/GTDB/coding",
        nodes = "analysis/rhodopsin_cds/GTDB/nodes.dmp",
        names = "analysis/rhodopsin_cds/GTDB/names.dmp",
        mapping = "analysis/rhodopsin_cds/GTDB/mapping",
        merged = "analysis/rhodopsin_cds/GTDB/merged.dmp",
        delnodes = "analysis/rhodopsin_cds/GTDB/delnodes.dmp"
    output:
        "analysis/rhodopsin_cds/GTDB/coding_taxonomy"
    conda:
        "envs/mmseqs2.yaml"
    shell:
        "mmseqs createtaxdb {input.db} tmp --ncbi-tax-dump $(dirname {input.nodes}) --tax-mapping-file {input.mapping}"
