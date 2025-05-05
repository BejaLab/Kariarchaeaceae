# Predict genes in genomes (not filtered)
rule prodigal:
    input:
        "analysis/genomes/{genome}.fna"
    output:
        faa = "analysis/genomes/{genome}.faa",
        cds = "analysis/genomes/{genome}.cds",
        gff = "analysis/genomes/{genome}.gff"
    params:
        proc = lambda w: 'single' if w.genome in species_genomes else 'meta'
    shadow:
        "minimal"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal -i {input} -a {output.faa} -d {output.cds} -f gff -o {output.gff} -p {params.proc}"

# Collect predicted proteins (pangenome-filtered for the target species, all for the rest)
# for phylophlan
rule phylophlan_input:
    input:
        expand("analysis/pangenome/genomes/{genome}.faa", genome = species_genomes),
        expand("analysis/genomes/{genome}.faa", genome = other_genomes)
    output:
        directory("analysis/proteins")
    shell:
        "mkdir {output} && cp {input} {output}/"

# Generate customized phylophaln config based on the
# template
rule phylophlan_config:
    input:
        "metadata/phylophlan.cfg"
    output:
        "analysis/proteins_phylophlan/phylophlan.cfg"
    params:
        model = "PROTCATAUTO",
        trimal = "gappyout",
        mafft = "auto",
        evalue = 1e-10,
        seed = 123,
        n_boot = 1000
    shell:
        "MODEL={params.model} SEED={params.seed} N_BOOT={params.n_boot} TRIMAL={params.trimal} MAFFT={params.mafft} EVALUE={params.evalue} envsubst < {input} > {output}"

# Run phylophlan
rule phylophlan:
    input:
        work_dir = "analysis/proteins",
        cfg = "analysis/proteins_phylophlan/phylophlan.cfg",
        maas = "metadata/phylophlan.tsv"
    output:
        info = "analysis/proteins_phylophlan/RAxML_info.proteins.tre",
        aln = "analysis/proteins_phylophlan/proteins_concatenated.aln",
        msas = directory("analysis/proteins_phylophlan/tmp/msas"),
        tree = "analysis/proteins_phylophlan/RAxML_bipartitions.proteins.tre"
    params:
        diversity = "low",
        min_entries = ceil(len(species_genomes + other_genomes) * 0.6)
    conda:
        "envs/phylophlan.yaml"
    threads:
        workflow.cores
    shell:
        "phylophlan -i {input.work_dir} -t a -f {input.cfg} --diversity {params.diversity} --subsample full -d phylophlan --nproc {threads} --output_folder analysis --verbose --maas {input.maas} --min_num_entries {params.min_entries}"

# Run prottest for pplacer
rule prottest:
    input:
        "analysis/pplacer/pplacer_references.fasta",
    output:
        "analysis/pplacer/pplacer_prottest.txt"
    params:
        flags = [ "-F", "-G" ],
        models = [ "-LG", "-WAG" ]
    conda:
        "envs/prottest3.yaml"
    shell:
        "prottest3 -i {input} -o {output} {params.flags} {params.models}"

# Evaluate the species phylogeny using the prottest model and the
# additional non-maker genes (for pplacer)
rule raxml_evaluate:
    input:
        aln = "analysis/pplacer/pplacer_references.fasta",
        tree = "analysis/proteins_phylophlan/RAxML_bipartitions.proteins.tre",
        prottest = "analysis/pplacer/pplacer_prottest.txt"
    output:
        tree = "analysis/pplacer/RAxML_result.pplacer",
        info = "analysis/pplacer/RAxML_info.pplacer"
    params:
        suffix = "pplacer"
    conda:
        "envs/phylophlan.yaml"
    shell:
        """
        model=$(grep -m1 '^Best model according to BIC' {input.prottest} | cut -f6 -d" " | sed -e 's/+G//' -e 's/+F/F/')
        out_dir=$(realpath $(dirname {output.tree}))
        raxmlHPC-PTHREADS-SSE3 -f e -t {input.tree} -s {input.aln} -m "PROTGAMMA${{model}}" -n {params.suffix} -w "$out_dir"
        """

# Create ref package for pplacer
rule taxit:
    input:
        tree = "analysis/pplacer/RAxML_result.pplacer",
        info = "analysis/pplacer/RAxML_info.pplacer",
        aln = "analysis/pplacer/pplacer_references.fasta"
    output:
        directory("analysis/pplacer/RAxML.refpkg")
    conda:
        "envs/pplacer.yaml"
    shell:
        "taxit create -l locus_tag -P {output} --tree-file {input.tree} --aln-fasta {input.aln} --tree-stats {input.info}"

# Run pplacer
rule pplacer:
    input:
        refpkg = "analysis/pplacer/RAxML.refpkg",
        references = "analysis/pplacer/pplacer_references.fasta",
        queries = "analysis/pplacer/pplacer_queries.fasta"
    output:
        "analysis/pplacer/pplacer.jplace"
    conda:
        "envs/pplacer.yaml"
    shell:
        "pplacer -o {output} -c {input.refpkg} -r {input.references} {input.queries}"

# Run blastp to search homologs of genes on the unbinned scaffold (for pplacer)
rule blastp:
    input:
        query = "analysis/genomes/{scaffold}.faa".format(scaffold = scaffold_to_pplace),
        db = "analysis/genomes/{genome}.faa",
        pdb = "analysis/genomes/{genome}.faa.pdb"
    output:
        "analysis/pplacer/blastp/{genome}.txt"
    conda:
        "envs/blast.yaml"
    shell:
        "blastp -query {input.query} -db {input.db} -outfmt 6 -out {output}"

# Create input for pplacer based on homologs of genes on the unbinned scaffold
rule pplacer_input:
    input:
        fasta = expand("analysis/genomes/{genome}.faa", genome = species_genomes + other_genomes),
        blast = expand("analysis/pplacer/blastp/{genome}.txt", genome = species_genomes + other_genomes),
        query = "analysis/genomes/{scaffold}.faa".format(scaffold = scaffold_to_pplace),
        msas = "analysis/proteins_phylophlan/tmp/msas",
        exclude = "analysis/blastp/all_genomes-HeimdallR.txt"
    output:
        queries = "analysis/pplacer/pplacer_queries.fasta",
        references = "analysis/pplacer/pplacer_references.fasta",
        genes = "analysis/pplacer/pplacer_input.txt"
    params:
        genomes = species_genomes + other_genomes,
        scaffold = scaffold_to_pplace,
        evalue = 1e-15,
        min_genes = 3
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/pplacer_input.py"

# Convert jplace pplacer file to newick
rule gappa:
    input:
        "analysis/pplacer/pplacer.jplace"
    output:
        "analysis/pplacer/pplacer.newick"
    conda:
        "envs/gappa.yaml"
    shell:
        "gappa examine graft --jplace-path {input} --out-dir $(dirname {output})"

# Transfer RAxML support values to the newick after pplacer
rule transfer_supports:
    input:
        "analysis/proteins_phylophlan/RAxML_bipartitions.proteins.tre",
        "analysis/pplacer/pplacer.newick"
    output:
        "analysis/pplacer/pplacer_with_supports.newick"
    conda:
        "envs/dendropy.yaml"
    script:
        "scripts/transfer_supports.py"

# Plot the species tree
rule plot_species_tree:
    input:
        jplace = "analysis/pplacer/pplacer.jplace",
        newick = "analysis/pplacer/pplacer_with_supports.newick",
        metadata = "metadata/genomes.xlsx",
        genes = "analysis/pplacer/pplacer_input.txt"
    output:
        "output/species.svg"
    params:
        genomes = species_genomes
    conda:
        "envs/r.yaml"
    script:
        "scripts/plot_species_tree.R"

# Download checkm reference data
rule checkm_dload:
    output:
        directory("analysis/checkm_database")
    params:
        url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
    shell:
        "mkdir -p {output} && wget -O- {params.url} | tar xvz -C {output}"

# Run checkm for the individual genomes
rule checkm_genomes:
    input:
        data = "analysis/checkm_database",
        fasta = expand("analysis/genomes/{genome}.fna", genome = species_genomes + other_genomes)
    output:
        dir = directory("analysis/checkm"),
        txt = "analysis/checkm.txt"
    params:
        domain = "Archaea"
    conda:
        "envs/checkm.yaml"
    shadow:
        "minimal"
    threads:
        10
    shell:
        "CHECKM_DATA_PATH={input.data} checkm taxonomy_wf domain {params.domain} $(dirname {input.fasta[0]}) {output.dir} -f {output.txt} -t {threads}"
