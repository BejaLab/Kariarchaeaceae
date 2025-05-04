rule rhodopsins_cat:
    input:
        expand("rhodopsins/{rhodopsin_set}.faa", rhodopsin_set = rhodopsin_sets)
    output:
        "analysis/rhodopsins/proteins.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit grep -vsrp XXX -o {output} {input}"

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

rule rhodopsins_mafft:
    input:
        "analysis/rhodopsins/proteins.faa"
    output:
        "analysis/rhodopsins/mafft.faa"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "cat {input} | mafft --auto --reorder - > {output}"

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

rule rhodopsins_raxml:
    input:
        "analysis/rhodopsins/trimal.faa"
    output:
        "analysis/rhodopsins/RAxML_info.txt",
        "analysis/rhodopsins/RAxML_bipartitions.txt",
        "analysis/rhodopsins/RAxML_bestTree.txt"
    params:
        model = "PROTGAMMAAUTO",
        seed = 123,
        bootstrap = 1000
    conda:
        "envs/phylophlan.yaml"
    threads:
        20
    shell:
        "raxmlHPC-PTHREADS-SSE3 -f a -p {params.seed} -x {params.seed} -# {params.bootstrap} -m {params.model} -T {threads} -s {input} -n txt -w $(dirname $(realpath {output}))"
        # "raxmlHPC-PTHREADS-SSE3 -p {params.seed} -m {params.model} -T {threads} -s {input} -n txt -w $(dirname $(realpath {output}))"

rule rhodopsins_plot:
    input:
        # newick = "analysis/rhodopsins/RAxML_bipartitions.txt",
        metadata = "metadata/rhodopsins.xlsx"
    output:
        plot = "output/rhodopsins.svg",
        jtree = "output/rhodopsins.jtree"
    conda:
        "envs/r.yaml"
    script:
        "scripts/plot_protein_tree.R"

# not available from conda
rule splitstree_run:
    input:
        fasta = "analysis/rhodopsins/trimal.faa",
        nexus = "metadata/splitstree.nex"
    output:
        "analysis/rhodopsins/splitstree.nex"
    shell:
        "xvfb-run -a SplitsTreeCMD -x 'IMPORT FILE={input.fasta} DATATYPE=PROTEIN; EXECUTE FILE={input.nexus}; UPDATE; SAVE FILE={output}; QUIT;'"

rule plot_rhodopsins:
    input:
        residues = "analysis/rhodopsins/residues.txt",
        network = "analysis/rhodopsins/splitstree.nex",
        metadata = "metadata/rhodopsins.xlsx"
    output:
        "output/rhodopsin_network.svg"
    conda:
        "envs/r.yaml"
    script:
        "scripts/plot_protein_network.R"

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

rule blastp_gtdb:
    input:
        query = "rhodopsins/{family}.faa",
        pdb = "analysis/rhodopsin_cds/GTDB/coding.faa.pdb"
    output:
        "analysis/blastp/GTDB-{family}.txt"
    params:
        db = "analysis/rhodopsin_cds/GTDB/coding.faa",
        evalue = 1e-10,
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qseq', 'sseq', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' ]
    conda:
        "envs/blast.yaml"
    threads:
        20
    shell:
        "blastp -query {input.query} -db {params.db} -outfmt '6 {params.cols}' -out {output} -max_target_seqs 1000000000 -num_threads {threads} -evalue {params.evalue}"

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

rule blastp_faa:
    input:
        "analysis/blastp/{database}_merged.csv"
    output:
        "analysis/blastp/{database}_merged.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "csvcut -c sseqid,stitle,sseq {input} | csvformat -T | sed 1d | sed 's/\\t/ /' | seqkit tab2fx | seqkit seq -g -o {output}"

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

rule all_clades_csv:
    input:
        expand("analysis/blastp/{database}_clades.csv", database = [ 'JGI_IMG_unrestricted.faa', 'OM-RGC_v2_orfs.faa', 'all_genomes', 'GTDB' ])

rule subclades_poseidoniia:
    input:
        expand("analysis/blastp/{database}-{{clade}}.faa", database = [ 'JGI_IMG_unrestricted.faa', 'OM-RGC_v2_orfs.faa', 'GTDB' ])
    output:
        "analysis/subclades/poseidoniia-{clade}/sequences.fasta"
    params:
        m = 210
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -m {params.m} -o {output} {input}"

rule subclades_heidallR:
    input:
        expand("analysis/blastp/{database}-{clade}.faa", database = [ 'all_genomes', 'JGI_IMG_unrestricted.faa', 'OM-RGC_v2_orfs.faa' ], clade = [ 'HeimdallR' ])
    output:
        "analysis/subclades/kariarchaeum-heimdallR/sequences.fasta"
    params:
        m = 210
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -m {params.m} {input} | seqkit rmdup -s -o {output}"

rule subclades_poseidoniia_cdhit:
    input:
        "analysis/subclades/poseidoniia-{clade}/sequences.fasta"
    output:
        "analysis/subclades/poseidoniia-{clade}/cdhit.fasta"
    params:
        c = 0.9
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c {params.c} -d 0"

rule subclades_kariarchaeum_cdhit:
    input:
        "analysis/subclades/kariarchaeum-heimdallR/sequences.fasta"
    output:
        "analysis/subclades/kariarchaeum-heimdallR/cdhit.fasta"
    params:
        c = 1.0
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c {params.c} -d 0"

rule subclades_mafft:
    input:
        "analysis/subclades/{taxon}-{clade}/cdhit.fasta"
    output:
        "analysis/subclades/{taxon}-{clade}/mafft.fasta"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "mafft --reorder --auto {input} > {output}"

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

rule link_gtdb_coding:
    input:
        "databases/rhodopsins/GTDB_coding.fna"
    output:
        "analysis/rhodopsin_cds/GTDB/coding.fna"
    shell:
        "cp {input} {output}"

rule gtdb_coding_translate:
    input:
        "analysis/rhodopsin_cds/GTDB/coding.fna"
    output:
        "analysis/rhodopsin_cds/GTDB/coding.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -g {input} | seqkit translate -o {output}"

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

rule gtdb_codon_createdb:
    input:
        "analysis/rhodopsin_cds/GTDB/coding.faa"
    output:
        "analysis/rhodopsin_cds/GTDB/coding"
    conda:
        "envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb {input} {output}"

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

rule tblastn_gtdb:
    input:
        query = "rhodopsins/{clade}.faa",
        ndb = "databases/gtdb/gtdb_genomes_reps_r214.fna.ndb"
    output:
        "analysis/tblastn/gtdb_{clade}.txt"
    params:
        db = "databases/gtdb/gtdb_genomes_reps_r214.fna",
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' ],
        max_seq = 1000,
        evalue = 1e-5
    conda:
        "envs/blast.yaml"
    threads:
        40
    shell:
        "tblastn -query {input.query} -db {params.db} -out {output} -outfmt '6 {params.cols}' -num_threads {threads} -max_target_seqs {params.max_seq} -evalue {params.evalue}"

rule tblastn_gtdb_entries:
    input:
        tblastn = "analysis/tblastn/gtdb_{clade}.txt",
        ndb = "databases/gtdb/gtdb_genomes_reps_r214.fna.ndb"
    output:
        "analysis/tblastn/gtdb_{clade}.fna"
    params:
        db = "databases/gtdb/gtdb_genomes_reps_r214.fna"
    conda:
        "envs/blast.yaml"
    shell:
        "cut -f2 {input.tblastn} | sort -u | blastdbcmd -entry_batch - -db {params.db} > {output}"

rule gtdb_matches:
    input:
        blastp = "analysis/blastp/GTDB-{clade}.txt",
        fasta = "analysis/rhodopsin_cds/GTDB/coding.faa",
        ref = "rhodopsins/{clade}.faa"
    output:
        "analysis/recomb/GTDB-{clade}.faa"
    params:
        evalue_col = 13,
        sseqid_col = 2,
        evalue = 1e-20
    conda:
        "envs/kits.yaml"
    shell:
        "awk '${params.evalue_col}<{params.evalue}' {input.blastp} | cut -f{params.sseqid_col} | seqkit grep -f- {input.fasta} | cat {input.ref} - | seqkit rmdup | seqkit rmdup -s -o {output}"

rule gtdb_combine_matches:
    input:
        expand("analysis/recomb/GTDB-{clade}.faa", clade = [ "HeimdallR", "ACA", "ACB" ])
    output:
        "analysis/recomb/GTDB_combined.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit rmdup -o {output} {input}"

rule gtdb_matches_cdhit:
    input:
        "analysis/recomb/GTDB_combined.faa"
    output:
        fasta = "analysis/recomb/GTDB_combined.cdhit",
        clstr = "analysis/recomb/GTDB_combined.cdhit.clstr"
    params:
        c = 0.9,
        n = 2
    conda:
        "envs/cd-hit.yaml"
    threads:
        10
    shell:
        "cd-hit -i {input} -o {output.fasta} -c {params.c} -n {params.n} -d 0 -T {threads}"

rule gtdb_cluster_reps:
    input:
        refs = expand("rhodopsins/{clade}.faa", clade = [ "HeimdallR", "ACA", "ACB" ]),
        clstr = "analysis/recomb/GTDB_combined.cdhit.clstr",
        fasta = "analysis/recomb/GTDB_combined.faa"
    output:
        "analysis/recomb/GTDB_combined_reps.faa"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/cluster_reps.py"

rule gtdb_heimdallr_matches_select:
    input:
        # additional = "analysis/recomb/recomb-additional.faa",
        matches = "analysis/recomb/GTDB-HeimdallR.fna",
        clstr = "analysis/recomb/GTDB-HeimdallR.cdhit.clstr",
        blast = "analysis/blastp/GTDB-HeimdallR.txt"
    output:
        "analysis/recomb/GTDB-HeimdallR-selected.fna"
    params:
        max_ident = 50,
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qseq', 'sseq', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' ]
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/select_matches.py"

rule gtdb_heimdallr_matches_translate:
    input:
        "analysis/recomb/GTDB-HeimdallR-selected.fna"
    output:
        "analysis/recomb/GTDB-HeimdallR-selected.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit translate -o {output} {input}"

rule tblastn_filter_hits:
    input:
        outfmt6 = "analysis/tblastn/{database}.txt",
        metadata = "metadata/rhodopsins.xlsx"
    output:
        "analysis/tblastn/{database}.csv"
    params:
        id = 40,
        length = 180,
        evalue = 1e-10
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/extract_hits.py"

rule ref_genomes_cat:
    input:
        expand("genomes/{genome}.fna", genome = species_genomes),
        "metadata/pplacer.fna"
    output:
        "analysis/ref_genomes/genomes.fna"
    shell:
        "cat {input} > {output}"

rule ref_genomes_makeblastdb:
    input:
        "analysis/ref_genomes/genomes.fna"
    output:
        "analysis/ref_genomes/genomes.fna.ndb"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule ref_genomes_tblastn:
    input:
        query = "rhodopsins/HeimdallR.faa",
        ndb = "analysis/ref_genomes/genomes.fna.ndb",
        db = "analysis/ref_genomes/genomes.fna"
    output:
        "analysis/tblastn/ref_genomes.txt"
    params:
        cols = [ 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'stitle' ],
        evalue = 1e-30
    conda:
        "envs/blast.yaml"
    shell:
        "tblastn -query {input.query} -db {input.db} -out {output} -outfmt '6 {params.cols}' -evalue {params.evalue}"

