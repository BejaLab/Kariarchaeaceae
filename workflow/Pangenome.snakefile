
rule superpang_scaffolds:
    input:
        txt = "analysis/pangenome/superpang/core_origin.txt",
        fna = "analysis/genomes/{genome}.fna",
        fai = "analysis/genomes/{genome}.fna.fai"
    output:
        "analysis/pangenome/genomes/{genome}.fna"
    conda:
        "envs/kits.yaml"
    shell:
        "grep ^{wildcards.genome} {input.txt} | cut -f2 | xargs seqkit faidx {input.fna} > {output}"

rule prodigal_pangenome_scaffolds:
    input:
        "analysis/pangenome/genomes/{genome}.fna"
    output:
        faa = "analysis/pangenome/genomes/{genome}.faa",
        cds = "analysis/pangenome/genomes/{genome}.cds",
        gff = "analysis/pangenome/genomes/{genome}.gff"
    shadow:
        "minimal"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal -i {input} -a {output.faa} -d {output.cds} -f gff -o {output.gff} -p single"

rule superpang:
    input:
        fasta = expand("genomes/{genome}.fna", genome = species_genomes),
        checkm = "analysis/checkm.txt"
    output:
        expand("analysis/pangenome/superpang/{file}", file = [ "assembly.fasta", "assembly.info.tsv" ]),
        expand("analysis/pangenome/superpang/{file}", file = [ "graph.fastg", "graph.NBP2origins.csv" ]),
        expand("analysis/pangenome/superpang/{file}", file = [ "NBP2origins.tsv", "NBPs.accessory.fasta", "NBPs.core.fasta", "NBPs.fasta" ])
    log:
        "analysis/pangenome/superpang.log"
    params:
        ident = 0.8,
        k = 101,
        gap = 10
    conda:
        "envs/superpang.yaml"
    threads:
        20
    shell:
        "SuperPang.py -f {input.fasta} -q {input.checkm} -o $(dirname {output[0]}) -t {threads} -k {params.k} -i {params.ident} -b {params.ident} --gap-size {params.gap} --force-overwrite &> {log}"

rule superpang_core_origin:
    input:
        info = "analysis/pangenome/superpang/assembly.info.tsv",
        origins = "analysis/pangenome/superpang/NBP2origins.tsv"
    output:
        "analysis/pangenome/superpang/core_origin.txt"
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/core_origin.py"

rule superpang_core_fasta:
    input:
        info = "analysis/pangenome/superpang/assembly.info.tsv",
        fasta = "analysis/pangenome/superpang/assembly.fasta"
    output:
        "analysis/pangenome/superpang/assembly.core.fna"
    conda:
        "envs/kits.yaml"
    shell:
        "csvgrep -t -c regions -m core {input.info} | csvcut -c id | sed 1d | seqkit grep -f- {input.fasta} -o {output}"

rule superpang_makeblastdb:
    input:
        "analysis/pangenome/superpang/assembly.core.fna"
    output:
        "analysis/pangenome/superpang/assembly.core.fna.ndb"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule copy_assembly_core:
    input:
        "analysis/pangenome/superpang/assembly.core.fna"
    output:
        "analysis/pangenome/pgap/superpang.fna"
    shell:
        "cp {input} {output}"

rule copy_submol:
    input:
        "metadata/pgap.yaml"
    output:
        "analysis/pangenome/pgap/{genome}_submol.yaml"
    shell:
        "cp {input} {output}"

rule make_yaml:
    input:
        fasta  = "analysis/pangenome/pgap/{genome}.fna",
        submol = "analysis/pangenome/pgap/{genome}_submol.yaml"
    output:
        "analysis/pangenome/pgap/{genome}.yaml"
    wildcard_constraints:
        genome = '.*(?<!_submol)'
    run:
        sets = dict(
            fasta  = { "class": "File", "location": path.basename(input.fasta)  },
            submol = { "class": "File", "location": path.basename(input.submol) }
        )
        with open(output[0], 'w') as fd:
            yaml.dump(sets, fd)

rule pgap:
    input:
        yaml = "analysis/pangenome/pgap/{genome}.yaml",
        submol = "analysis/pangenome/pgap/{genome}_submol.yaml"
    output:
        gbk = "analysis/pangenome/pgap/{genome}/annot.gbk"
    params:
        dirname = directory("analysis/pangenome/pgap/{genome}/")
    threads:
        8
    shell:
        """
        rm -fr {params.dirname}
        pgap.py --no-internet --no-self-update --debug --ignore-all-errors --report-usage-true --cpu {threads} -o {params.dirname} {input.yaml}
        """

rule checkm_pangenome:
    input:
        data = "analysis/checkm_database",
        fasta = "analysis/pangenome/superpang/assembly.core.fna"
    output:
        directory("analysis/pangenome/checkm")
    params:
        domain = "Archaea"
    conda:
        "envs/checkm.yaml"
    threads:
        20
    shell:
        "CHECKM_DATA_PATH={input.data} checkm taxonomy_wf domain {params.domain} $(dirname {input.fasta}) {output} -t {threads}"

rule non_singleton:
    input:
        fasta = "analysis/pangenome/superpang/assembly.fasta",
        tsv = "analysis/pangenome/superpang/assembly.info.tsv"
    output:
        "analysis/pangenome/superpang/assembly.non-singleton.fna"
    conda:
        "envs/kits.yaml"
    params:
        regex = 'core|aux'
    shell:
        "grep -E '{params.regex}' {input.tsv} | cut -f1 | seqkit grep -f- {input.fasta} > {output}"

rule minimap2_index:
    input:
        "analysis/pangenome/superpang/assembly.core.fna"
    output:
        "analysis/pangenome/superpang/assembly.core.fna.mmi"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -d {output} {input}"

rule map_logan:
    input:
        index = "analysis/pangenome/superpang/assembly.core.fna.mmi",
        contigs = "databases/logan/marine_metagenome/{sra}.contigs.fa.zst"
    output:
        "analysis/pangenome/logan/minimap2/{sra}.bam"
    params:
        SM = "LOGAN"
    conda:
        "envs/minimap2.yaml"
    threads:
        3
    shell:
        "zstd -fcd {input.contigs} | minimap2 {input.index} - --sam-hit-only -a -t {threads} | samtools addreplacerg -r ID:{wildcards.sra} -r SM:{params.SM} - | samtools sort -o {output} -@{threads}"

rule merge_logan:
    input:
        expand("analysis/pangenome/logan/minimap2/{sra}.bam", sra = logan_sra)
    output:
        "analysis/pangenome/logan/merge.bam"
    params:
        min_cov = 1000
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/bam_merge.py"

rule rename_fasta:
    input:
        "analysis/pangenome/superpang/assembly.core.fna"
    output:
        "analysis/pangenome/superpang/assembly.core_renamed.fna"
    shell:
        "sed s/length=/length_/ < {input} > {output}"

#rule rename_bam:
#    input:
#        "analysis/pangenome/logan/merge.bam"
#    output:
#        "analysis/pangenome/logan/merge_renamed.bam"
#    conda:
#        "envs/samtools.yaml"
#    shell:
#        "samtools view -h {input} | sed s/length=/length_/ | samtools view -b -o {output}"

rule freebayes:
    input:
        fasta = "analysis/pangenome/superpang/assembly.core.fna",
        bam = "analysis/pangenome/logan/merge.bam"
    output:
        "analysis/pangenome/logan/freebayes.vcf"
    params:
        skip_cov = 1000,
        F = 0.05,
        G = 2
    conda:
        "envs/freebayes.yaml"
    shell:
        "freebayes --skip-coverage {params.skip_cov} --pooled-continuous -F {params.F} -G {params.G} -p 1 -f {input.fasta} {input.bam} > {output}"

rule call_variants:
    input:
        fasta = "analysis/pangenome/superpang/assembly.core_renamed.fna",
        fai = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai",
        dict = "analysis/pangenome/superpang/assembly.core_renamed.dict",
        bam = "analysis/pangenome/logan/merge_renamed.bam",
        bai = "analysis/pangenome/logan/merge_renamed.bam.bai"
    output:
        "analysis/pangenome/logan/mutect2.vcf.gz"
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I {input.bam} -O {output}"

rule cactus_seqfile:
    input:
        pangenome = "analysis/pangenome/superpang/assembly.core.fna",
        genomes = expand("genomes/{genome}.fna", genome = species_genomes)
    output:
        "analysis/pangenome/cactus/genomes.txt"
    params:
        genome_names = species_genomes
    run:
        with open(str(output), 'w') as file:
            file.write(f"pangenome\t{input.pangenome}\n")
            for i, genome in enumerate(params.genome_names):
                file.write(f"{genome}\t{input.genomes[i]}\n")

rule cactus_minigraph:
    input:
        pangenome = "analysis/pangenome/superpang/assembly.core.fna",
        genomes = expand("genomes/{genome}.fna", genome = species_genomes),
        txt = "analysis/pangenome/cactus/genomes.txt"
    output:
        gfa = "analysis/pangenome/cactus/minigraph.gfa",
        job = directory("analysis/pangenome/cactus/minigraph_job")
    log:
        "analysis/pangenome/cactus/minigraph.log"
    conda:
        "envs/cactus.yaml"
    threads:
        20
    shell:
        "cactus-minigraph {output.job} {input.txt} {output.gfa} --reference pangenome --mgCores {threads} --clean never &> {log}"

rule cactus_graphmap:
    input:
        pangenome = "analysis/pangenome/superpang/assembly.core.fna",
        genomes = expand("genomes/{genome}.fna", genome = species_genomes),
        txt = "analysis/pangenome/cactus/genomes.txt",
        gfa = "analysis/pangenome/cactus/minigraph.gfa"
    output:
        job = directory("analysis/pangenome/cactus/graphmap_job"),
        txt = "analysis/pangenome/cactus/graphmap.txt",
        paf = "analysis/pangenome/cactus/graphmap.paf",
        fasta = "analysis/pangenome/cactus/graphmap.fasta"
    log:
        "analysis/pangenome/cactus/graphmap.log"
    conda:
        "envs/cactus.yaml"
    threads:
        20
    shell:
        "cp {input.txt} {output.txt} && cactus-graphmap {output.job} {output.txt} {input.gfa} {output.paf} --reference pangenome --mapCores {threads} --outputFasta {output.fasta} --clean never &> {log}"

rule cactus_align:
    input:
        pangenome = "analysis/pangenome/superpang/assembly.core.fna",
        genomes = expand("genomes/{genome}.fna", genome = species_genomes),
        txt = "analysis/pangenome/cactus/graphmap.txt",
        paf = "analysis/pangenome/cactus/graphmap.paf",
        fasta = "analysis/pangenome/cactus/graphmap.fasta"
    output:
        job = directory("analysis/pangenome/cactus/align_job"),
        hal = "analysis/pangenome/cactus/cactus.hal",
        vg = "analysis/pangenome/cactus/cactus.vg"
    log:
        "analysis/pangenome/cactus/align.log"
    conda:
        "envs/cactus.yaml"
    threads:
        20
    shell:
        "cactus-align {output.job} {input.txt} {input.paf} {output.hal} --reference pangenome --pangenome --outVG --clean never &> {log}"

rule cactus_graphmap_join:
    input:
        pangenome = "analysis/pangenome/superpang/assembly.core.fna",
        genomes = expand("genomes/{genome}.fna", genome = species_genomes),
        txt = "analysis/pangenome/cactus/graphmap.txt",
        vg = "analysis/pangenome/cactus/cactus.vg"
    output:
        job = directory("analysis/pangenome/cactus/graphmap_join_job"),
        gfa = "analysis/pangenome/cactus/cactus.gfa.gz",
        stats = "analysis/pangenome/cactus/cactus.stats.tgz",
        odgi = "analysis/pangenome/cactus/cactus.full.og"
    log:
        "analysis/pangenome/cactus/graphmap_join.log"
    conda:
        "envs/cactus.yaml"
    threads:
        20
    shell:
        "cactus-graphmap-join {output.job} --vg {input.vg} --outDir $(dirname {output.gfa}) --outName cactus --reference pangenome --maxCores {threads} --odgi --clean never &> {log}"

rule hal2maf:
    input:
        "analysis/pangenome/cactus/cactus.hal"
    output:
        job = directory("analysis/pangenome/cactus/hal2maf"),
        maf = "analysis/pangenome/cactus/cactus.maf"
    log:
        "analysis/pangenome/cactus/hal2maf.log"
    threads:
        4
    conda:
        "envs/cactus.yaml"
    shell:
        "cactus-hal2maf {output.job} {input} {output.maf} --dupeMode single --refGenome pangenome --chunkSize 100000000 --noAncestors --batchCores {threads} --clean never &> {log}"

rule maf_to_fasta:
    input:
        "analysis/pangenome/cactus/cactus.maf"
    output:
        fasta = "analysis/pangenome/cactus/cactus.fasta",
        bed = "analysis/pangenome/cactus/cactus.fasta.bed"
    params:
        genomes = species_genomes + [ 'pangenome' ]
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/maf_to_fasta.py"

rule window_clusters:
    input:
        "analysis/pangenome/cactus/cactus.fasta"
    output:
        "analysis/pangenome/cactus/cactus.fasta.clust"
    params:
        max_gaps = 0.5,
        dist_threshold = 0.05,
        window = 1000
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/window_clusters.py"
