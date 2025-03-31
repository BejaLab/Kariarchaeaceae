
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
        gbk = "analysis/pangenome/pgap/{genome}/annot.gbk",
        gff = "analysis/pangenome/pgap/{genome}/annot.gff"
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
        contigs = "databases/logan/contigs/{sra}.contigs.fa.zst"
    output:
        "analysis/pangenome/logan/minimap2/{sra}.bam"
    wildcard_constraints:
        sra = '[A-Z0-9]+'
    conda:
        "envs/minimap2.yaml"
    threads:
        3
    shell:
        "zstd -fcd {input.contigs} | minimap2 {input.index} - --sam-hit-only -a -t {threads} | samtools addreplacerg -r ID:{wildcards.sra} - | samtools sort -o {output}"

rule filter_bam:
    input:
        "analysis/pangenome/logan/minimap2/{sra}.bam"
    output:
        "analysis/pangenome/logan/minimap2/{sra}_filt.bam"
    params:
        max_id = 0.2
    conda:
        "envs/bedtools.yaml"
    shell:
        "samtools view -b -e '[de]<{params.max_id}' {input} > {output}"

rule merge_bam_logan_all:
    input:
        expand("analysis/pangenome/logan/minimap2/{sra}_filt.bam", sra = logan_all.index)
    output:
        "analysis/pangenome/logan/merged_all.bam"
    params:
        quals = 'j',
        sort = True,
        trim = '_length=\\d+',
        platform = 'ILLUMINA',
        sample = 'LOGAN'
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/merge_bam.py"

rule merge_bam_logan:
    input:
        expand("analysis/pangenome/logan/minimap2/{sra}_filt.bam", sra = logan_mg.index)
    output:
        "analysis/pangenome/logan/merged_mg.bam"
    params:
        quals = 'j',
        sort = True,
        trim = '_length=\\d+',
        platform = 'ILLUMINA',
        sample = 'LOGAN'
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/merge_bam.py"

rule rename_fasta:
    input:
        "analysis/pangenome/superpang/assembly.core.fna"
    output:
        "analysis/pangenome/superpang/assembly.core_renamed.fna"
    params:
        trim = '_length=\\d+'
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit replace -p '{params.trim}' -o {output} {input}"

rule metagenomic_coverage:
    input:
        fasta = "analysis/pangenome/superpang/assembly.core_renamed.fna",
        fai = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai",
        dict = "analysis/pangenome/superpang/assembly.core_renamed.dict",
        bam = "analysis/pangenome/logan/merged_mg.bam",
        bai = "analysis/pangenome/logan/merged_mg.bam.bai",
        bed = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai.bed"
    output:
        "analysis/pangenome/logan/depth_of_coverage.csv"
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk DepthOfCoverage -R {input.fasta} -I {input.bam} -L {input.bed} -O {output}"

rule call_variants_all:
    input:
        fasta = "analysis/pangenome/superpang/assembly.core_renamed.fna",
        fai = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai",
        dict = "analysis/pangenome/superpang/assembly.core_renamed.dict",
        bam = "analysis/pangenome/logan/merged_all.bam",
        bai = "analysis/pangenome/logan/merged_all.bam.bai"
    output:
        "analysis/pangenome/logan/mutect2_all.vcf.gz"
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I {input.bam} -O {output}"

rule call_variants:
    input:
        fasta = "analysis/pangenome/superpang/assembly.core_renamed.fna",
        fai = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai",
        dict = "analysis/pangenome/superpang/assembly.core_renamed.dict",
        bam = "analysis/pangenome/logan/merged_mg.bam",
        bai = "analysis/pangenome/logan/merged_mg.bam.bai"
    output:
        "analysis/pangenome/logan/mutect2.vcf.gz"
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I {input.bam} -O {output}"

rule filter_variants_all:
    input:
        vcf = "analysis/pangenome/logan/mutect2_all.vcf.gz",
        fasta = "analysis/pangenome/superpang/assembly.core_renamed.fna",
        fai = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai",
        dict = "analysis/pangenome/superpang/assembly.core_renamed.dict"
    output:
        "analysis/pangenome/logan/mutect2_all_filtered.vcf.gz"
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk FilterMutectCalls -V {input.vcf} -R {input.fasta} -O {output}"

rule filter_variants:
    input:
        vcf = "analysis/pangenome/logan/mutect2.vcf.gz",
        fasta = "analysis/pangenome/superpang/assembly.core_renamed.fna",
        fai = "analysis/pangenome/superpang/assembly.core_renamed.fna.fai",
        dict = "analysis/pangenome/superpang/assembly.core_renamed.dict"
    output:
        "analysis/pangenome/logan/mutect2_filtered.vcf.gz"
    params:
        min_fract = 0.05
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk FilterMutectCalls -V {input.vcf} -R {input.fasta} -O {output} --min-allele-fraction {params.min_fract}"

rule vfc_to_bed:
    input:
        "analysis/pangenome/logan/mutect2_filtered.vcf.gz"
    output:
        "analysis/pangenome/logan/mutect2_filtered.bed"
    params:
        filter_out = 'map_qual|base_qual|contamination|weak_evidence|low_allele_frac|normal_artifact|panel_of_normals',
        max_var_width = 9
    conda:
        "envs/bedtools.yaml"
    shell:
        "gzip -cd {input} | awk 'length($4)<={params.max_var_width}&&$7!~/{params.filter_out}/' | bedtools merge > {output}"

rule plot_variation:
    input:
        bed = "analysis/pangenome/logan/mutect2_filtered.bed",
        cov = "analysis/pangenome/logan/depth_of_coverage.csv",
        gff = "analysis/pangenome/pgap/superpang/annot.gff"
    output:
        scaffold = "output/scaffold.svg",
        genome = "output/genome.svg",
        scaffold_csv = "output/scaffold.csv",
        genome_csv = "output/genome.csv"
    params:
        min_scaf_len = 20000,
        min_cds_len = 100 * 3,
        target = "000140",
        win_size = 200,
        discard_dep_percentile = 0.99,
        filter_out = 'map_qual|base_qual|contamination|weak_evidence|low_allele_frac|normal_artifact|panel_of_normals'
    conda:
        "envs/r-plot_variation.yaml"
    script:
        "scripts/plot_variation.R"

rule coverage:
    input:
        "analysis/pangenome/logan/minimap2/{sra}_filt.bam"
    output:
        "analysis/pangenome/logan/minimap2/{sra}_cov.txt"
    params:
        mapq = 0
    conda:
        "envs/bedtools.yaml"
    shell:
        "samtools view -q {params.mapq} -u {input} | bedtools genomecov -ibam - > {output}"

rule coverage_per_sra:
    input:
        expand("analysis/pangenome/logan/minimap2/{sra}_cov.txt", sra = logan_all.index)
    output:
        "analysis/pangenome/logan/coverage_per_sra.csv"
    params:
        metadata = logan_all
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/coverage_per_sra.py"

rule plot_map_logan:
    input:
        genomes = "metadata/genomes.xlsx",
        shape = "analysis/maps/ne_110m_land",
        logan = "analysis/pangenome/logan/coverage_per_sra.csv"
    output:
        plot = "output/map_logan.svg",
        data = "output/map_logan.csv"
    params:
        min_coverage = 0.005,
        taxon = "g__Kariarchaeum"
    conda:
        "envs/r-map.yaml"
    script:
        "scripts/plot_logan_map.R"
