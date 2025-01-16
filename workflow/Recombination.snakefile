rule recomb_rhodopsin_sets:
    input:
        "rhodopsins/{rhodopsin_set}.faa"
    output:
        "analysis/recomb/{rhodopsin_set}.cdhit"
    params:
        c = 0.9,
        n = 3
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c {params.c} -n {params.n} -d 0"

rule recomb_all_sequences:
    input:
        expand("analysis/recomb/{rhodopsin_set}.cdhit", rhodopsin_set = rhodopsin_sets),
        "analysis/recomb/GTDB-HeimdallR.cdhit"
    output:
        "analysis/recomb/all_sequences.faa"
    shell:
        "cat {input} > {output}"

rule recomb_cdhit:
    input:
        "analysis/recomb/all_sequences.faa"
    output:
        "analysis/recomb/all_sequences.cdhit"
    params:
        c = 0.99
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output} -c {params.c} -d 0"

rule recomb_mafft:
    input:
        "analysis/recomb/GTDB_combined_reps.faa"
    output:
        "analysis/recomb/GTDB_combined_reps.mafft"
    conda:
        "envs/phylophlan.yaml"
    threads:
        10
    shell:
        "mafft --thread {threads} --auto --reorder {input} > {output}"

rule recomb_backtrans:
    input:
        fna = "analysis/recomb/GTDB-HeimdallR-selected.fna",
        aln = "analysis/recomb/GTDB-HeimdallR-selected.mafft"
    output:
        "analysis/recomb/GTDB-HeimdallR-selected.mafft.fna"
    conda:
        "envs/phylophlan.yaml"
    shell:
        "trimal -backtrans {input.fna} -in {input.aln} -out {output} -splitbystopcodon"

rule recomb_trim:
    input:
        "analysis/recomb/GTDB_combined_reps.mafft"
    output:
        txt = "analysis/recomb/GTDB_combined_reps_trimmed_all.txt",
        fas = "analysis/recomb/GTDB_combined_reps_trimmed_all.faa"
    params:
        gt = 0.4
    conda:
        "envs/phylophlan.yaml"
    shell:
        "trimal -gt {params.gt} -in {input} -colnumbering -out {output.fas} > {output.txt}"

rule recomb_drop_short:
    input:
        "analysis/recomb/GTDB_combined_reps_trimmed_all.faa"
    output:
        "analysis/recomb/GTDB_combined_reps_trimmed.faa"
    params:
        m = 200
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -nigm {params.m} {input} | seqkit grep -f- {input} > {output}"

rule trim_cds_aln:
    input:
        fasta = "analysis/recomb/all_sequences.mafft",
        trimal = "analysis/recomb/all_sequences.mafft.trimmed.txt"
    output:
        "analysis/recomb/all_sequences.trimmed_ends.fna"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_aln.py"

rule gard_hyphy:
    input:
        "analysis/recomb/GTDB-HeimdallR-selected.trimmed.faa"
    output:
        best_gard = "analysis/recomb/gard/best-gard",
        fit_bf = "analysis/recomb/gard/best-gard.fit.bf",
        json = "analysis/recomb/gard/GARD.json"
    params:
        max_breakpoints = 5,
        rate_classes = 3,
        rv = "GDD",
        type = "amino-acid"
    threads:
        20
    log:
        "analysis/recomb/gard/gard.log"
    shadow:
        "minimal"
    conda:
        "envs/hyphy.yaml"
    shell:
        """
        hyphy GARD --mode Faster --type {params.type} --rv {params.rv} --rate-classes {params.rate_classes} --max-breakpoints {params.max_breakpoints} --alignment {input} CPU={threads} &> {log}
        mv {input}.best-gard {output.best_gard}
        mv {input}.best-gard.fit.bf {output.fit_bf}
        mv {input}.GARD.json {output.json}
        """

rule geneconv:
    input:
        "analysis/recomb/GTDB_combined_reps_trimmed.faa"
    output:
        "analysis/recomb/GTDB_combined_reps_geneconv.frags"
    params:
        seed = 123,
        prefix = "analysis/recomb/GTDB_combined_reps_trimmed"
    shadow:
        "minimal"
    log:
        "analysis/recomb/GTDB_combined_reps_geneconv.frags.log"
    shell:
        """
        geneconv {input} -nolog -Seed={params.seed} -Gscale=1 -WideCols &> {log}
        mv {params.prefix}.frags {output}
        """
