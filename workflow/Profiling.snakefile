rule dload_map:
    output:
        directory("analysis/maps/ne_110m_land")
    params:
        url = "https://naciscdn.org/naturalearth/110m/physical/ne_110m_land.zip"
    shell:
        "mkdir -p {output} && wget -qO- '{params.url}' | bsdtar -xvf- -C {output}"

rule om_rgc_profiles_PF01036:
    input:
        clades = "analysis/hmmsearch/OM-RGC_v2_orfs.faa-PF01036.csv",
        profile_metaG = "databases/faa/OM-RGC_v2_gene_profile_metaG.tsv",
        profile_metaT = "databases/faa/OM-RGC_v2_gene_profile_metaT.tsv",
        metadata = "databases/faa/Salazar_et_al_2019_Suppl_Info.xlsx"
    output:
        profiles = "analysis/profiles/OM-RGC_v2_orfs_PF01036_profiles.csv",
        matches = "analysis/profiles/OM-RGC_v2_orfs_PF01036_matches.csv",
        samples = "analysis/profiles/OM-RGC_v2_orfs_PF01036_samples.csv"
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/om_rgc_profiles.py"

rule jgi_img_profiles_PF01036:
    input:
        clades = "analysis/hmmsearch/JGI_IMG_unrestricted.faa-PF01036.csv",
        profile = "databases/faa/JGI_IMG_unrestricted.gene_info.txt"
    output:
        profiles = "analysis/profiles/JGI_IMG_unrestricted_PF01036_profiles.csv",
        matches = "analysis/profiles/JGI_IMG_unrestricted_PF01036_matches.csv",
        samples = "analysis/profiles/JGI_IMG_unrestricted_PF01036_samples.csv"
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/jgi_img_profiles.py"

rule om_rgc_profiles:
    input:
        clades = "analysis/blastp/OM-RGC_v2_orfs.faa_clades.csv",
        profile_metaG = "databases/faa/OM-RGC_v2_gene_profile_metaG.tsv",
        profile_metaT = "databases/faa/OM-RGC_v2_gene_profile_metaT.tsv",
        metadata = "databases/faa/Salazar_et_al_2019_Suppl_Info.xlsx"
    output:
        profiles = "analysis/profiles/OM-RGC_v2_orfs_profiles.csv",
        matches = "analysis/profiles/OM-RGC_v2_orfs_matches.csv",
        samples = "analysis/profiles/OM-RGC_v2_orfs_samples.csv"
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/om_rgc_profiles.py"

rule jgi_img_profiles:
    input:
        clades = "analysis/blastp/JGI_IMG_unrestricted.faa_clades.csv",
        profile = "databases/faa/JGI_IMG_unrestricted.gene_info.txt"
    output:
        profiles = "analysis/profiles/JGI_IMG_unrestricted_profiles.csv",
        matches = "analysis/profiles/JGI_IMG_unrestricted_matches.csv",
        samples = "analysis/profiles/JGI_IMG_unrestricted_samples.csv"
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/jgi_img_profiles.py"

rule plot_map_clades:
    input:
        shape = "analysis/maps/ne_110m_land",
        jgi_img_profiles = "analysis/profiles/JGI_IMG_unrestricted_profiles.csv",
        jgi_img_matches = "analysis/profiles/JGI_IMG_unrestricted_matches.csv",
        jgi_img_samples = "analysis/profiles/JGI_IMG_unrestricted_samples.csv",
        om_rgc_profiles = "analysis/profiles/OM-RGC_v2_orfs_profiles.csv",
        om_rgc_matches = "analysis/profiles/OM-RGC_v2_orfs_matches.csv",
        om_rgc_samples = "analysis/profiles/OM-RGC_v2_orfs_samples.csv"
    output:
        plot = "output/map.svg",
        all_profiles = "output/map_all_profiles.csv",
        om_rgc = "output/map_om_rgc.csv",
        jgi_img = "output/map_jgi_img.csv"
    params:
        clade = "HeimdallR",
        max_depth = 200,
        size_breaks = list(range(5, 21, 5))
    conda:
        "envs/r-map.yaml"
    script:
        "scripts/plot_map.R"

rule plot_map_pfam:
    input:
        shape = "analysis/maps/ne_110m_land",
        jgi_img_profiles = [ "analysis/profiles/JGI_IMG_unrestricted_profiles.csv", "analysis/profiles/JGI_IMG_unrestricted_PF01036_profiles.csv" ],
        jgi_img_matches = [ "analysis/profiles/JGI_IMG_unrestricted_matches.csv", "analysis/profiles/JGI_IMG_unrestricted_PF01036_matches.csv" ],
        jgi_img_samples = [ "analysis/profiles/JGI_IMG_unrestricted_samples.csv", "analysis/profiles/JGI_IMG_unrestricted_PF01036_samples.csv" ],
        om_rgc_profiles = [ "analysis/profiles/OM-RGC_v2_orfs_profiles.csv", "analysis/profiles/OM-RGC_v2_orfs_PF01036_profiles.csv" ],
        om_rgc_matches = [ "analysis/profiles/OM-RGC_v2_orfs_matches.csv", "analysis/profiles/OM-RGC_v2_orfs_PF01036_matches.csv" ],
        om_rgc_samples = [ "analysis/profiles/OM-RGC_v2_orfs_samples.csv", "analysis/profiles/OM-RGC_v2_orfs_PF01036_samples.csv" ]
    output:
        plot = "output/map_pfam.svg",
        all_profiles = "output/map_all_profiles_pfam.csv",
        om_rgc = "output/map_om_rgc_pfam.csv",
        jgi_img = "output/map_jgi_img_pfam.csv"
    params:
        clade = "HeimdallR",
        total_clade = "Bac_rhodopsin",
        max_depth = 200,
        size_breaks = [ 0.05, 0.1, 0.2, 0.4 ]
    conda:
        "envs/r-map.yaml"
    script:
        "scripts/plot_map.R"

rule pfam_dload:
    output:
        "analysis/pfam/{profile}.hmm"
    params:
        url = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{profile}?annotation=hmm"
    shell:
        "wget -O- {params.url} | gzip -cd > {output}"

rule pfam_hmmsearch:
    input:
        hmm = "analysis/pfam/{profile}.hmm",
        pdb = "databases/faa/{database}.pdb"
    output:
        "analysis/hmmsearch/{database}-{profile}.txt"
    params:
        db = "databases/faa/{database}",
    conda:
        "envs/hmmer.yaml"
    threads:
        4
    shell:
        "blastdbcmd -db {params.db} -entry all | hmmsearch --cut_ga -o {output} --cpu {threads} {input.hmm} -"

rule pfam_csv:
    input:
        "analysis/hmmsearch/{database}-{profile}.txt"
    output:
        "analysis/hmmsearch/{database}-{profile}.csv"
    conda:
        "envs/bioformats.yaml"
    script:
        "scripts/hmmsearch_to_csv.py"
