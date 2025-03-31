SELECT
    accession, lat_parsed, lon_parsed,
    library_strategy IN ('RNA-Seq') OR library_source IN ('METATRANSCRIPTOMIC','TRANSCRIPTOMIC') OR library_selection IN ('cDNA','RT-PCR','Oligo-dT') AS is_transcriptome
FROM metadata WHERE passed AND lat_parsed IS NOT NULL AND lon_parsed IS NOT NULL AND
(
    name IN ('marine metagenome', 'seawater metagenome', 'estuary metagenome', 'lagoon metagenome', 'marine plankton metagenome', 'marsh metagenome') OR 

    env_broad_scale LIKE 'ocean%biome%' OR
    env_broad_scale LIKE 'marine%biome%' OR
    env_broad_scale LIKE 'estuarine%biome%' OR
    env_broad_scale LIKE 'ocean%water%' OR
    env_broad_scale LIKE 'saline%water%' OR
    env_broad_scale LIKE 'estuarine%water%' OR
    env_broad_scale LIKE 'sea%water%' OR
    env_broad_scale LIKE '%pelagic%' OR
    env_broad_scale LIKE 'sea' OR
    env_broad_scale LIKE 'marine' OR
    env_broad_scale LIKE 'ocean' OR
    env_broad_scale LIKE 'estuary' OR

    env_local_scale LIKE 'ocean%biome%' OR
    env_local_scale LIKE 'marine%biome%' OR
    env_local_scale LIKE 'estuarine%biome%' OR
    env_local_scale LIKE 'ocean%water%' OR
    env_local_scale LIKE 'saline%water%' OR
    env_local_scale LIKE 'estuarine%water%' OR
    env_local_scale LIKE 'sea%water%' OR
    env_local_scale LIKE '%pelagic%' OR
    env_local_scale LIKE 'sea' OR
    env_local_scale LIKE 'marine' OR
    env_local_scale LIKE 'ocean' OR
    env_local_scale LIKE 'estuary' OR

    env_medium LIKE 'ocean%biome%' OR
    env_medium LIKE 'marine%biome%' OR
    env_medium LIKE 'estuarine%biome%' OR
    env_medium LIKE 'ocean%water%' OR
    env_medium LIKE 'saline%water%' OR
    env_medium LIKE 'estuarine%water%' OR
    env_medium LIKE 'sea%water%' OR
    env_medium LIKE '%pelagic%' OR
    env_medium LIKE 'sea' OR
    env_medium LIKE 'marine' OR
    env_medium LIKE 'ocean' OR
    env_medium LIKE 'estuary'

)
