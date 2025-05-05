# Kariarchaeaceae

This repository contains the [Snakemake](https://snakemake.readthedocs.io/) workflow used for the bioinformatic analyses in Tzlil et al (2025) [Structural insights into light harvesting by antenna-containing rhodopsins in marine Asgard archaea](http://doi.org/10.1038/s41564-025-02016-5) Nat Microbiol.

The workflow was run using Snakemake v. 8.18.2. Most of the dependencies are managed with conda (`--use-conda`). SplitsTreeCMD v. 4 can be downloaded [here](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).

All of the output files are specified in rule `all`. For the data and their description refer to [Structural insights into light harvesting by antenna-containing rhodopsins in marine Asgard archaea. Bioinformatic Data](https://doi.org/10.6084/m9.figshare.26906593).
