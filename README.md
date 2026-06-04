# Scripts and infrastructure for the EMERGE temporal paper

This is the public-facing version of the private repository <https://github.com/emerge-bii/temporal_paper>, which was originally copied on 2025-05-30 @ 12:59 EDT, retaining all the original contents except for the exclusion of data files with sharing restrictions (those paths are listed in this repository's .gitignore).

The following people contributed to the original code:

* [Samuel Aroney](https://github.com/AroneyS)
* [Hannah Holland-Moritz](https://github.com/hhollandmoritz)
* [Suzanne Hodgkins](https://github.com/shodgkins)
* [Derek Smith](https://github.com/thiolovinglife)
* [Benjamin Woodcroft](https://github.com/wwood)
* [Sarah Bagby](https://github.com/scbagby)

## Setup

### Create conda environment

``` bash
git clone git@github.com:emerge-bii/temporal_paper.git
cd temporal_paper
./install_dependencies.sh
conda activate temporal_paper_vX
```

If you require packages not installed in the above environment, add them to temporal_paper.yml and update the version number.

### Prepare data

Download from Zenodo at <doi:10.5281/zenodo.13901514> and <doi:10.5281/zenodo.20536110>.
Create symlinks in `data` to the location of these files in your system.

1. Emerge_MAGs_v1
2. SingleM_otu_tables_v4
3. EMERGE_distillate_v9
4. DRAM_distillate_v2
5. Manual_methanogen_calls_v1
6. DRAM_annotations_v4
7. Emerge_metaTs_v5
8. Emerge_metaTs_processed_v6
9. AA_frequencies_v1

Example:

``` bash
cd data
ln -s <path/to/Emerge_MAGs_v1> Emerge_MAGs_v1
```

Large files or folders should also be added here and to `.gitignore`.

### Create your analysis file

Source `setup.R` for all inputs.

### Output results from your analysis into results folder

Separate analysis subfolders for each type.

### Cazyme Scraper

If you plan to use the cazyme scraper, it will need to run in it's own conda environment. Instructions for installation of this environment and use of the script can be found in the [cazyme_scraper repository](https://github.com/hhollandmoritz/cazyme_scraper/tree/09c0d9d1fa99cfc12d075597c77ff3b4f56768ba).
