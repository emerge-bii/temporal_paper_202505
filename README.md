# Setup

## Create conda environment

``` bash
git clone git@github.com:emerge-bii/temporal_paper.git
cd temporal_paper
./install_dependencies.sh
conda activate temporal_paper_vX
```

If you require packages not installed in the above environment, add them to temporal_paper.yml and update the version number.

## Create symlinks in `data` to the location of these files in your system

1. Emerge_MAGs_v1 (globus: 20220420_MAGs)
2. SingleM_otu_tables_v4 (globus: SingleM)
3. EMERGE_distillate_v9 (globus: temporal_paper)
4. DRAM_distillate_v2 (globus: temporal_paper)
5. Manual_methanogen_calls_v1 (globus: temporal_paper)
6. DRAM_annotations_v4 (globus: temporal_paper)
7. Emerge_metaTs_v5 (globus: temporal_paper)
8. Emerge_metaTs_processed_v6 (globus: temporal_paper)
9. AA_frequencies_v1 (globus: temporal_paper)

Example:

``` bash
cd data
ln -s <path/to/Emerge_MAGs_v1> Emerge_MAGs_v1
```

Large files or folders should also be added here and to `.gitignore`.

## Create your analysis file

Source `setup.R` for all inputs.

## Output results from your analysis into results folder

Separate analysis subfolders for each type.

## Cazyme Scraper

If you plan to use the cazyme scraper, it will need to run in it's own conda environment. Instructions for installation of this environment and use of the script can be found in the [cazyme_scraper repository](https://github.com/hhollandmoritz/cazyme_scraper/tree/09c0d9d1fa99cfc12d075597c77ff3b4f56768ba).
