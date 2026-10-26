# IBS7048 Advanced Metagenomics - 2026 Fall Midterm

Code and outputs for the take-home midterm exam (Problems 1-6).
Each problem has its own folder. R Markdown files (`.Rmd`) are knitted to GitHub
Markdown (`.md`), so the code and its output can be read directly on GitHub.

## Repository structure

```
.
├── P1_bash/                      Problem 1 (UNIX)
│   ├── make_dummy_fastq.sh       dummy FASTQ generator given in the exam
│   ├── sort_fastq_bash.sh        sorts raw/*.fastq into sorted_bash/<condition>/<read>/
│   └── check_placement.sh        checks that every file is in the folder matching its name
├── P2_R/                         Problem 2 (same sorting with R functions)
│   ├── sort_fastq_R.R
│   ├── P2_sort_R.Rmd
│   └── P2_sort_R.md
├── P3_sequencing/                Problem 3 (samples and microbial Gb within USD 600)
│   ├── P3_sequencing_budget.R
│   ├── P3_sequencing_budget.Rmd
│   └── P3_sequencing_budget.md
├── P4_dada2/                     Problem 4 (DADA2, chimera removal, NJ and ML trees)
│   ├── P4_dada2.Rmd
│   ├── P4_dada2.md
│   └── P4_output/                quality profiles, read tracking tables, top 5 ASVs, Newick trees
├── P5_host_depletion/            Problem 5c (microbial reads and fold change)
│   ├── P5_host_depletion.Rmd
│   └── P5_host_depletion.md
└── P6_classroom_microbiome/      Problem 6b (OR, CASI and probability for a 0.4 score change)
    ├── P6_classroom_microbiome.Rmd
    └── P6_classroom_microbiome.md
```

The dummy FASTQ files are not uploaded. They are created by the script in `P1_bash/`.

## Requirements

- Bash (macOS, Linux or WSL on Windows)
- R (>= 4.3) with `rmarkdown`
- Problem 4: `dada2`, `DECIPHER`, `Biostrings` (Bioconductor) and `phangorn`

```r
install.packages(c("BiocManager", "rmarkdown", "phangorn"))
BiocManager::install(c("dada2", "DECIPHER", "Biostrings"))
```

Problems 2, 3, 5 and 6 use base R only.

## How to run

All commands are run from the repository root unless stated otherwise.

### Problem 1 (Bash)

```bash
bash P1_bash/make_dummy_fastq.sh          # creates IBS7048_2026/raw (40 files)
cd IBS7048_2026
bash ../P1_bash/sort_fastq_bash.sh        # creates sorted_bash/ and prints file counts
bash ../P1_bash/check_placement.sh        # prints "Misplaced files: 0"
cd ..
```

### Problem 2 (R)

Run after Problem 1. The script must be run with `IBS7048_2026` as the working directory.

```bash
cd IBS7048_2026
R -q --no-save -f ../P2_R/sort_fastq_R.R
```

or knit `P2_R/P2_sort_R.Rmd` in RStudio (it sets the working directory to `IBS7048_2026`).

### Problems 3-6 (R Markdown)

Open each `.Rmd` file in RStudio and click **Knit**, or run:

```r
rmarkdown::render("P3_sequencing/P3_sequencing_budget.Rmd")
rmarkdown::render("P4_dada2/P4_dada2.Rmd")              # uses the example data in the dada2 package
rmarkdown::render("P5_host_depletion/P5_host_depletion.Rmd")
rmarkdown::render("P6_classroom_microbiome/P6_classroom_microbiome.Rmd")
```

Problem 4 uses `sam1F/R` and `sam2F/R` in `system.file("extdata", package = "dada2")`,
so no download is needed. `optim.pml` uses `set.seed(2026)` for reproducibility.
In Problems 5 and 6, only the parts that need calculation (5c and 6b) are included.
