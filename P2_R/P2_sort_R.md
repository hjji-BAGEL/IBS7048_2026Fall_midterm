IBS7048 Midterm Problem 2 - Sorting FASTQ files with R
================
hjji-BAGEL
2026-10-25

The dummy FASTQ files and the Bash result (`sorted_bash`) from Problem 1
must exist in `IBS7048_2026/` at the repository root before knitting
this file.

``` r
basename(getwd())   # working directory is IBS7048_2026
```

    ## [1] "IBS7048_2026"

## Sorting with R functions

``` r
# =============================================================================
# Problem 2. Sort the dummy FASTQ files with R (same organization as Problem 1)
# - Input : raw/*.fastq (40 files made in Problem 1)
# - Output: sorted_R/<condition>/<read>/ (files are copied, raw is kept)
# - Run this script with IBS7048_2026 as the working directory
# =============================================================================

# Remove sorted_R from a previous run so the result is always made from scratch.
# raw/ and sorted_bash/ are not touched.
unlink("sorted_R", recursive = TRUE)

# Loop over the 2 conditions and 2 read directions (2 * 2 = 4 folders),
# the same structure as the Bash script in Problem 1.
for (condition in c("Control", "Treated")) {
  for (read in c("R1", "R2")) {

    # Build the destination path, e.g. "sorted_R/Control/R1"
    dest <- file.path("sorted_R", condition, read)

    # Create the folder (and its parent folder) if it does not exist yet.
    # This works like "mkdir -p" in Bash.
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    # Pick the files that match this condition and read direction.
    # ".*" stands for any text, so it covers the replicate number,
    # e.g. "_Control_.*_R1.fastq" matches S01_Control_1_R1.fastq ... S05_Control_2_R1.fastq
    files <- list.files("raw", pattern = paste0("_", condition, "_.*_", read, "\\.fastq$"))

    # Copy the files (not move), so raw/ keeps all 40 files.
    # file.copy() returns TRUE/FALSE for each file; stop if any copy failed.
    stopifnot(all(file.copy(file.path("raw", files), dest)))
  }
}

# =============================================================================
# Verification: compare the R result (sorted_R) with the Bash result (sorted_bash)
# =============================================================================

# List every file under each folder as a relative path, e.g. "Control/R1/S01_Control_1_R1.fastq"
bash_files <- sort(list.files("sorted_bash", recursive = TRUE))
r_files    <- sort(list.files("sorted_R", recursive = TRUE))

# 1) raw/ should still have 40 files because we only copied
length(list.files("raw"))
```

    ## [1] 40

``` r
# 2) Same folder paths and filenames in both results -> TRUE
identical(bash_files, r_files)
```

    ## [1] TRUE

``` r
# 3) Number of files in each folder, Bash vs R (10 each)
cbind(bash = table(dirname(bash_files)),
      R    = table(dirname(r_files)))
```

    ##            bash  R
    ## Control/R1   10 10
    ## Control/R2   10 10
    ## Treated/R1   10 10
    ## Treated/R2   10 10

``` r
# 4) Same file contents: compare MD5 checksums of every file pair -> TRUE
all(tools::md5sum(file.path("sorted_bash", bash_files)) ==
    tools::md5sum(file.path("sorted_R", r_files)))
```

    ## [1] TRUE
