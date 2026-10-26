IBS7048 Midterm Problem 5 - Host DNA depletion
================
hjji-BAGEL
2026-10-25

# 5c. Microbial reads, reads without Taxon_X, and fold change

``` r
options(scipen = 999)   # print numbers without scientific notation

# Given values (made up; Taxon_X fraction is among microbial reads)
condition        <- c("Untreated", "Treated")
total_reads      <- c(20000000, 20000000)
host_fraction    <- c(0.99, 0.90)
taxon_X_fraction <- c(0.10, 0.91)

# Microbial reads = total reads * non-host fraction
microbial_reads <- total_reads * (1 - host_fraction)

# Microbial reads that are not Taxon_X
reads_without_X <- microbial_reads * (1 - taxon_X_fraction)

result <- data.frame(condition, total_reads, host_fraction, taxon_X_fraction,
                     microbial_reads, reads_without_X)
result
```

    ##   condition total_reads host_fraction taxon_X_fraction microbial_reads
    ## 1 Untreated    20000000          0.99             0.10          200000
    ## 2   Treated    20000000          0.90             0.91         2000000
    ##   reads_without_X
    ## 1          180000
    ## 2          180000

``` r
# Treated / Untreated fold change for each measure
fold_change <- c(microbial_reads = microbial_reads[2] / microbial_reads[1],
                 reads_without_X = reads_without_X[2] / reads_without_X[1])
fold_change
```

    ## microbial_reads reads_without_X 
    ##              10               1

## Verification

``` r
# 1) Host + Taxon_X + other microbial reads add up to the total reads
host_reads   <- total_reads * host_fraction
taxon_X_reads <- microbial_reads * taxon_X_fraction
all.equal(host_reads + taxon_X_reads + reads_without_X, total_reads)
```

    ## [1] TRUE

``` r
# 2) Same values as the hand calculation
all.equal(microbial_reads, c(20e6 * 0.01, 20e6 * 0.10))
```

    ## [1] TRUE

``` r
all.equal(reads_without_X, c(20e6 * 0.01 * 0.90, 20e6 * 0.10 * 0.09))
```

    ## [1] TRUE

``` r
# 3) Fold changes: 10 for microbial reads, 1 for reads without Taxon_X
all.equal(unname(fold_change), c(10, 1))
```

    ## [1] TRUE

# 5d, 5e. Code excerpts from the exam (not run)

``` r
# Identify contaminants with decontam (combined = frequency + prevalence)
isContaminant(ps, method = "combined", neg = "is.neg",
              conc = "DNA_bac_ng_uL", threshold = 0.1)
```

``` r
# Rarefaction curves of each sample, coloured by treatment
phyloseq.extended::ggrare(ps, step = 1000,
                          color = "treatment", label = "subject_id", se = TRUE)
```
