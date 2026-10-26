IBS7048 Midterm Problem 6 - Classroom microbiome score and asthma
symptoms
================
hjji-BAGEL
2026-10-26

# 6b. Effect of a 0.20 to 0.60 increase in the microbiome score

``` r
# Given values (per one-unit increase in the microbiome score)
score <- c(0.20, 0.60)
OR    <- c(estimate = 0.47, lower = 0.22, upper = 0.97)
CASI  <- c(estimate = -0.92, lower = -1.51, upper = -0.33)
p0    <- 0.40

# Size of the change in score
delta <- diff(score)
delta
```

    ## [1] 0.4

``` r
# Log-odds scale is linear in the score: OR for a change of delta = OR^delta
OR_delta <- OR^delta
round(OR_delta, 3)
```

    ## estimate    lower    upper 
    ##    0.739    0.546    0.988

``` r
# CASI is linear in the score: difference for a change of delta = coefficient * delta
CASI_delta <- CASI * delta
round(CASI_delta, 3)
```

    ## estimate    lower    upper 
    ##   -0.368   -0.604   -0.132

``` r
# New probability for a child with initial probability p0
# (covariates and random effects held fixed): probability -> odds -> x OR -> probability
odds0 <- p0 / (1 - p0)
odds1 <- odds0 * OR_delta
p1 <- odds1 / (1 + odds1)
round(p1, 3)
```

    ## estimate    lower    upper 
    ##    0.330    0.267    0.397

``` r
# Absolute and relative change in the probability (point estimate)
c(absolute_change = unname(p1["estimate"] - p0),
  relative_change = unname((p1["estimate"] - p0) / p0))
```

    ## absolute_change relative_change 
    ##      -0.0698428      -0.1746070

## Verification

``` r
# 1) Same OR as working on the log scale: exp(delta * log(OR))
all.equal(OR_delta, exp(delta * log(OR)))
```

    ## [1] TRUE

``` r
# 2) Back-transforming the new probability gives the expected odds ratio
all.equal(unname((p1 / (1 - p1)) / odds0), unname(OR_delta))
```

    ## [1] TRUE

``` r
# 3) Point estimates lie inside their confidence intervals
c(OR   = OR_delta[["lower"]] < OR_delta[["estimate"]] & OR_delta[["estimate"]] < OR_delta[["upper"]],
  CASI = CASI_delta[["lower"]] < CASI_delta[["estimate"]] & CASI_delta[["estimate"]] < CASI_delta[["upper"]])
```

    ##   OR CASI 
    ## TRUE TRUE
