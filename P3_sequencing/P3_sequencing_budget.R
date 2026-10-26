# =============================================================================
# Problem 3. Number of whole samples and total microbial Gb per platform
# - Budget: USD 600, one platform at a time
# - Non-host bases are treated as microbial
# =============================================================================

# Given values (from the exam)
sequencing <- data.frame(
  platform = c("Illumina", "Oxford Nanopore", "PacBio HiFi"),
  yield_Gb = c(18, 30, 12),
  host_fraction = c(0.20, 0.70, 0.10),
  cost_USD = c(150, 200, 300)
)
budget <- 600

# Number of whole samples within the budget (round down, no partial samples)
sequencing$n_samples <- floor(budget / sequencing$cost_USD)

# Microbial Gb per sample = usable Gb * fraction of non-host bases
sequencing$microbial_Gb_per_sample <- sequencing$yield_Gb * (1 - sequencing$host_fraction)

# Total microbial Gb = number of samples * microbial Gb per sample
sequencing$total_microbial_Gb <- sequencing$n_samples * sequencing$microbial_Gb_per_sample

# Result table and the platform with the largest total microbial yield
sequencing
sequencing$platform[which.max(sequencing$total_microbial_Gb)]

# =============================================================================
# Verification
# =============================================================================

# 1) Money spent must not exceed the budget -> TRUE for all
sequencing$n_samples * sequencing$cost_USD <= budget

# 2) One more sample must go over the budget, so n_samples is the maximum -> TRUE for all
(sequencing$n_samples + 1) * sequencing$cost_USD > budget

# 3) Same values as the hand calculation -> TRUE
all.equal(sequencing$total_microbial_Gb,
          c(4 * 18 * 0.8, 3 * 30 * 0.3, 2 * 12 * 0.9))
