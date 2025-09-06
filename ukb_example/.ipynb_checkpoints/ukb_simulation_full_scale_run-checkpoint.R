#!/usr/bin/env Rscript
library("data.table")
library("tidyverse")

# GWAS result file downloaded from Neale lab UKB results (https://www.nealelab.is/uk-biobank) with command:
# wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O 50_raw.gwas.imputed_v3.both_sexes.tsv.bgz
# Checksums available at: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?usp=sharing
f <- "/nfs/research/birney/users/saul/raw_data/ukb_gwas_results_nealelab/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz"
f_md5 <- tools::md5sum(f) |> unname()
stopifnot(f_md5 == "1c98d9301533ff1b822e7066fc9c4b64")

df <- fread(f)
df_signif <- df[pval < 5e-8]

#nrows <- nrow(df_signif)
nrows <- 1000
set.seed(1)
make_sim <- function(i, df){
    message(sprintf("Progress: %d/%d, %d%%", i, nrows, i/nrows))
    the_snp <- df[i]
    variant <- the_snp[["variant"]]
    n <- the_snp[["n_complete_samples"]]
    maf <- the_snp[["minor_AF"]]
    beta <- the_snp[["beta"]]
    beta_se <- the_snp[["se"]]
    # var(x) computed from the binomial definition as n * p * (1 - p) with n = 2 and p = maf
    var_x <- 2 * maf * (1 - maf)
    # Since: beta_se = resid_se / sqrt(n * var(x))
    # I can derive the standard deviation of the residuals:
    # beta_se^2 = var(resid) / (n * var(x))
    # var(resid) = beta_se^2 * n * var(x)
    # resid_sd = sqrt(beta_se^2 * n * var(x))
    resid_sd <- sqrt(beta_se^2 * n * var_x)
    e <- rnorm(n = n, mean = 0, sd = resid_sd)
    x <- rbinom(n = n, size = 2, prob = maf)
    y <- beta * x + e
    coef <- lm(y ~ x) |> summary() |> coef()
    pval <- coef["x", "Pr(>|t|)"]
    tstat <- coef["x", "t value"]
    ret <- data.table(variant = variant, tstat_sim = tstat, pval_sim = pval)
    return(ret)
}

sim_res <- lapply(1:nrows, make_sim, df = df_signif) |> rbindlist()
comp <- merge(df_signif, sim_res)
fwrite(comp, "/ukb_height_sim_result.csv.gz")