#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Simulation of discovery power under different models for loci with the medaka
heart rate architecture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

process read_pheno_covar {
    label "r_tidyverse_datatable"

    input:
        path pheno
        path covar

    output:
        path "pheno_covar.csv.gz"

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")
        library("tidyverse")

        pheno <- fread("${pheno}")
        covar <- fread("${covar}")

        pheno <- pheno[, .(individual = IID, heart_rate_avg_21C, heart_rate_avg_28C, heart_rate_avg_35C)]
        pheno <- melt(
            pheno,
            id.vars = "individual",
            measure.vars = c("heart_rate_avg_21C", "heart_rate_avg_28C", "heart_rate_avg_35C"),
            value.name = "heart_rate",
            variable.name = "temperature"
        )
        pheno[, temperature_cat := str_remove(temperature, "heart_rate_avg_")]
        pheno[, temperature_cont := str_remove(temperature_cat, "C\$") |> as.numeric()]
        heart_rate_trans <- function(x){${params.heart_rate_trans}}
        pheno[, heart_rate_trans := heart_rate_trans(heart_rate)]
        temp_trans <- function(x){${params.temp_trans}}
        pheno[, temperature_trans := temp_trans(temperature_cont)]

        # the kronoecker product creates id of the kind individual:temperature
        pheno[, full_id := sprintf("%s:%s", individual, temperature_cat)]

        covar <- covar[
            , .(
                individual = IID,
                phenotyping_plate_id,
                cross_id,
                chr15_qtl,
                chr21_qtl
            )
        ]

        df <- merge(pheno, covar, by = "individual")
        fwrite(df, "pheno_covar.csv.gz")
        """
}

process get_qtls_and_models {
    label "r_tidyverse_datatable"

    input:
        path qtls
        path formulas

    output:
        path "qtls.csv.gz", emit: qtls
        path "qtl_models.csv.gz", emit: models

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")

        qtls <- fread("${qtls}")[
            !(locus_id %in% c("chr15_qtl", "chr21_qtl"))
        ]
        formulas <- readRDS("${formulas}")

        qtl_models <- merge(
            qtls[, k := ""],
            data.table(model = names(formulas), k = ""),
            allow.cartesian = TRUE,
            by = "k"
        )

        qtls[, k := NULL]
        qtl_models[, k := NULL]

        fwrite(qtls, "qtls.csv.gz")
        fwrite(qtl_models, "qtl_models.csv.gz")
        """
}

process get_formulas {
    label "r_tidyverse_datatable"

    output:
        path "formulas.rds"

    script:
        """
        #!/usr/bin/env Rscript

        formulas <- list(
            model_gxg = formula(
                heart_rate_trans ~
                1 +
                cross_id*temperature_trans +
                phenotyping_plate_id +
                chr15_qtl*temperature_trans +
                chr21_qtl*temperature_trans +
                snp*temperature_trans +
                dominance*temperature_trans +
                chr15_qtl*snp*temperature_trans +
                chr21_qtl*snp*temperature_trans
            ),
            model_no_gxg = formula(
                heart_rate_trans ~
                1 +
                cross_id*temperature_trans +
                phenotyping_plate_id +
                chr15_qtl*temperature_trans +
                chr21_qtl*temperature_trans +
                snp*temperature_trans +
                dominance*temperature_trans +
                chr15_qtl*temperature_trans +
                chr21_qtl*temperature_trans
            )
        )

        saveRDS(formulas, "formulas.rds")
        """
}

process make_freq {
    label "plink2"

    input:
        path vcf

    output:
        path "freq.afreq.zst"

    script:
        """
        plink2 \\
            --threads ${task.cpus} \\
            --memory ${task.memory.getMega()} \\
            --vcf ${vcf} \\
            --min-alleles 2 \\
            --max-alleles 2 \\
            --set-missing-var-ids @_#_\\\$r_\\\$a \\
            --new-id-max-allele-len 10 missing \\
            --out freq \\
            --chr-set ${params.n_chr} \\
            --freq zs
        """
}

process make_pgen {
    label "plink2"

    input:
        path vcf

    output:
        tuple(
            path("pgen.pgen"),
            path("pgen.psam"),
            path("pgen.pvar.zst")
        )

    script:
        """
        plink2 \\
            --threads ${task.cpus} \\
            --memory ${task.memory.getMega()} \\
            --set-all-var-ids @_#_\\\$r_\\\$a \\
            --min-alleles 2 \\
            --max-alleles 2 \\
            --out pgen \\
            --chr-set ${params.n_chr} \\
            --vcf ${vcf} dosage=DS \\
            --make-pgen vzs fill-missing-from-dosage erase-dosage
        """
}

process make_grm {
    label "plink2"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(pgen),
            path(psam),
            path(pvar),
            path(freq)
        )

    output:
        tuple(
            val(meta),
            path("grm_loco_${meta.id}.rel.id"),
            path("grm_loco_${meta.id}.rel.bin")
        )

    script:
        """
        plink2 \
            --pgen $pgen \\
            --psam $psam \\
            --pvar $pvar \\
            --threads ${task.cpus} \\
            --memory ${task.memory.getMega()} \\
            --not-chr "15,21,${meta.chr}" \\
            --chr-set ${params.n_chr} \\
            --read-freq $freq \\
            --out "grm_loco_${meta.id}" \\
            --maf 0.01 \\
            --make-rel bin
        """
}

process get_qtl_matrices {
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(grm_id),
            path(grm_bin),
            path(pheno_covar),
            path(pgen),
            path(psam),
            path(pvar),
            path(formulas)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}.qtl_matrices.rds"),
            optional: true
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")
        library("tidyverse")

        get_snp <- function(var_id){
            pvar <- pgenlibr::NewPvar("${pvar}")
            pgen <- pgenlibr::NewPgen("${pgen}", pvar = pvar)
            ret <- fread("${psam}")[, .(individual = `#IID`)]
            ret[["snp"]] <- pgenlibr::Buf(pgen)
            var_num <- pgenlibr::GetVariantsById(pvar, var_id)
            pgenlibr::Read(pgen, ret[["snp"]], var_num, allele_num = 2L)
            ret[["dominance"]] = as.numeric(ret[["snp"]] == 1)
            return(ret)
        }

        get_design_matrix <- function(the_formula, model_frame){
            X <- model.matrix(the_formula, data = model_frame)
            return(X)
        }

        read_K <- function(grm_id, grm_bin){
            samples_K <- read.table(grm_id, header = FALSE, check.names = FALSE)[,1]
            K <- matrix(
                readBin(grm_bin, what = "numeric", n = length(samples_K) ** 2),
                ncol = length(samples_K)
            )
            colnames(K) <- samples_K
            rownames(K) <- samples_K
            stopifnot(sum(is.na(K)) == 0)
            return(K)
        }

        get_K_ind <- function(df, model_frame){
            # expand the factors in a dummy encoded matrix with n_samples rows and n_groups columns
            # needs to be df and not model_frame because the latter does not contain indivudual
            Z <- model.matrix(~ 0 + individual, data = df)
            # get a block diagonal square matrix with n_samples rows and columns and
            # 1 for samples in the same group, 0 for samples in different groups
            K <- Z %*% t(Z)

            match_v <- match(rownames(model_frame), rownames(K))
            K <- K[match_v, match_v]

            return(K)
        }

        get_K_list <- function(grm_id, grm_bin, model_frame, df){
            K_grm <- read_K(grm_id, grm_bin)
            possible_temps <- df[["temperature_cat"]] |> unique()
            n_temps <- length(possible_temps)
            # reletedness among environments is diagonal
            Kt_indep <- diag(n_temps)
            rownames(Kt_indep) <- possible_temps
            colnames(Kt_indep) <- possible_temps
            # relatedness among environment is full
            Kt_full <- matrix(1, n_temps, n_temps)
            rownames(Kt_full) <- possible_temps
            colnames(Kt_full) <- possible_temps
            # kronecker products of the relationships among environments and genetic relatedness
            # order does not matter since I match samples afterwards (matters for the row and colnames)
            K_grm_indep <- kronecker(K_grm, Kt_indep, make.dimnames = TRUE) # this is used to model s2_gxe
            K_grm_full <- kronecker(K_grm, Kt_full, make.dimnames = TRUE) # this is used to model s2_g
            # s2_g and s2_gxe together create the "compound symmetry" model for GxE where an overall variance s2_g is estimated and a within-environment variance s2_gxe
            # more complex models are possible but probably overkill
            match_v <- match(rownames(model_frame), rownames(K_grm_indep))
            # one match is sufficient, Kt_indep and Kt_full have the same names
            ret <- list(
                K_grm_indep = K_grm_indep[match_v, match_v],
                K_grm_full = K_grm_full[match_v, match_v],
                K_ind = get_K_ind(df, model_frame) # correlation among sample duplicates
            )
            return(ret)
        }

        pheno_covar <- fread("${pheno_covar}")
        formulas <- readRDS("${formulas}")
        the_formula <- formulas[["${meta.model}"]]
        snp <- get_snp("${meta.lead_snp_id}")
        df <- merge(pheno_covar, snp, by = "individual")
        df <- as.data.frame(df)
        rownames(df) <- df[["full_id"]]
        model_frame <- model.frame(the_formula, data = df)
        K_list <- get_K_list("${grm_id}", "${grm_bin}", model_frame, df)
        y <- model.response(model_frame)
        X <- get_design_matrix(the_formula, model_frame)
        res <- list(
            y = y,
            X = X,
            K_list = K_list
        )

        saveRDS(res, "${meta.id}.qtl_matrices.rds")
        """
}

process fit_mixed_model {
    label "r_gaston"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(qtl_matrices)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}.mm_fit.rds"),
        )

    script:
        """
        #!/usr/bin/env Rscript

        qtl_matrices <- readRDS("${qtl_matrices}")
        X <- qtl_matrices[["X"]]
        y <- qtl_matrices[["y"]]
        K_list <- qtl_matrices[["K_list"]]
        fit <- gaston::lmm.aireml(
            Y = y,
            X = X,
            K = K_list,
            verbose = TRUE
        )

        saveRDS(fit, "${meta.id}.mm_fit.rds")
        """
}

process decorrelate_matrices {
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(qtl_matrices),
            path(mixed_model)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}.mm_matrices.rds")
        )

    script:
        """
        #!/usr/bin/env Rscript

        qtl_matrices <- readRDS("${qtl_matrices}")
        fit <- readRDS("${mixed_model}")

        X <- qtl_matrices[["X"]]
        y <- qtl_matrices[["y"]]
        K_list <- qtl_matrices[["K_list"]]
        sigma2 <- fit[["sigma2"]]
        tau <- fit[["tau"]]
        names(tau) <- names(K_list)
        sigma2_tot <- sigma2 + sum(tau)
        var_explained <- lapply(tau, function(x){x / sigma2_tot})
        # residual variance-covariance matrix
        V <- Reduce(
            "+", lapply(1:length(K_list), function(i){tau[[i]] * K_list[[i]]})
        ) + diag(
            sigma2, dim(K_list[[1]])
        )
        L <- t(chol(V)) # Rlang returns the upper Cholesky triangle
        colnames(L) <- colnames(K_list[[1]])
        rownames(L) <- rownames(K_list[[1]])
        # remove covariance structure
        y.mm <- forwardsolve(L, y)
        X.mm <- forwardsolve(L, X)
        names(y.mm) <- names(y)
        rownames(X.mm) <- rownames(X)
        colnames(X.mm) <- colnames(X)
        res <- list(
            X.mm = X.mm,
            y.mm = y.mm,
            # residual independentent SD of the errors
            sigma = sqrt(sigma2)
        )

        saveRDS(res, "${meta.id}.mm_matrices.rds")
        """
}

process fit_lm {
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(mm_matrices_full)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}.lm_fit.rds")
        )

    script:
        """
        #!/usr/bin/env Rscript

        mm_mat_full <- readRDS("${mm_matrices_full}")
        X.mm <- mm_mat_full[["X.mm"]]
        y.mm <- mm_mat_full[["y.mm"]]
        sigma <- mm_mat_full[["sigma"]]
        stopifnot(all.equal(rownames(X.mm), names(y.mm)))
        fit <- lm.fit(y = y.mm, x = X.mm)
        ret <- list(
            locus_id = "${meta.locus_id}",
            model = "${meta.model}",
            fit = fit,
            sigma = sigma
        )

        saveRDS(ret, "${meta.id}.lm_fit.rds")
        """
}

process simulate {
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(lm_fit),
            path(temperature_file),
            val(seed)
        )
        val noise_levels
        val n
        val r2
        val allele_freq

    output:
        tuple(
            val(meta),
            path("*.simulated_data.rds"),
            emit: fits,
            optional: true
        )
        tuple(
            val(meta),
            path("${meta.id}.csv.gz"),
            emit: csv
        )

    script:
        """
        #!/usr/bin/env Rscript

        set.seed(${seed})

        library("data.table")

        add_reciprocal_p <- function(p){
            # add reciprocal allele_freq to a vector of allele_freq
            # e.g. 0.1 -> 0.1, 0.9
            # unique() avoids 0.5 -> 0.5, 0.5
            c(p, 1 - p) |> unique()
        }

        fit_model <- function (model, vars, y) {
            fit <- lm(formula(model[["model"]]), data = vars)
            fit_null <- lm(formula(model[["null_model"]]), data = vars)

            # Likelihood ratio test
            ll_fit <- logLik(fit)
            ll_null <- logLik(fit_null)
            lrt_chisq <- as.numeric(2 * (ll_fit - ll_null))
            df <- attributes(ll_fit)[["df"]] - attributes(ll_null)[["df"]]
            pval <- pchisq(lrt_chisq, df = df, lower.tail = FALSE)

            ret <- list(
                fit = fit,
                fit_null = fit_null,
                lrt_chisq = lrt_chisq,
                df = df,
                pval = pval
            )
            return(ret)
        }

        extract_pvals <- function(model_name, model_fits) {
            ret <- model_fits[[model_name]][c("lrt_chisq", "df", "pval")]
            ret[["model"]] <- model_name
            return(ret)
        }

        generate_correlated_snp <- function(snp1, r2, snp1_p = mean(snp1)/2) {
            # generate a SNP with arbitrary r2 to another SNP

            # probabilities for all possible pairs for 2 correlated bernoulli variables
            # see https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables
            # assuming p == q (same SNP frequencies)
            p <- snp1_p

            a <- sqrt(r2) * (p * (1 - p)) + (1 - p)^2
            prob_bernoulli <- c(
              p00 = a,
              p10 = 1 - p - a,
              p01 = 1 - p - a,
              p11 = a + 2 * p - 1
            )
            stopifnot(all.equal(sum(prob_bernoulli), 1))

            # extend the bernoulli probabilities to the conditional probabilities for
            # a pair of binomials of size 2 given the result of the first binomial

            # first compute the joint probabilities for the pair of binomials
            prob_binomial_joint <- c(
              p00 = prob_bernoulli[["p00"]]^2,
              p22 = prob_bernoulli[["p11"]]^2,
              p02 = prob_bernoulli[["p01"]]^2,
              p20 = prob_bernoulli[["p10"]]^2,
              p01 = prob_bernoulli[["p00"]] * prob_bernoulli[["p01"]] * 2,
              p10 = prob_bernoulli[["p00"]] * prob_bernoulli[["p10"]] * 2,
              p21 = prob_bernoulli[["p11"]] * prob_bernoulli[["p10"]] * 2,
              p12 = prob_bernoulli[["p11"]] * prob_bernoulli[["p01"]] * 2,
              p11 = prob_bernoulli[["p01"]] * prob_bernoulli[["p10"]] * 2 + prob_bernoulli[["p00"]] * prob_bernoulli[["p11"]] * 2
            )
            stopifnot(all.equal(sum(prob_binomial_joint), 1))

            # prior for the binomial
            prior_binom <- c(
              p0 = p^2, p1 = 2 * p * (1 - p), p2 = (1 - p)^2
            )
            stopifnot(all.equal(sum(prior_binom), 1))

            # now finally the conditional of the result of the second binomial given the first
            prob_binomial_conditional <- c(
              p00 = prob_binomial_joint[["p00"]] / prior_binom[["p0"]],
              p01 = prob_binomial_joint[["p01"]] / prior_binom[["p1"]],
              p02 = prob_binomial_joint[["p02"]] / prior_binom[["p2"]],
              p10 = prob_binomial_joint[["p10"]] / prior_binom[["p0"]],
              p11 = prob_binomial_joint[["p11"]] / prior_binom[["p1"]],
              p12 = prob_binomial_joint[["p12"]] / prior_binom[["p2"]],
              p20 = prob_binomial_joint[["p20"]] / prior_binom[["p0"]],
              p21 = prob_binomial_joint[["p21"]] / prior_binom[["p1"]],
              p22 = prob_binomial_joint[["p22"]] / prior_binom[["p2"]]
            )

            sample_snp2 <- function(x, probs){
              # subset the desired conditional
              probs_curr <- probs[sprintf("p%s%s", x, 0:2)]
              sample(0:2, size = 1, prob = probs_curr)
            }

            snp2 <- sapply(snp1, sample_snp2, probs = prob_binomial_conditional)

            return(snp2)
        }

        simulate <- function(i){
            # unpack simulation parameters
            n <- params_df[i, n]
            p <- params_df[i, p]
            r2 <- params_df[i, r2]
            noise_lev <- params_df[i, noise_lev]

            message(sprintf(
                "Simulation %s of %s, n = %s, p = %s, r2 = %s, noise = %s",
                i, i_max, n, p, r2, noise_lev
            ))

            # simulated model variables
            temp_raw <- sample(temperature_vec, size = n, replace = TRUE)
            vars <- list(
                intercept = rep(1, n),
                snp = rbinom(n = n, size = 2, prob = p),
                chr15_qtl = rbinom(n = n, size = 2, prob = chr15_qtl_p),
                chr21_qtl = rbinom(n = n, size = 2, prob = chr21_qtl_p),
                temp = temp_trans(temp_raw)
            )
            vars[["snp:temp"]] <- vars[["snp"]] * vars[["temp"]]
            vars[["dominance"]] <- as.numeric(vars[["snp"]] == 1)
            vars[["dominance:temp"]] <- vars[["dominance"]] * vars[["temp"]]
            vars[["snp:chr15_qtl"]] <- vars[["snp"]] * vars[["chr15_qtl"]]
            vars[["snp:chr21_qtl"]] <- vars[["snp"]] * vars[["chr21_qtl"]]
            vars[["chr15_qtl:temp"]] <- vars[["snp"]] * vars[["chr15_qtl"]]
            vars[["chr21_qtl:temp"]] <- vars[["snp"]] * vars[["chr21_qtl"]]
            vars[["snp:chr15_qtl:temp"]] <- vars[["snp"]] * vars[["temp"]] * vars[["chr15_qtl"]]
            vars[["snp:chr21_qtl:temp"]] <- vars[["snp"]] * vars[["temp"]] * vars[["chr21_qtl"]]
            vars[["chr15_qtl:temp"]] <- vars[["temp"]] * vars[["chr15_qtl"]]
            vars[["chr21_qtl:temp"]] <- vars[["temp"]] * vars[["chr21_qtl"]]

            # residual error of the same variance as in the model fit
            vars[["resid"]] <- rnorm(n, sd = sigma)

            # phenotype vector created by summing all the effects and the
            # residual
            bX <- do.call(
                "cbind",
                lapply(
                    names(vars),
                    function(var_name){vars[[var_name]] * betas[[var_name]]}
                )
            )
            y <- rowSums(bX)

            # add noise to the temperature if needed
            if (noise_lev != 0) {
                # noise added to the raw measure
                noise_sd <- sqrt(var(temp_raw) * noise_lev)
                vars[["temp"]] <- temp_trans(
                    temp_raw + rnorm(n, sd = noise_sd)
                )
            }


            # simulate a snp in ld if needed
            if (!is.na(r2)) {
                vars[["snp"]] <- generate_correlated_snp(
                    vars[["snp"]], r2 = r2, snp1_p = p
                )
                vars[["chr15_qtl"]] <- generate_correlated_snp(
                    vars[["chr15_qtl"]], r2 = r2, snp1_p = chr15_qtl_p
                )
                vars[["chr21_qtl"]] <- generate_correlated_snp(
                    vars[["chr21_qtl"]], r2 = r2, snp1_p = chr15_qtl_p
                )
            }

            # regenerate interactions with the new terms
            vars[["snp:temp"]] <- vars[["snp"]] * vars[["temp"]]
            vars[["dominance"]] <- as.numeric(vars[["snp"]] == 1)
            vars[["dominance:temp"]] <- vars[["dominance"]] * vars[["temp"]]
            vars[["snp:chr15_qtl"]] <- vars[["snp"]] * vars[["chr15_qtl"]]
            vars[["snp:chr21_qtl"]] <- vars[["snp"]] * vars[["chr21_qtl"]]
            vars[["chr15_qtl:temp"]] <- vars[["snp"]] * vars[["chr15_qtl"]]
            vars[["chr21_qtl:temp"]] <- vars[["snp"]] * vars[["chr21_qtl"]]
            vars[["snp:chr15_qtl:temp"]] <- vars[["snp"]] * vars[["temp"]] * vars[["chr15_qtl"]]
            vars[["snp:chr21_qtl:temp"]] <- vars[["snp"]] * vars[["temp"]] * vars[["chr21_qtl"]]
            vars[["chr15_qtl:temp"]] <- vars[["temp"]] * vars[["chr15_qtl"]]
            vars[["chr21_qtl:temp"]] <- vars[["temp"]] * vars[["chr21_qtl"]]

            # fit all the models
            model_fits <- lapply(models, fit_model, vars = vars, y = y)

            # save fit examples for only one rep and a specific n
            if ( ${meta.rep} == 1 & n == 10000 ) {
                saveRDS(
                    model_fits,
                    sprintf(
                        "${meta.id}_p%s_noise_lev%s_rsq%s.simulated_data.rds",
                        p, noise_lev, r2
                    )
                )
            }

            df <- lapply(
                names(model_fits), extract_pvals, model_fits = model_fits
            ) |> rbindlist()
            df[, n := n]
            df[, p := p]
            df[, noise_lev := noise_lev]
            df[, r2 := r2]
            return(df)
        }

        # allele frequencies from wild data
        chr15_qtl_p <- 0.03
        chr21_qtl_p <- 0.005

        # variables to iterate over across simulations
        sim_params <- list(
            n = c(${n.join(', ')}),
            p = c(${allele_freq.join(', ')}) |> add_reciprocal_p(),
            r2 = c(${r2.join(', ')}),
            noise_lev = c(${noise_levels.join(', ')})
        )
        # all the combinations of parameters
        params_df <- expand.grid(sim_params) |> as.data.table()

        # fit on the real data to take effect sizes
        l <- readRDS("${lm_fit}")
        fit_orig <- l[["fit"]]
        sigma <- l[["sigma"]]
        temperature_vec <- fread("${temperature_file}")[
            ${params.temperature_filter}, temperature_celsius
        ]
        # perform same transformation as in the discovery
        temp_trans <- function(x){${params.temp_trans}}

        # linear coefficients from the original fit
        fit_coef <- coef(fit_orig)

        for (
            v in c(
                "chr15_qtl:snp",
                "chr21_qtl:snp",
                "temperature_trans:chr15_qtl:snp",
                "temperature_trans:chr21_qtl:snp"
            )
        ) {
            if (!(v %in% names(fit_coef))) fit_coef[[v]] <- 0
        }
        stopifnot(length(fit_coef) == 66)


        betas <- list(
            intercept = fit_coef[["(Intercept)"]],
            snp = fit_coef[["snp"]],
            dominance = fit_coef[["dominance"]],
            temp = fit_coef[["temperature_trans"]],
            chr15_qtl = fit_coef[["chr15_qtl"]],
            chr21_qtl = fit_coef[["chr21_qtl"]],
            `chr15_qtl:temp` = fit_coef[["temperature_trans:chr15_qtl"]],
            `chr21_qtl:temp` = fit_coef[["temperature_trans:chr21_qtl"]],
            `snp:temp` = fit_coef[["temperature_trans:snp"]],
            `dominance:temp` = fit_coef[["temperature_trans:dominance"]],
            `snp:chr15_qtl` = fit_coef[["chr15_qtl:snp"]],
            `snp:chr21_qtl` = fit_coef[["chr21_qtl:snp"]],
            `snp:chr15_qtl:temp` = fit_coef[["temperature_trans:chr15_qtl:snp"]],
            `snp:chr21_qtl:temp` = fit_coef[["temperature_trans:chr21_qtl:snp"]],
            # just to make the code simpler - the residual does not have an
            # effect of course
            resid = 1
        )

        # replace missing betas with 0
        betas[is.na(betas)] <- 0

        # models to fit to the simulated data - I always test whatever the
        # model is versus an environment and interacting loci-only model. When
        # the environment or the intera loci are
        # not in the model I test against an intercept only model.
        models <- list(
            gxgxe_gxg_gxe_dominance = list(
                model = "y ~ temp*snp*chr15_qtl + temp*snp*chr21_qtl + dominance*temp",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            gxgxe_gxg_gxe = list(
                model = "y ~ temp*snp*chr15_qtl + temp*snp*chr21_qtl",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            gxg_gxe = list(
                model = "y ~ temp*snp + snp*chr15_qtl + snp*chr21_qtl + temp*chr15_qtl + temp*chr21_qtl",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            gxe_and_dominance = list(
                model = "y ~ temp*snp + temp*dominance + temp*chr15_qtl + temp*chr21_qtl",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            gxe = list(
                model = "y ~ temp*snp + temp*chr15_qtl + temp*chr21_qtl",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            dominance = list(
                model = "y ~ snp + dominance + temp*chr21_qtl + temp*chr15_qtl",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            linear = list(
                model = "y ~ snp + temp*chr21_qtl + temp*chr15_qtl",
                null_model = "y ~ temp*chr15_qtl + temp*chr21_qtl"
            ),
            no_temp = list(
                model = "y ~ snp + chr15_qtl + chr21_qtl",
                null_model = "y ~ chr15_qtl + chr21_qtl"
            )
        )

        # run the simulations and save the output
        i_max <- nrow(params_df)
        df <- lapply(1:i_max, simulate) |> rbindlist()
        df[, replicate := ${meta.rep}]
        df[, current_seed := ${seed}]
        df[, model := ${meta.model}]
        df[, simulated_locus := "${meta.locus_id}"]
        df[, simulated_locus := "${meta.locus_id}"]
        fwrite(df, "${meta.id}.csv.gz")
        """
}

process aggregate_results {
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(csv)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_aggregated.csv.gz")
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")

        read_file <- function(i, the_files = the_files, n_files = n_files) {
            message(sprintf("%s/%s", i, n_files))
            fread(the_files[[i]])
        }

        the_files <- list.files(pattern = ".*.csv.gz")
        stopifnot(length(the_files) == ${csv.size()})
        n_files <- length(the_files)
        df <- lapply(seq_along(the_files), read_file, the_files = the_files, n_files = n_files) |> rbindlist()
        fwrite(df, "${meta.id}_aggregated.csv.gz")
        """
}

workflow {
    read_pheno_covar ( params.pheno, params.covar )
    get_formulas ()
    get_qtls_and_models ( params.qtls, get_formulas.out )
    get_qtls_and_models.out.qtls
        .splitCsv ( header: true )
        .map {
            it.id = it.locus_id
            return ( it )
        }
        .set { qtls }
    get_qtls_and_models.out.models
        .splitCsv ( header: true )
        .map {
            it.id = it.locus_id + "_" + it.model
            return ( it )
        }
        .set { qtl_models }
    make_freq ( params.freq )
    make_pgen ( params.vcf )
    qtls.combine ( make_pgen.out ).combine ( make_freq.out ).set { make_grm_in_ch }
    make_grm ( make_grm_in_ch )
    make_grm.out
        .map { meta, id, bin -> [meta.locus_id, meta, id, bin] }
        .combine ( qtl_models.map { meta -> [meta.locus_id, meta] }, by: 0 )
        .map { match_tuple, meta1, id, bin, meta2 -> [meta2, id, bin] }
        .combine ( read_pheno_covar.out )
        .combine ( make_pgen.out )
        .combine ( get_formulas.out )
        .set { get_qtl_matrices_in_ch }
    get_qtl_matrices ( get_qtl_matrices_in_ch )
    fit_mixed_model ( get_qtl_matrices.out )
    fit_mixed_model.out
        .map { meta, mm -> [meta.locus_id, mm] }
        .combine ( get_qtl_matrices.out.map { meta, qtl_mat -> [meta.locus_id, meta, qtl_mat] }, by: 0 )
        .map { match_tuple, mm, meta, qtl_mat -> [meta, qtl_mat, mm] }
        .set { decorrelate_matrices_in_ch }
    decorrelate_matrices ( decorrelate_matrices_in_ch )
    fit_lm ( decorrelate_matrices.out )
    Channel.of ( ( 1..params.n_rep ) * params.random_seed ).flatten().set { seed }
    Channel.fromPath ( params.temperature_file ).set { temperature_file }
    fit_lm.out
        .combine ( temperature_file )
        .combine ( seed )
        .map {
            meta, lm_fit, temp_file, rep_random_seed ->
            def new_meta = meta.clone()
            def rep = rep_random_seed / params.random_seed
            new_meta.id = "${meta.id}_rep${rep}"
            new_meta.rep = rep
            [new_meta, lm_fit, temp_file, rep_random_seed]
        }
        .set { simulate_in_ch }
    simulate (
        simulate_in_ch,
        params.noise_levels,
        params.n_samples,
        params.r2,
        params.allele_freq
    )
    simulate.out.csv.map {
        meta, csv ->
        def new_meta = [
            id: meta.locus_id,
            locus_id: meta.locus_id,
            chr: meta.chr,
            lead_snp_id: meta.lead_snp_id,
        ]
        [new_meta, csv]
    }
    .groupTuple ( by: 0 )
    .set { aggregate_results_in_ch }
    aggregate_results ( aggregate_results_in_ch )
}
