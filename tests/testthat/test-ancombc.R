test_that("ancombc works correctly", {
    if (FALSE) {
        # result of ancombc package, example from the package vignette
        # atalas1006 from microbiome package
        pseq <- subset_samples(atlas1006, time == 0)
        sample_data(pseq)$bmi_group <- recode(
            sample_data(pseq)$bmi_group,
            underweight = "lean",
            lean = "lean",
            overweight = "overweight",
            obese = "obese",
            severeobese = "obese",
            morbidobese = "obese"
        )
        sample_data(pseq)$nation <- dplyr::recode(
            sample_data(pseq)$nationality,
            Scandinavia = "NE",
            UKIE = "NE",
            SouthEurope = "SE",
            CentralEurope = "CE",
            EasternEurope = "EE"
        )

        phylum_data <- aggregate_taxa(pseq, "Phylum")
        tax_table(phylum_data) <- tax_table(phylum_data)[, 1]

        out <- ANCOMBC::ancombc(
            phyloseq = phylum_data,
            formula = "nation",
            p_adj_method = "holm",
            prv_cut = 0.10,
            lib_cut = 1000,
            group = "nation",
            struc_zero = FALSE,
            neg_lb = FALSE,
            tol = 1e-5,
            max_iter = 100,
            conserve = FALSE,
            alpha = 0.05,
            global = TRUE
        )
        group_lvls <- levels(phyloseq::sample_data(phylum_data)[["nation"]])
        ef <- out$res$lfc

        # extract enrich groups according the effect size
        #
        # https://github.com/FrederickHuangLin/ANCOMBC/issues/8
        # The first level (CE) of a categorical value will be set as a
        # reference level by default in R.
        # Yes, the W statistics from the primary result (res), which aims to
        # test for the effects of covariates of interest, can be either positive
        # or negative and its sign indicates the direction;
        get_enrich_group <- function(ef, group_lvls) {
            if (all(ef < 0)) {
                enrich_group <- group_lvls[1]
            } else {
                enrich_group <- group_lvls[which.max(ef) + 1]
            }

            enrich_group
        }

        enrich_group <- apply(ef, 1, get_enrich_group, group_lvls)
        global_res <- out$res_global
        global_res$enrich_group <- enrich_group
        global_res <- global_res[global_res$diff_abn, ]
        res <- data.frame(
            feature = paste0("p__", rownames(global_res)),
            enrich_group = global_res$enrich_group,
            ef_W = global_res$W,
            pvalue = global_res$p_val,
            padj = global_res$q_val
        )
        rownames(res) <- paste0("marker", 1:nrow(res))


        out2 <- run_ancombc(
            phylum_data,
            p_adjust = "holm",
            prv_cut = 0.10,
            lib_cut = 1000,
            group = "nation",
            struc_zero = FALSE,
            neg_lb = FALSE,
            tol = 1e-5,
            max_iter = 100,
            conserve = FALSE,
            pvalue_cutoff = 0.05
        )
        marker <- data.frame(marker_table(out2))

        # TRUE
        identical(marker, res)
    }
})
