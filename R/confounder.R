## reference
# https://github.com/biomedbigdata/namco/blob/647d3108a281eb0e36af31c44f5bf38d0c70c07d/app/R/utils.R#L480-L561
# https://github.com/biomedbigdata/namco/blob/647d3108a281eb0e36af31c44f5bf38d0c70c07d/app/R/server/confounding_server.R#L71-L100
# vagan: https://fromthebottomoftheheap.net/slides/advanced-vegan-webinar-2020/advanced-vegan

## adjust confounders: 
  # deseq2, edger, voom, metagenomeSeq:  ~ covariates + interested_var (condition)
  # ancom: adj_formula
  # ancombc: formula
  # maAsLin fixed effect https://forum.biobakery.org/t/confounding-factors/154

#' Confounder analysis
#' 
#' Confounding variables may mask the actual differential features. This
#' function utilizes constrained correspondence analysis (CCA) to measure the
#' confounding factors.
#' 
#' @param ps  a [`phyloseq::phyloseq-class`] object.
#' @param target_var character, the variable of interest 
#' @param norm  norm the methods used to normalize the microbial abundance data. See
#'   [`normalize()`] for more details.
#' @param confounders the confounding variables to be measured, if `NULL`, all 
#'   variables in the meta data will be analyzed.
#' @param permutations the number of permutations, see [`vegan::anova.cca()`].
#' @param ... extra arguments passed to [`vegan::anova.cca()`].
#' 
#' @return a `data.frame` contains three variables: confounder,  
#' pseudo-F and p value.
#' 
#' @examples 
#' data(caporaso)
#' confounder(caporaso, "SampleType", confounders = "ReportedAntibioticUsage")
#' 
#' @importFrom vegan cca
#' @importFrom phyloseq t taxa_are_rows
#' @importFrom stats anova
#' @export
confounder <- function(ps, 
                       target_var, 
                       norm = "none", 
                       confounders = NULL,
                       permutations = 999,
                       ...) {
    stopifnot(inherits(ps, "phyloseq"))
    abd <- otu_table(ps)
    abd <- normalize(abd, method = norm)
    if (taxa_are_rows(abd)) {
        abd <- as(t(abd), "matrix")
    }
    
    meta <- data.frame(sample_data(ps))
    confounders <- check_confounder(ps, target_var, confounders)
    confounders_meta <- meta[confounders]
    cca_out <- cca(abd ~ ., data = confounders_meta)
    
    cca_sig <- anova(cca_out, by = "terms", permutations = permutations, ...)
    cca_sig <- cca_sig[confounders, ]
    pseudo_F <- cca_sig$F
    pvalue <- cca_sig$`Pr(>F)`
    sig <- data.frame(
        confounder = row.names(cca_sig),
        pseudo_F = pseudo_F, 
        pvalue = pvalue
    )
    
    sig
}

