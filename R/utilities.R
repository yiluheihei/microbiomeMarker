#' phyloseq quality control, remove otu/asv of which abundance is zero
#' @noRd
phyloseq_qc <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}

#' get the means abundance of each class for each feature
#' @noRd
get_feature_enrich_group <- function(class, feature) {
  feature$class <- class
  feature_mean <- group_by(feature, class) %>%
    summarise(colMeans(cur_data(), ) %>% bind_rows())
  feature_enrich_index <- select(feature_mean, -class) %>%
    purrr::map_dbl(which.max)

  feature_enrich_group <- feature_mean$class[feature_enrich_index]
  names(feature_enrich_group) <- names(feature)[names(feature) != "class"]
  feature_max_mean <- purrr::map2_dbl(
    select(feature_mean, -class),
    feature_enrich_index,
    ~ .x[.y]
  )
  # log10(max(feature_max_mean), 1)
  feature_max_mean[feature_max_mean < 1] <- 1

  return(list(
    group = feature_enrich_group,
    log_max_mean = log10(feature_max_mean)
  ))
}

#' check bootstrap sample, not contast within class and too fewer in each
#' class
#'
#' @param feature_abudance data.frame, significant feature abundance,
#'   where columns represents samples and rows represents features
#' @param sample_indx numeric vector, sample index be bootstraped
#' @param sample_min integer, min samples in each class
#' @param class character vector, class of all samples
#'
#' @noRd
check_bootstrap_sample <- function(feature_abudance,
                                   sample_indx,
                                   sample_min,
                                   class) {
  feature_abudance <- feature_abudance[, sample_indx]
  class_n <- length(unique(class))
  class <- class[sample_indx]

  if (length(unique(class)) < class_n) {
    return(FALSE)
  }

  for (cls in unique(class)) {
    if (sum(class == cls) < sample_min) {
      return(FALSE)
    }

    # sig feature smaller than min sample count
    cls_abundance <- feature_abudance[, class == cls]

    for (i in seq_along(nrow(cls_abundance))) {
      unique_abd <- length(unique(cls_abundance[i, ]))

      if ((unique_abd <= sample_min && sample_min > 1) || (unique_abd <= 1 && sample_min == 1)) {
        return(FALSE)
      }
    }
  }

  return(TRUE)
}

#' bootstrap iteration of samples for lda analysis
#' @importFrom dplyr group_by summarize cur_data
#' @noRd
bootstap_lda <- function(feature_abundance,
                            boot_n,
                            class,
                            sample_fract) {
  ldas <- purrr::rerun(
    boot_n,
    bootstap_lda_one(
      feature_abundance,
      class,
      sample_fract
    )
  ) %>%
    purrr::transpose() %>%
    purrr::map(bind_rows) %>%
    bind_rows()

  mean_lds <- colMeans(ldas)
  mean_lds <- sign(mean_lds) * log10(1 + abs(mean_lds))

  mean_lds
}

bootstap_lda_one <- function(feature_abundance,
                            class,
                            sample_fract) {
  sample_groups <- unique(class)
  class_count <- table(class)
  feature_abundance$class <- class

  feature_abundance <- preprocess_feature_all(feature_abundance, class)

  sample_n <- nrow(feature_abundance)
  random_n <- floor(sample_n * sample_fract)
  class_n <- length(sample_groups)
  sample_min <- floor(min(class_count) * sample_fract * sample_fract * 0.5) %>%
    max(1)

  # class vs class
  pairs <- utils::combn(sample_groups, 2, simplify = FALSE) %>%
    purrr::map(sort, decreasing = TRUE)

  while (TRUE) {
    # random select samples using bootstrap method
    sample_indx <- sample(sample_n, random_n, replace = TRUE)

    is_checked <- check_bootstrap_sample(
      feature_abundance,
      sample_indx,
      sample_min,
      class
    )
    if (is_checked) {
      break
    }
  }

  lda <- purrr::map(
    pairs,
    ~ cal_pair_lda(feature_abundance, sample_indx, .x)
  )
  names(lda) <- purrr::map(pairs, paste, collapse = " -VS- ")

  lda
}

#' calculate lda score of single pair groups
#' @noRd
cal_pair_lda <- function(feature_abundance,
                         sample_indx,
                         pair) {
  sample_feature_abundance <- feature_abundance[sample_indx, ]

  # reference lefse.py in lefse
  lda_res <- suppressWarnings(
    MASS::lda(
      class ~ .,
      data = sample_feature_abundance,
      tol = 1.0e-10
    )
  )
  w <- lda_res$scaling[, 1]
  w_unit <- w/sqrt(sum(w^2))
  feature_remove_class <- sample_feature_abundance[-1]

  # not support subclass and subject argument in lefse

  ld <- as.matrix(feature_remove_class)%*%w_unit
  group1_indx <- sample_feature_abundance$class == pair[1]
  effect_size <- abs(mean(ld[group1_indx]) - mean(ld[-group1_indx]))
  wfinal <- w_unit * effect_size
  lda_means <- lda_res$means
  lda_row_nms <- row.names(lda_means)
  feature_n <- ncol(lda_means)
  coeff <- ifelse(is.nan(wfinal), 0, abs(wfinal))

  res <- purrr::map(
    pair,
    function(x) {
      if (x %in% lda_row_nms) {
        lda_means[x, ]
      } else {
        rep(0, feature_n)
      }
    }
  )
  names(res) <- pair

  feature <- names(feature_remove_class)
  lda_score <- purrr::map_dbl(
    seq_along(feature),
    function(i) {
      gm <- abs(res[[1]][i] - res[[2]][i])
      return(gm + coeff[i] * 0.5)
    }
  )
  names(lda_score) <- feature

  lda_score

}


#' feature abundance preprocess
#' @noRd
preprocess_feature_all <- function(x, class) {
  res <- group_by(
    x,
    class
  ) %>%
    summarise(purrr::map_df(cur_data(), preprocess_feature))

  res
}

preprocess_feature <- function(x) {
  if (length(unique(x)) <= max(length(x)*0.5, 4)) {
    x <- purrr::map_dbl(x, ~ abs(.x + rnorm(1, 0, max(.x*0.05, 0.01))))
  }

  x
}
