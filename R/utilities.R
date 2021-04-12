#' check whether tax abundance table is summarized or not
#' @noRd
check_tax_summarize <- function(ps) {
  taxa <- row.names(otu_table(ps))
  # taxa2 <- tax_table(ps)@.Data[, 1]

  # whether taxa is separated by `|`,
  # may be required to add extra separate strings in the future
  has_separate <- any(grepl("[|]", taxa))

  has_separate
}

#' check whether all names of taxonomic ranks include in available_ranks
#' @noRd
check_rank_names <- function(ps) {
  summarized <- check_tax_summarize(ps)
  if (summarized) {
    return(TRUE)
  }

  ps_ranks <- rank_names(ps)
  if (!all(ps_ranks %in% available_ranks)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' only first letter in lower case
#' @noRd
upper_firstletter <- function(x){
  paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
}

#' add prefix of taxonomic ranks for summarized data construct from original
#' lefse (galaxy lefse or python app) input, p__, k__
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param ranks character vector, prefix of ranks to add, e.g. "p", "c"
#' @importFrom phyloseq taxa_names<-
#' @noRd
add_prefix_summarized <- function(ps, ranks_prefix, sep = "|") {
  tax <- tax_table(ps)@.Data[, 1]
  tax_split <- strsplit(tax, split = sep, fixed = TRUE)

  if (max(lengths(tax_split)) != length(ranks_prefix)) {
    stop("The length of `ranks_prefix` muste be equal to number of taxonomic ranks.")
  }

  # ensure the ranks_prefix is contained in available_ranks
  # and in descending order
  available_prefix <- get_available_prefix(available_ranks)
  if (!all(ranks_prefix %in% available_prefix)) {
    stop("all elements of ranks_prefix must be contained in available_ranks")
  }
  ranks_prefix <- keep_prefix_desc(ranks_prefix, type = "ranks_prefix")

  tax_prefix <- purrr::map(
    tax_split,
    ~ paste0(ranks_prefix[1:length(.x)], "__", .x) %>%
      paste0(collapse = sep))
  tax_prefix <- do.call(rbind, tax_prefix)
  # row.names(tax_prefix) <- tax_prefix
  colnames(tax_prefix) <- paste0(ranks_prefix, collapse = sep)
  tax_table(ps) <- tax_table(tax_prefix)

  taxa_names(ps) <- tax_prefix[, 1]

  ps
}

# extract the first letter of taxonomic ranks as the prefixes of the ranks
get_available_prefix <- function(ranks) {
  substr(ranks, 1, 1) %>%
    tolower()
}

# keep prefix in descending order: "k" "p" "c" "o" "f" "g" "s"
keep_prefix_desc <- function(ranks_prefix, type = c("ranks", "ranks_prefix")) {
  type <- match.arg(type, choices = c("ranks", "ranks_prefix"))
  available_prefix <- get_available_prefix(available_ranks)
  idx_desc <- sort(match(ranks_prefix, available_prefix))

  if (type == "ranks") {
    return(available_ranks[idx_desc])
  } else {
    return(available_prefix[idx_desc])
  }
}

# check whether var in sample meta data, or raise an error
check_var_in_meta <- function(var, sample_meta) {
  stopifnot(inherits(sample_meta, "sample_data"))
  meta_nms <- names(sample_meta)
  if (!var %in% meta_nms) {
    stop(var, " must be one of variable of `sample_meta`", call. = FALSE)
  }
}

################################################################################
## preprocessing ps object
################################################################################

# preprocess of phyloseq object, including keep taxa in rows,
# filter taxa whose abundance is zero, fix duplicated tax
# filter samples whose library size is zero
#' @importFrom phyloseq prune_samples
preprocess_ps <- function(ps) {
  zero_sample <- check_samples(ps)
  if (!is.null(zero_sample)) {
    warning(
      "The library size of sample(s): ",
      paste(zero_sample, collapse = ", "),
      " is/are zero, and will be removed in the subsequent analysis."
    )

    keep <- setdiff(sample_names(ps), zero_sample)
    ps <- prune_samples(keep, ps)
  }

  # keep taxa in rows
  ps <- keep_taxa_in_rows(ps)
  # filter the taxa whose abundance is zero
  ps <- phyloseq_qc(ps)
  # fix duplicated tax
  ps <- fix_duplicate_tax(ps)

  ps
}

#' phyloseq quality control, remove otu/asv of which abundance is zero
#' @noRd
phyloseq_qc <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}

#' Transpose the phyloseq object to ensure taxa are in rows
#' @param ps a [phyloseq::phyloseq-class] object
#' @importMethodsFrom phyloseq t
#' @keywords internal
keep_taxa_in_rows <- function(ps) {
  if (!taxa_are_rows(ps)) {
    ps <- t(ps)
  }

  ps
}

#' Duplicated taxa: e.g. maybe multiple species (s__uncultured)
#' belong to different genera. append the upper level taxa to the taxa to
#' distinguish this duplicated taxa
#' @param ps [phyloseq::phyloseq-class] object or [phyloseq::taxonomyTable-class]
#'   object
#' @importFrom phyloseq tax_table<-
#' @keywords internal
#' @noRd
#' @references https://github.com/lch14forever/microbiomeViz/blob/94cbfe452a735aadf88733b27b8221a03f450a55/R/utils.R#L68-L86
fix_duplicate_tax <- function(ps) {
  # convert na to Unknown first
  ps <- fix_na_tax(ps)

  tax <- tax_table(ps)
  if (ncol(tax) == 1) {
    return(ps)
  }

  for(i in 2:ncol(tax)) {
    tax_uniq <-  unique(tax[, i])
    for(j in 1:length(tax_uniq)) {
      if(is.na(tax_uniq[j])) next
      ind <-  which(tax[, i] == as.character(tax_uniq[j]))
      if(length(unique(tax[ind, i - 1])) > 1) {
        tax[ind,i] <- paste(tax[ind, i - 1], tax[ind, i], sep = "_")
      }
    }
  }

  tax_table(ps) <- tax

  ps
}

#' set NA (missing) tax to its prefix, e.g. s__ (or s__unknown?)
#' @keywords internal
#' @noRd
fix_na_tax <- function(ps) {
  tax <- as.data.frame(tax_table(ps))

  tax_fixed <- purrr::imap_dfc(
    tax,
    ~ ifelse(is.na(.x), get_prefix(.y), .x)) %>%
    as.matrix()
  row.names(tax_fixed) <- taxa_names(ps)
  tax_table(ps) <- tax_fixed

  ps
}

# extract the prefix of taxonomic ranks
get_prefix <- function(ranks) {
  diff_ranks <- setdiff(ranks, available_ranks)
  if (length(diff_ranks) != 0) {
    stop(
      "ranks must be one of Kingdom, Phylum,",
      " Class, Order, Family, Genus, Species",
      call. = FALSE
    )
  }
  prefix <- substr(ranks, 1, 1) %>%
    tolower() %>%
    paste("__", sep = "")

  prefix
}


# `metagenomeSeq::cumNormStatFast()` requires counts of all samples at least
# have two non zero features. Thus, if there are samples with only one non-zer
# features, `cumNormStat()` is taken to compute the pth quantile.
# This function was used to select the function to calculate the quantile used
# for CSS norm factors estimation in metagenomeSeq.
select_quantile_func <- function(counts) {
  if (sum(colSums(counts > 0) > 1) < ncol(counts)) {
    fun_p <- metagenomeSeq::cumNormStat
  }
  else {
    fun_p <- metagenomeSeq::cumNormStatFast
  }

  fun_p
}


get_norm_method <- function(norm) {
  new_norm <- ifelse(
    is.numeric(norm),
    paste("per-sample normalized (sum of all taxa) to", norm),
    norm
  )

  new_norm
}


# check samples in ps, make sure at least one non zero features in a sample
check_samples <- function(ps) {
  if (!taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  lib_size <- colSums(otu_table(ps))
  zero_ind <- which(lib_size == 0)

  if (length(zero_ind) == 0) {
    return(NULL)
  }

  return(sample_names(ps)[zero_ind])
}


# remove samples with missing values for any variable specified in the group
remove_na_samples <- function(ps, group_var) {
  groups <- sample_data(ps)[[group_var]]
  na_idx <- is.na(groups)
  if (all(!na_idx)) {
    return(ps)
  }

  sample_nms <- sample_names(ps)
  warning(
    "remove sample(s): ", paste(sample_nms[na_idx], collapse = ","),
    ", which with missing value in `", group_var, "`"
  )
  ps <- phyloseq::prune_samples(sample_nms[!na_idx], ps)

  ps
}
