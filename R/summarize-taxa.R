#' Summarize taxa into a taxonomic level within each sample
#'
#' Provides summary information of the representation of a taxonomic levels
#' within each sample.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param level taxonomic level to summarize, default the top level rank of the
#'  `ps`.
#' @param absolute logical, whether return the absolute abundance or
#'   relative abundance, default `FALSE`.
#' @param sep a character string to separate the taxonomic levels.
#' @return a [`phyloseq::phyloseq-class`] object, where each row represents a
#'   taxa, and each col represents the taxa abundance of each sample.
#' @export
#' @examples
#' data(enterotypes_arumugam)
#' summarize_taxa(enterotypes_arumugam)
summarize_taxa <- function(ps,
    level = rank_names(ps)[1],
    absolute = TRUE,
    sep = "|") {

    # return ps if it has been summarized
    summarized <- check_tax_summarize(ps)
    if (summarized) {
        otu_summarized <- otu_table(ps) %>%
            add_missing_levels()
        tax_summarized <- row.names(otu_summarized) %>%
            matrix() %>%
            tax_table()
        row.names(tax_summarized) <- row.names(otu_summarized)
        return(phyloseq(otu_summarized, tax_summarized, sample_data(ps)))
    }

    if (!has_prefix(ps)) {
        ps <- add_prefix(ps)
    }

    ps_ranks <- rank_names(ps)
    if (!level %in% ps_ranks) {
        stop("`level` must in the ranks of `ps` (rank_names(ps))")
    }

    ind <- match(level, ps_ranks)
    levels <- ps_ranks[ind:length(ps_ranks)]
    res <- purrr::map(
        levels,
        ~ .summarize_taxa_level(
            ps,
            rank = .x,
            absolute = absolute,
            sep = sep
        )
    )
    tax_nms <- purrr::map(res, row.names) %>% unlist()
    res <- bind_rows(res)
    row.names(res) <- tax_nms

    otu_summarized <- otu_table(res, taxa_are_rows = TRUE)
    tax_summarized <- row.names(otu_summarized) %>%
        matrix() %>%
        tax_table()
    row.names(tax_summarized) <- row.names(otu_summarized)
    row.names(otu_summarized) <- row.names(tax_summarized)

    # To ensure the rank of the summarized object is valid (one of  "domain"  
    # "kingdom" "phylum"  "class"   "order"   "family"  "genus"   "species"),  
    # set it (column names of tax_summarized) as the top level rank in the ps
    # object.
    #  
    # rank_prefix <- extract_prefix(ps_ranks)
    # colnames(tax_summarized) <- paste0(rank_prefix, collapse = sep)
    colnames(tax_summarized) <- ps_ranks[1]

    return(phyloseq(otu_summarized, tax_summarized, sample_data(ps)))
}

#' extract prefix of names of the taxonomic ranks
#' @noRd
extract_prefix <- function(ranks) {
    if (inherits(ranks, "phyloseq")) {
        ranks <- rank_names(ranks)
    }

    tolower(substr(ranks, 1, 1))
}

#' Summarize the taxa for the specific rank
#' @noRd
.summarize_taxa_level <- function(ps,
    rank_name,
    absolute = TRUE,
    sep = "|") {
    if (!absolute) {
        ps <- transform_sample_counts(ps, function(x) x / sum(x))
    }

    otus <- otu_table(ps)
    otus_extend <- slot(otus, ".Data") %>%
        tibble::as_tibble()

    taxas <- tax_table(ps)@.Data %>%
        tibble::as_tibble()

    ranks <- setdiff(available_ranks, "Summarize")
    rank_level <- match(rank_name, ranks)
    select_ranks <- intersect(ranks[seq_len(rank_level)], rank_names(ps))

    consensus <- taxas[, select_ranks] %>%
        purrr::pmap_chr(paste, sep = sep)
    otus_extend$consensus <- consensus

    taxa_summarized <- group_split(otus_extend, consensus) %>%
        purrr::map(.sum_consensus)
    taxa_summarized <- do.call(rbind, taxa_summarized)
    # filter taxa of which abundance is zero
    ind <- rowSums(taxa_summarized) != 0
    taxa_summarized <- taxa_summarized[ind, ]

    taxa_summarized
}

#' sum all otus which belongs to the same taxa
#' @noRd
.sum_consensus <- function(x) {
    consensus <- unique(x$consensus)
    if (length(consensus) != 1) {
        stop("consensus in the same group muste be the same")
    }

    x$consensus <- NULL
    res <- as.data.frame(t(colSums(x)))
    row.names(res) <- consensus

    return(res)
}

# check whether taxa has prefix or not
has_prefix <- function(ps) {
    sample_tax <- tax_table(ps)[1, 1]
    if (substr(sample_tax, 2, 3) == "__") {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# add ranks prefix, e.g k__, p__, only worked for unsummarized data
add_prefix <- function(ps) {
    tax <- as(tax_table(ps), "matrix") %>%
        as.data.frame()
    lvl <- colnames(tax)
    prefix <- get_prefix(lvl)

    tax_new <- mapply(function(x, y) paste0(x, y), 
        prefix, tax, SIMPLIFY = FALSE)
    tax_new <- do.call(cbind, tax_new)
    row.names(tax_new) <- row.names(tax)
    colnames(tax_new) <- lvl
    tax_table(ps) <- tax_new

    ps
}
