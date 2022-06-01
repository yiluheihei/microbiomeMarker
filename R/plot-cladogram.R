##  codes for cladogram plot are modified from microbiomeViz
##  https://github.com/lch14forever/microbiomeViz

#' @title plot cladogram of micobiomeMaker results
#'
#' @param mm a [microbiomeMarker-class] object
#' @param color a color vector, used to highlight the clades of microbiome
#'   biomarker. The values will be matched in order (usually alphabetical) with
#'   the groups. If this is a named vector, then the colors will be matched
#'   based on the names instead.
#' @param only_marker logical, whether show all the features or only 
#'   markers in the cladogram, default `FALSE`.
#' @param branch_size numeric, size of branch, default `0.2`
#' @param alpha alpha parameter for shading, default `0.2`
#' @param clade_label_level max level of taxa used to label the clade, other
#' level of taxa  will be shown on the side.
#' @param clade_label_font_size font size of the clade label, default 4.
#' @param node_size_scale the parameter 'a' controlling node size:
#' `node_size=a*log(relative_abundance) + b`
##' @param node_size_offset the parameter 'b' controlling node size:
##' `node_size=a*log(relative_abundance) + b`
##' @param annotation_shape shape used for annotation, default `22`
##' @param annotation_shape_size size used for annotation shape, default `5`
##' @param  group_legend_param,marker_legend_param a list specifying
##'   extra parameters of group legend and marker legend, such as `direction` (
##'   the direction of the guide), `nrow` (the desired number of rows of
##'   legends). See [`ggplot2::guide_legend()`] for more details.
#' @return a ggtree object
#' @importFrom tidytree treedata
#' @importFrom ggplot2 geom_point theme element_blank geom_rect guides
#' guide_legend aes_ scale_shape_manual
#' @importFrom ggtree ggtree geom_hilight geom_point2 geom_cladelabel
#' @author Chenhao Li, Guangchuang Yu, Chenghao Zhu, Yang Cao
#' @seealso [ggtree::ggtree()]
#' @export
#' @references This function is modified from `clada.anno` from microbiomeViz.
#' @examples
#' data(kostic_crc)
#' kostic_crc_small <- phyloseq::subset_taxa(
#'     kostic_crc,
#'     Phylum %in% c("Firmicutes")
#' )
#' mm_lefse <- run_lefse(
#'     kostic_crc_small,
#'     wilcoxon_cutoff = 0.01,
#'     group = "DIAGNOSIS",
#'     kw_cutoff = 0.01,
#'     multigrp_strat = TRUE,
#'     lda_cutoff = 4
#' )
#' plot_cladogram(mm_lefse, color = c("darkgreen", "red"))
plot_cladogram <- function(mm,
    color,
    only_marker = FALSE,
    branch_size = 0.2,
    alpha = 0.2,
    node_size_scale = 1,
    node_size_offset = 1,
    clade_label_level = 4,
    clade_label_font_size = 4,
    annotation_shape = 22,
    annotation_shape_size = 5,
    group_legend_param = list(),
    marker_legend_param = list()) {
    ps <- create_ps_from_mm(mm, only_marker = only_marker)
    tree <- get_treedata_phyloseq(ps) %>%
        generate_taxa_tree(size = branch_size)

    annotation <- generate_cladogram_annotation(
        mm@marker_table,
        color = color
    )

    # background highlight
    annotation_info <- dplyr::left_join(
        annotation,
        tree$data,
        by = c("node" = "label")
    ) %>%
        mutate(
            label = .data$node,
            id = .data$node.y,
            level = as.numeric(.data$node_class)
        )

    hilight_para <- dplyr::transmute(
        annotation_info,
        node = .data$id,
        fill = .data$color,
        alpha = alpha,
        extend = get_offset(.data$level)
    )
    hilights_g <- purrr::pmap(hilight_para, geom_hilight)
    tree <- purrr::reduce(hilights_g, `+`, .init = tree)

    # hilight legend
    # default: colors were matched for in alphabetical of groups, which requires
    # arrange hilights_df according to enrich_group
    hilights_df <- dplyr::distinct(
        annotation_info,
        .data$enrich_group,
        .data$color
    ) %>%
        arrange(.data$enrich_group)
    hilights_df$x <- 0
    hilights_df$y <- 1
    group_legend_param <- c(
        group_legend_param,
        list(
            title = NULL,
            order = 1,
            override.aes = list(fill = hilights_df$color)
        )
    )
    group_lgd <- do.call(guide_legend, group_legend_param)
    tree <- tree +
        geom_rect(
            aes_(xmin = ~x, xmax = ~x, ymax = ~y, ymin = ~y, 
                fill = ~enrich_group),
            data = hilights_df, inherit.aes = FALSE
        ) +
        guides(fill = group_lgd)

    # set nodes color and size
    nodes_colors <- rep("white", nrow(tree$data))
    nodes_colors[annotation_info$id] <- annotation_info$color
    node_size <- node_size_scale * log(tree$data$abd) + node_size_offset
    tree$data$node_size <- node_size
    tree <- tree +
        geom_point2(aes(size = I(node_size)), fill = nodes_colors, shape = 21)

    ## add clade labels
    clade_label <- dplyr::transmute(
        annotation_info,
        node = .data$id,
        offset = get_offset(.data$level) - 0.4,
        angle = purrr::map_dbl(.data$id, get_angle, tree = tree) + 90,
        label = .data$label,
        fontsize = clade_label_font_size,
        barsize = 0,
        hjust = 0.5,
        level = .data$level
    ) %>%
        dplyr::arrange(desc(.data$level))
    ind <- clade_label$level < clade_label_level
    short_label <- get_short_label_id(clade_label, clade_label_level)
    clade_label_para <- mutate(
        clade_label,
        label = c(.data$label[!ind], short_label),
        level = NULL
    )
    clade_label_g <- purrr::pmap(clade_label_para, geom_cladelabel)
    tree <- purrr::reduce(clade_label_g, `+`, .init = tree)

    ## add guide labels
    guide_label <- clade_label[ind, ] %>%
        mutate(
            label2 = paste0(short_label, ": ", .data$label),
            color = annotation_info$color[
                match(.data$label, annotation_info$label)]
        )

    # marker annotation, legend
    marker_legend_param <- c(
        marker_legend_param,
        list(
            p = tree,
            color = guide_label$color,
            label = guide_label$label2,
            shape = annotation_shape,
            size = annotation_shape_size
        )
    )
    p <- do.call(set_marker_annotation, marker_legend_param) +
        theme(legend.position = "right", legend.title = element_blank())

    p
}

#' Get short label id
#' @keywords internal
#' @noRd
get_short_label_id <- function(clade_label, clade_label_level) {
    ind <- clade_label$level < clade_label_level
    unique_id <- get_unique_id(sum(ind))
    short_label <- unique_id[seq_len(sum(ind))]

    short_label
}

#' Get unique id for short label annotation
#' {so}/questions/21681785/repeating-vector-of-letters/21689613#21689613
#' @keywords internal
#' @noRd
get_unique_id <- function(n, depth = 1) {
    args <- lapply(seq_len(depth), FUN = function(x) letters)
    x <- do.call(expand.grid, args = list(args, stringsAsFactors = FALSE))
    x <- x[, rev(names(x)), drop = FALSE]
    x <- do.call(paste0, x)
    if (n <= length(x)) {
        return(x[seq_len(n)])
    }

    return(c(x, get_unique_id(n - length(x), depth = depth + 1)))
}

#' Generate tree data from phyloseq object
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param sep character, separate between different levels of taxa, default `|`
#' @author Yang Cao
#' @return a [`tidytree::treedata-class`] object
#' @keywords internal
get_treedata_phyloseq <- function(ps, sep = "|") {
    if (!taxa_are_rows(ps)) {
        stop("Requires taxa in rows of phyloseq")
    }

    taxa <- tax_table(ps)
    otu <- otu_table(ps)
    row.names(otu) <- taxa@.Data[, 1]
    taxa_nms <- row.names(otu)

    tree_table <- data.frame(
        taxa = taxa_nms,
        abd = rowMeans(otu),
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            taxa = paste("r__Root", .data$taxa, sep = "|"),
            abd = .data$abd / max(.data$abd) * 100
        )

    taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
    nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
    # add root node
    nodes <- c("r__Root", nodes)

    ## data may not contain all the seven ranks of the taxa, such as
    ## enterotypes_arumugam only contains Phylum and Genus ranks
    taxa_deepest <- taxa_split[[which.max(lengths(taxa_split))]]
    prefix <- vector("character", length(taxa_deepest))
    for (i in seq_along(taxa_deepest)) {
        if (!grepl("__$", taxa_deepest[i])) {
            prefix[i] <- gsub("(.*)__.*", "\\1", taxa_deepest[i])
        } else {
            pos <- nchar(taxa_deepest[i]) - 2
            prefix[i] <- substr(taxa_deepest[i], pos, pos)
        }
    }

    levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
        factor(levels = rev(prefix))

    nodes_parent <- purrr::map_chr(
        taxa_split,
        ~ .x[length(.x) - 1]
    )
    # root must be a parent node
    nodes_parent <- c("root", nodes_parent)

    ## tips comes first ?
    is_tip <- !nodes %in% nodes_parent
    index <- vector("integer", length(is_tip))
    index[is_tip] <- seq_len(sum(is_tip))
    index[!is_tip] <- (sum(is_tip) + 1):length(is_tip)

    edges <- cbind(
        parent = index[match(nodes_parent, nodes)],
        child = index
    )
    edges <- edges[!is.na(edges[, 1]), ]

    # not label the tips
    node_label <- nodes[!is_tip]

    phylo <- structure(
        list(
            edge = edges,
            node.label = node_label,
            tip.label = nodes[is_tip],
            edge.length = rep(1, nrow(edges)),
            Nnode = length(node_label)
        ),
        class = "phylo"
    )

    mapping <- data.frame(
        node = index,
        abd = c(100, tree_table$abd),
        node_label = nodes,
        stringsAsFactors = FALSE
    )
    mapping$node_class <- levels

    tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
}

#' generate taxa hierarchy tree
#' @noRd
generate_taxa_tree <- function(treedata,
    size = 0.2,
    layout = "circular") {
    ggtree::ggtree(treedata, size = size, layout = layout)
}

#' generate annotaion data for cladogram plot
#' @param marker data.frame
#' @param color a color vector, used to highlight the clades of ggtree
#' @param sep seprator between different of levels of taxa
#' @noRd
generate_cladogram_annotation <- function(marker,
    color,
    sep = "|") {
    enrich_group <- marker$enrich_group
    if (length(color) != length(unique(enrich_group))) {
        stop("the number of colors must be equal to ",
            "the number of enriched groups.")
    }

    feature <- marker$feature
    label <- strsplit(feature, split = sep, fixed = TRUE) %>%
        purrr::map_chr(utils::tail, n = 1)
    label_level <- lengths(strsplit(feature, sep, fixed = TRUE))

    # may be no marker are identified enriched in some groups
    # drop the levels of this groups if the enrich_group is a factor
    if (inherits(enrich_group, "factor")) {
        enrich_group <- droplevels(enrich_group)
    }

    # named colors: set the colors based on the matched names to groups
    if (is.vector(color) && !is.null(names(color))) {
        if (!all(names(color) %in% enrich_group)) {
            stop("names of `color` muste be contained in enriched groups")
        }
        color <- color[match(enrich_group, names(color))]
    } else {
        # colors will be matched in order (usually alphabetical) with the groups
        names(color) <- sort(unique(enrich_group))
        color <- color[match(enrich_group, names(color))]
    }

    annotation <- data.frame(
        node = label,
        color = color,
        enrich_group = enrich_group,
        stringsAsFactors = FALSE
    )

    annotation
}

#' get clade background offset
#' @noRd
get_offset <- function(x) {
    (x * 0.2 + 0.2)^2
}

#' get the mean angle of a clade
#' @noRd
get_angle <- function(tree, node) {
    if (length(node) != 1) {
        stop("The length of `node` must be 1")
    }
    tree_data <- tree$data
    sp <- tidytree::offspring(tree_data, node)$node
    sp2 <- c(sp, node)
    sp.df <- tree_data[match(sp2, tree_data$node), ]
    mean(range(sp.df$angle))
}

#' set legend for multiple geom_cladelabel layers
#'
#' This function can be used to set the microbiome marker annotations
#'
#' @param p a ggtree object
#' @param color a color vector
#' @param label a character vector, with the same length with `color`
#' @param shape shape of label, default `22`
#' @param size size of shape, default `5`
#' @param ... extra arguments passed to  [ggplot2::guide_legend()],
#' e.g. `ncol`, more details see [ggplot2::guide_legend()].
#' @seealso [ggplot2::guide_legend()]
#' @return an updated `ggtree` object
#' @importFrom ggplot2 geom_point aes_ scale_shape_manual guides guide_legend
#' @noRd
set_marker_annotation <- function(p,
    color,
    label,
    size = 5,
    shape = 22,
    ...) {
    dat <- data.frame(
        color = color,
        label = label,
        stringsAsFactors = FALSE
    )

    # suppress warning: The shape palette can deal with a maximum of 6 discrete
    # values because more than 6 becomes difficult to discriminate; you have 18.
    # Consider specifying shapes manually if you must have them.
    # using scale_shape_manual
    p <- p +
        geom_point(
            data = dat, inherit.aes = FALSE,
            aes_(x = 0, y = 0, shape = ~label),
            size = 0, stroke = 0,
        ) +
        scale_shape_manual(values = rep(shape, nrow(dat)), limits = dat$label) +
        guides(
            shape = guide_legend(
                override.aes = list(
                    size = size,
                    shape = shape,
                    fill = dat$color
                ),
                order = 2,
                ...
            )
        )

    p
}
