##  codes for lefse cladogram plot are modified from microbiomeViz
##  https://github.com/lch14forever/microbiomeViz

#' @title plot cladogram of lefse results
#'
#' @param mm a [microbiomeMarker-class] object
#' @param color a color vector, used to highlight the clades of micribome
#' biomaker
#' @param branch_size numberic, size of branch, default `0.2`
#' @param alpha alpha parameter for shading, default `0.2`
#' @param clade_label_level max level of taxa used to label the clade, other
#' level of taxa  will be shown on the side
#' @param node_size_scale the parameter 'a' controlling node size:
#' `node_size=a*log(relative_abundance) + b`
##' @param node_size_offset the parameter 'b' controlling node size:
##' `node_size=a*log(relative_abundance) + b`
#' @return a ggtree object
#' @importFrom tidytree treedata
#' @importFrom ggplot2 geom_point theme element_blank
#' @importFrom ggtree ggtree geom_hilight geom_point2 geom_cladelabel
#' @author Chenhao Li, Guangchuang Yu, Chenghao Zhu, Yang Cao
#' @seealso [ggtree::ggtree()]
#' @export
#' @references This function is scratch from `clada.anno` from microbiomeViz.
#' \url{https://github.com/lch14forever/microbiomeViz/blob/master/R/visualizer.R}
#' @description annotate a ggtree plot to highlight certain clades
lefse_cladogram <- function(mm,
                            color =
                              scales::hue_pal()(length(unique(mm@microbiome_marker$enrich_group))),
                            branch_size = 0.2,
                            alpha = 0.2,
                            node_size_scale = 1,
                            node_size_offset = 1,
                            clade_label_level = 4){
  ps <- phyloseq(mm@otu_table, mm@tax_table)
  tree <- get_treedata_phyloseq(ps) %>%
    generate_taxa_tree(size = branch_size)


  annotation <- generate_cladogram_annotation(
    mm@microbiome_marker,
    color = color
  )

  # backgroup hilight
  annotation_info <- dplyr::left_join(
    annotation,
    tree$data,
    by = c("node" = "label")) %>%
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

  # set nodes color and size
  nodes_colors <- rep("white", nrow(tree$data))
  nodes_colors[annotation_info$id] <- annotation_info$color
  node_size <- node_size_scale*log(tree$data$abd) + node_size_offset
  tree$data$node_size <- node_size
  tree <- tree +
    geom_point2(aes(size = I(node_size)), fill = nodes_colors, shape = 21)


  ## add clade labels
  clade_label <-dplyr::transmute(
    annotation_info,
    node = .data$id,
    offset = get_offset(.data$level) - 0.4,
    angle = purrr::map_dbl(.data$id, get_angle, tree = tree) + 90,
    label = .data$label,
    fontsize = 1.5 + sqrt(.data$level),
    barsize = 0,
    hjust = 0.5,
    level = .data$level
  ) %>%
    dplyr::arrange(desc(.data$level))
  ind <- clade_label$level < clade_label_level
  short_label <- letters[1:sum(ind)]
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
      anno_shape = purrr::map_int(short_label, utf8ToInt),
      label2 = paste0(short_label, ": ", .data$label))

  tree +
    geom_point(
      data = guide_label,
      aes(x = 0, y = 0, shape = factor(.data$label2)),
      size = 0,
      stroke = 0) +
    theme(legend.position = c(1,0.5),
          legend.title = element_blank())
}

#' Generate tree data from phyloseq object
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param sep character, separate between different levels of taxa, default `|`
#' @author Yang Cao
#' @return a [`tidytree::treedata-class`] object
get_treedata_phyloseq <- function(ps, sep = "|") {
  if (!taxa_are_rows(ps)) {
    stop("Requires taxa in rows of phyloseq")
  }

  taxa <- tax_table(ps)
  feature <- otu_table(ps)

  is_summarized <- check_tax_summarize(ps)
  if (is_summarized) {
    row.names(feature) <- taxa@.Data[, 1]
    feature <- add_missing_levels(feature)
  } else {
    feature <- summarize_taxa(ps, sep = sep)
  }

  taxa_nms <- row.names(feature)
  has_prefix <- check_tax_prefix(taxa_nms)
  if (!has_prefix) {
    taxa_nms <- add_tax_level(taxa_nms, sep = sep)
  }

  tree_table <- data.frame(
    taxa = taxa_nms,
    abd = rowMeans(feature),
    stringsAsFactors = FALSE) %>%
    mutate(
      taxa =  paste("r__Root", .data$taxa, sep = "|"),
      abd = .data$abd/max(.data$abd)*100
    )

  taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
  nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
  # add root node
  nodes <- c("r__Root", nodes)

  # levels used for extend of clade label
  levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
    factor(
    levels = rev(c("r" , "k", "p", "c", "o", "f", "g", "s"))
  )

  nodes_parent <- purrr::map_chr(
    taxa_split,
    ~ .x[length(.x) - 1]
  )
  # root must be a parent node
  nodes_parent <- c("root", nodes_parent)

  ## tips comes first ?
  is_tip <- !nodes %in% nodes_parent
  index <- vector("integer", length(is_tip))
  index[is_tip] <- 1:sum(is_tip)
  index[!is_tip] <- (sum(is_tip)+1):length(is_tip)

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
                               layout = 'circular'){
  ggtree::ggtree(treedata, size = size, layout = layout)
}

#' generate annotaion data for cladogram plot
#' @param lefse_out a data.frame
#' @param color a color vector, used to highlight the clades of ggtree
#' @param sep seprator between different of levels of taxa
#' @noRd
generate_cladogram_annotation <- function(lefse_out,
                                          color,
                                          sep = "|") {
  feature <- lefse_out$feature
  has_prefix <- check_tax_prefix(feature)
  if (!has_prefix) {
    feature <- add_tax_level(feature)
  }
  label <- strsplit(feature, split = sep, fixed = TRUE) %>%
    purrr::map_chr(utils::tail, n =1)
  label_level <- lengths(strsplit(feature, sep, fixed = TRUE))
  color <- rep(color, times = table(lefse_out$enrich_group))

  annotation <- data.frame(
    node = label,
    color = color,
    stringsAsFactors = FALSE
  )

  annotation
}

#' get clade background offset
#' @noRd
get_offset <- function(x) {(x*0.2+0.2)^2}

#' get the mean angle of a clade
#' @noRd
get_angle <- function(tree, node){
  if (length(node) != 1) {
    stop("The length of `node` must be 1")
  }
  tree_data <- tree$data
  sp <- tidytree::offspring(tree_data, node)$node
  sp2 <- c(sp, node)
  sp.df <- tree_data[match(sp2, tree_data$node),]
  mean(range(sp.df$angle))
}

