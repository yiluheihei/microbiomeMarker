#' 16S rRNA data from "Moving pictures of the human microbiome"
#'
#' 16S read counts and phylogenetic tree file of 34 Illumina samples derived
#' from Moving Pictures of the Human Microbiome (Caporaso et al.) Group label:
#' gut, left palm, right palm, and tongue - indicating different sampled body
#' sites.
#'
#' @format a [phyloseq::phyloseq] object
#' @references
#' Caporaso, et al. Moving pictures of the human microbiome. Genome Biol 12, R50 (2011).
#'
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r50}
#'
#' @source Data was downloaded from https://www.microbiomeanalyst.ca:
#' \url{https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/treebiom.zip}
#'
#'
#' @name data-caporaso_phyloseq
#' @aliases caporaso_phyloseq
#' @docType data
#' @author Yang Cao
"caporaso_phyloseq"

#' 16S rRNA data of 94 patients from CID 2012
#'
#' Data from a cohort of 94 Bone Marrow Transplant patients previously published
#' on in CID
#'
#' @format a [phyloseq::phyloseq] object
#' @references
#' Ying, et al. Intestinal Domination and the Risk of Bacteremia in Patients
#' Undergoing Allogeneic Hematopoietic Stem Cell Transplantation,
#' Clinical Infectious Diseases, Volume 55, Issue 7, 1 October 2012,
#' Pages 905–914,
#'
#' \url{https://academic.oup.com/cid/article/55/7/905/428203}
#'
#' @source \url{https://github.com/ying14/yingtools2/tree/master/data}
#' @name data-cid.phy
#' @aliases cid.phy
#' @docType data
#' @author Yang Cao
"cid.phy"

#' IBD stool samples
#'
#' 43 pediatric IBD stool samples obtained from the Integrative Human Microbiome
#' Project Consortium (iHMP). Group label: CD and Controls.
#'
#' @format a [`phyloseq::phyloseq-class`] object
#' @source \url{https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/ibd_data.zip}
#' @name data-pediatric_ibd
#' @aliases pediatric_ibd
#' @docType data
#' @author Yang Cao
"pediatric_ibd"

#' Oxygen availability 16S dataset, of which taxa table has been summarized for
#' python lefse input
#'
#' A small subset of the HMP 16S dataset for finding biomarkers characterizing
#' different level of oxygen availability in different bodysites
#'
#' @format a [`phyloseq::phyloseq-class`] object
#' @source \url{http://huttenhower.sph.harvard.edu/webfm_send/129}
#' @name data-oxygen
#' @aliases oxygen
#' @docType data
#' @author Yang Cao
"oxygen"


