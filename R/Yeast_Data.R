#' Spellmans Yeast Cell Cycle data
#' From this study "Comprehensive Identification of Cell Cycle–regulated Genes of the Yeast Saccharomyces cerevisiae by Microarray HybridizationD"
#' by Spellman et al.
#' This microarray dataset is a comprehensive catalogue of yeast genes whose transcript levels
#' change periodically withinthe cell cycle. We use σ-factor based synchronisation data that containsexpression ofS. cerevisiaegenes sampled over 18 time points.
#' We filter out 786 genes from the entire genome, as those genes that have demonstrated consistent oscillating changes within the yeast cell cycle in Spellman et al.(1998).
#' This dataset has been copied from the Supplementary Materials of Lebre et al. paper "Inferring dynamic genetic networks with low order
#' independencies"

#' @docType data
#'
#' @usage data(yeast_data)
#'
#' @format An object of class data.frame. Each column represents a gene and each row represents the samples of gene expression measured over time.
#'
#'
#' @keywords datasets
#'
#' @references Spellman et al. (1998) Mol Biol Cell. 1998 Dec; 9(12): 3273–3297.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC25624/})
#' @references Lèbre et al. (2009) Statistical Applications in Genetics and Molecular Biology,8.
#' @source \href{http://icube-bfo.unistra.fr/fr/index.php/G1DBN_supplementary_material}
#'
#' @examples
#' data(yeast_data)
"yeast_data"
