#' Regulatory links verified by research and experimentation, downloaded from YEASTRACT database
#' First column represents the regulator genes and the second column represents the target genes
#' Each row represents a regulatory relationship validated by experimentation and curated in the YEASTRACT database


#' @docType data
#'
#' @usage data(yeasts_arcs_validated)
#'
#' @format An object of class data.frame. First column represents the regulator genes and the second column represents the target genes
#'
#' @keywords datasets
#'
#' @references Monteiro et al. (2020)  Nucleic Acids Research, 48(D1):D642-D649

#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC25624/})
#' @source \href{http://www.yeastract.com/index.php}
#'
#' @examples
#' data(yeasts_arcs_validated)
"yeasts_arcs_validated"
