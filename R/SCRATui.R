#' SCRATui
#' 
#' Launch the SCRAT user interface in local machine
#'
#' This function will automatically launch the SCRAT user interface in a web browser. 
#' The user interface can also be accessed by http://zhiji.shinyapps.io/SCRAT. Neither R nor any packages are required in this online version.
#' However, it is highly recommended that the user interface be launched locally for faster running speed.
#' 
#' @export
#' @import shiny GenomicAlignments ggplot2 reshape2 pheatmap scatterD3 
#' @author Zhicheng Ji, Weiqiang Zhou, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' \dontrun{
#'    SCRATui()
#' }

SCRATui <- function() {
      shiny::runApp(system.file("shiny",package="SCRAT"))
}
