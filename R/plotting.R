#' Generates the Manhattan plot for one set of results
#' 
#' One gene, one type of analysis
#'
#' @export
#' @import graphics
#' @param results.df dataframe containing results, needs SNP p-values and positions
#' @param gene.name gene name for which we want the Manhattan 
#' @param gene.info.file hdf5 file containing gene information
#' @param snp.name SNP of interest within the Manhattan plot (optional)
OneManhattan <- function(results.df, gene.name, gene.info.file,snp.name='') {
  res <- MakeManhattanDataFrame(results.df, gene.name, gene.info.file)
  chrom <- res[1,]$chrom
  # get best p-value to adjust ylim 
  best.pv <- min(res$p_value)
  y.max <- -log10(best.pv)
  bottom.rect <- -y.max/4
  top.rect <- bottom.rect + y.max/7
  # plot Mannhattan
  plot(res$pos, -log10(res$p_value), main = paste0(gene.name,", ",snp.name), 
       xlab = paste0("position on chromosome ",chrom), ylab = "-log10(pvalue)",
       frame.plot=FALSE, cex = 0.6, ylim = c(bottom.rect, y.max))
  # add axes
  axis(side = 2)
  snp.pos <- res[res$gdid == snp.name,]$pos
  # add line and diamond shape at SNP of interest
  points(x = snp.pos, -log10(res[res$pos == snp.pos,]$p_value), pch = 23, 
         cex = 1.8, col = "firebrick", bg="red")
  abline(v = snp.pos, col = "firebrick")
  segments(x0 = snp.pos, y0 = 0, x1 = snp.pos, y1 = 30, col = "firebrick")
  # add gene position
  gene.start <- res$gene_start
  gene.end <- res$gene_end
  rect(gene.start, bottom.rect, gene.end, top.rect, col = "cornflowerblue", border = "cornflowerblue")
}

#' Make dataframe for Manhattan plots
#' 
#' @export 
#' @param results.df dataframe containing results for all genes
#' @param gene.name gene name for which we want the Manhattan 
#' @param gene.info.file hdf5 file containing gene information
MakeManhattanDataFrame <- function(results.df, gene.name, gene.info.file) {
  res <- results.df[results.df$assoc_gene == gene.name,]
  gene.pos <- GetGenePosition(gene.name, gene.info.file)
  res$gene_start <-  gene.pos[[1]]
  res$gene_end <-  gene.pos[[2]]
  results.df <- as.data.frame(res)
  results.df
}

#' Get gene starting and ending positions
#' 
#' add gene start and gene end to dataframe for Manhattan plots
#' 
#' @export 
#' @import rhdf5
#' @param gene.name gene name for which we want the Manhattan 
#' @param gene.info.file hdf5 file containing gene information   
GetGenePosition <- function(gene.name, gene.info.file) {
  gene_info <- rhdf5::h5read(gene.info.file, name = "gene_info")
  geneID <- rhdf5::h5read(gene.info.file, name = "geneID")
  gene.start <- gene_info$start_position[geneID==gene.name]
  gene.end <- gene_info$end_position[geneID==gene.name]
  H5close()
  list(gene.start,gene.end)
}

#' Generates overlapping Manhattan plots for two or more sets of results
#' 
#' One gene, comparison between two or more types of analysis
#'
#' @export
#' @param results.dfs list of dataframes containing results, needs SNP p-values and positions
#' @param gene.name gene name for which we want the Manhattan 
#' @param gene.info.file hdf5 file containing gene information
#' @param snp.name SNP of interest within the Manhattan plot (optional)
TwoManhattan <- function(results.dfs, gene.name, gene.info.file,snp.name='') {
  res1 <- MakeManhattanDataFrame(results.dfs[[1]], gene.name, gene.info.file)
  res2 <- MakeManhattanDataFrame(results.dfs[[2]], gene.name, gene.info.file)
  chrom <- res1[1,]$chrom
  # get best p-value to adjust ylim 
  best.pv <- min(min(res1$p_value),min(res2$p_value))
  y.max <- -log10(best.pv)
  bottom.rect <- -y.max/4
  top.rect <- bottom.rect + y.max/7
  # plot first Mannhattan
  plot(res1$pos, -log10(res1$p_value), main = paste0(gene.name,", ",snp.name), 
       xlab = paste0("position on chromosome ",chrom), ylab = "-log10(pvalue)",
       frame.plot=FALSE, cex = 0.6, ylim = c(bottom.rect, y.max), col = "cornflowerblue")
  # plot second Mannhattan
  points(res2$pos, -log10(res2$p_value),col = "black",cex=0.6)
  # add axes
  axis(side = 2)
  snp.pos <- res1[res1$gdid == snp.name,]$pos
  # add line and diamond shape at SNP of interest
  points(x = snp.pos, -log10(res1[res1$pos == snp.pos,]$p_value), pch = 23, 
         cex = 1.8, col = "firebrick", bg="red")
  segments(x0 = snp.pos, y0 = -0.2, x1 = snp.pos, y1 = y.max+y.max*0.2, col = "firebrick")
  # add gene position
  gene.start <- res1$gene_start
  gene.end <- res1$gene_end
  rect(gene.start, bottom.rect, gene.end, top.rect, col = "forestgreen", border = "forestgreen")
}
#' Generates the qq plot for one set of results
#' All genes, one type of analysis
#' @export
#' @import graphics
#' @param results.df dataframe containing results, needs SNP p-values and permutations
PlotQQ <- function(results.df) {
  # plot qqplot
  x = sort(-log10(runif(dim(results.df)[1],min=0,max=1)))
  y1 = sort(-log10(results.df$p_value))
  plot(x, y1, xlab = "-log10(expected pvalues)", ylab = "-log10(observed pvalues)",
       frame.plot=FALSE, cex = 0.6, col = "cornflowerblue")
  # add axes
  axis(side = 2)
  # add permutation qqplot
  y2 = sort(-log10(results.df$permutation_0))
  points(x , y2, cex = 0.6, col = "black")
  lines(x = c(0,7), y = c(0,7),col='firebrick')
}