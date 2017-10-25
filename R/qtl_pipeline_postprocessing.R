#' Get results from qtl pipeline and summarize results 
#' 
#' All genes, one scan
#'
#' @export
#' @import rhdf5
#' @import qvalue
#' @import dplyr
#' @param results.folder full path to folder where results are stored
#' @param multiple.testing.global method, can be BF (Bonferroni) or ST (Storey) 
GetResults <- function(results.folder, multiple.testing.global = "ST") {
  observed.features <- 0
  results <- NULL
  files.to.read <- list.files(results.folder,pattern = "qtl_results.*h5", full.names = T)
  for ( i in files.to.read ) {
    base.name <- gsub(gsub(i, pattern = ".h5", replacement = ""), pattern = "./qtl_results", replacement = "")
    if ( length(list.files("./",pattern = base.name)) == 3 | length(list.files("./",pattern = base.name)) == 4 ) {
      tmp <- h5dump(file = i)
    } 
    if ( length(tmp) > 0 ) {
      for ( j in names(tmp) ) { tmp[[j]][["feature"]] <- j }
      observed.features = observed.features + length(tmp)
      df <- bind_rows(tmp)
      if ( nrow(df) > 0 ) { results = rbind(results,df) }
    }
  }
  H5close()
  colnames(results)[ which(colnames(results) == "corr_p_value") ] <- "feature_corr_p_value"
  if ( length(which(is.na(results$feature_corr_p_value))) != 0 ) {
    results <- results[-which(is.na(results$feature_corr_p_value)),]
  }
  ##Multiple testing
  if ( multiple.testing.global == "ST" ) {
    results["global_corr_p_value"] <- qvalue(results$feature_corr_p_value)$qvalues
  } else if ( multiple.testing.global == "BF" ) {
    results["global_corr_p_value"] <- results$feature_corr_p_value*observed.features
    results$global_corr_p_value[results$global_corr_p_value > 1] <- 1
  }
  results <- results[order(results$global_corr_p_value, decreasing = F),]
  results
}

#' Get permutations from qtl pipeline for plotting purposes 
#' @export
#' @param results.folder full path to folder where results are stored
GetPerms <- function(results.folder) {
  observed.features <- 0
  perms <- NULL
  files.to.read <- list.files(results.folder,pattern = "perm_results.*h5", full.names = T)
  for ( i in files.to.read ) {
    base.name <- gsub(gsub(i, pattern = ".h5", replacement = ""), pattern = "./perm_results", replacement = "")
    if ( length(list.files("./",pattern = base.name)) == 3 | length(list.files("./",pattern = base.name)) == 4 ) {
      tmp <- h5dump(file = i)
    } 
    if ( length(tmp) > 0 ) {
      for (j in names(tmp)) { tmp[[j]][["feature"]] <- j }
      observed.features = observed.features + length(tmp)
      df <- bind_rows(tmp)
      if ( nrow(df) > 0 ) { perms = rbind(perms,df) }
    }
  }
  H5close()
  perms
}

#' Merge results with one permutation - for plotting
#' @export
#' @param results.folder full path to folder where results are stored
#' @param multiple.testing.global method, can be BF (Bonferroni) or ST (Storey) 
MergePerms <- function(results.folder, multiple.testing.global = "ST"){
  results.df <- GetResults(results.folder, multiple.testing.global = "ST")
  perms.df <- GetPerms(results.folder)
  results.with.permutations <- inner_join(results.df, perms.df, by = c("feature", "snp_id"))
  results.with.permutations <- results.with.permutations[order(results.with.permutations$global_corr_p_value),]
  results.with.permutations
}

#' Merge results with gene info & snp info
#' @export
#' @param results.folder full path to folder where results are stored
#' @param gene.names.df dataframe containing snp info, e.g. rsID
#' @param snp.info.df dataframe containing snp info, e.g. rsID
#' @param multiple.testing.global method, can be BF (Bonferroni) or ST (Storey) 
MergeGeneInfo <- function(results.folder, gene.names.df, snp.info.df, multiple.testing.global = "ST"){
  results.perm <- MergePerms(results.folder, multiple.testing.global = "ST", significance.threshold = 0.1)
  results.genes <- left_join(gene.names.df, results.perm, by = "feature")
  results.snp.info <- inner_join(results.genes, snp.info.df, by = c("chrom", "pos","assoc_gene"))
  results.snp.info
}

#' Only save top lead per gene - optionally only significant 
#' @export
#' @param results.folder full path to folder where results are stored
#' @param gene.names.df dataframe containing snp info, e.g. rsID
#' @param snp.info.df dataframe containing snp info, e.g. rsID
#' @param multiple.testing.global method, can be BF (Bonferroni) or ST (Storey) 
#' @param significance.threshold at global level (optional)
LeadSnpsOnly <- function(results.folder, gene.names.df, snp.info.df, multiple.testing.global = "ST", significance.threshold = 1){
  results <- MergeGeneInfo(results.folder, gene.names.df, snp.info.df, multiple.testing.global = "ST")
  results <- results[order(results$global_corr_p_value),]
  results <- results[which(results$global_corr_p_value < significance.threshold),]
  res <- results[-which(duplicated(results$feature)),]
  res
}

