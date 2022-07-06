#' @name differentialDeviations_for_ArchR
#' Function to see whether deviations differ between groups
#' @param object ArchRProject object
#' @param groups column in colData of object 
#' with group information
#' @param alternative only used if there are two groups -- two.sided or one 
#' sided test
#' @param parametric use parametric test. alternatively will use kruskal wallace
#' @return data.frame with p value and adjusted p value
#' @export
#' @author Alicia Schep - all chromVAR code 
#' @author Angelo Pelonero - adapted to ArchR
#' @importFrom stats aggregate
#' @songs_used_for_dev "Gravy Train" by Lettuce & "Naptime" by Mark Lettieri

differentialDeviations_ArchR <- function(object, 
                                         groups,
                                         alternative = c("two.sided", "less", "greater"),
                                         parametric = TRUE) {
  stopifnot(is(object,"ArchRProject"))
  if (length(groups) == 1 && groups %in% colnames(ArchR_Project)) {
    groups_archr <- getCellColData(ArchRProj = ArchR_Project, groups, drop = FALSE)
  } else {
    print("'groups' found in ArchR Metadata?")
    print(table(groups %in% colnames(ArchR_Project)))
    stop("LAZY DEV ERROR: Please, please.. one group at a time. Support for ngroups > 1 may be added in future. \n Additionally, please ensure 'groups' %in% colnames(ArchR_Project) == TRUE")
  }
  
  ### force groups from ArchR into ordered factor vetor
  ## limited to one group per function call by use of index
  ## TODO by me or you: add support for multiple groups e.g. `for(i in length(groups)){ ~~~loop all this and return as list of group-named dataframes~~~ }`
  groups <- as.factor(groups_archr[,1])
  
  alternative <- match.arg(alternative)
  ## get deviations matrix from ArchR project
  if("MotifMatrix" %in% getAvailableMatrices(object)){
    archrMotifMatrix <- getMatrixFromProject(object, useMatrix = "MotifMatrix")
    inputs <- as.matrix(archrMotifMatrix@assays@data@listData$deviations)
  }
  else {
    errorCondition("MotifMatrix not found in ArchR Project - please run addDeviationsMatrix() before running stats on.. deviations.")
  }
  
  ### original chromVAR code:
  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      p_val <- apply(inputs, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(inputs, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(inputs, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(inputs, 1, kw_helper, groups)
    }
  }
  
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}


t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]],splitx[[2]],
                alternative = alternative, 
                paired = FALSE,
                var.equal = FALSE)$p.value)
}

anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- oneway.test(devs ~ groups, tmpdf, var.equal = FALSE)
  return(res$p.value)
}

kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
                     alternative = alternative, 
                     paired = FALSE)$p.value)
}

#' @name differentialVariability_for_ArchR
#'
#' Function to determine whether groups differ in variability
#' @param object ArchRProject object
#' @param groups name of column in colData of object
#' with group information
#' @param parametric use parametric test. alternatively will use kruskal wallace
#' @return data.frame with p value and adjusted p value
#' @author Alicia Schep - all chromVAR code 
#' @author Angelo Pelonero - adapted to ArchR
#' @export

differentialVariability_ArchR <- function(object, groups, parametric = TRUE) {
  stopifnot(is(object,"ArchRProject"))
  if (length(groups) == 1 && groups %in% colnames(ArchR_Project)) {
    groups_archr <- getCellColData(ArchRProj = ArchR_Project, groups, drop = FALSE)
  } else {
    print("'groups' found in ArchR Metadata?")
    print(table(groups %in% colnames(ArchR_Project)))
    stop("LAZY DEV ERROR: Please, please.. one group at a time. Support for ngroups > 1 may be added in future. \n Additionally, please ensure 'groups' %in% colnames(ArchR_Project) == TRUE")
  }

  ### force groups from ArchR into ordered factor vetor
  ## limited to one group per function call by use of index
  ## TODO by me or you: add support for multiple groups e.g. `for(i in length(groups)){ ~~~loop all this and return as list of group-named dataframes~~~ }`
  groups <- as.factor(groups_archr[,1])
  
  ## get deviations matrix from ArchR project
  if("MotifMatrix" %in% getAvailableMatrices(object)){
    archrMotifMatrix <- getMatrixFromProject(object, useMatrix = "MotifMatrix")
    inputs <- as.matrix(archrMotifMatrix@assays@data@listData$deviations)
  }
  else {
    errorCondition("MotifMatrix not found in ArchR Project - please run addDeviationsMatrix() before running stats on.. deviations.")
  }
  
  ### original chromVAR code:
  # Brown-Forsythe test
  if (parametric) {
    p_val <- apply(inputs, 1, bf_var_test, groups)
  } else {
    p_val <- apply(inputs, 1, bf_kw_var_test, groups)
  }
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}

bf_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(anova(lm(median_diff ~ groups))[1, 5])
}

bf_kw_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(kruskal.test(median_diff ~ groups)$p.value)
}