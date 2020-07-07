#' This function calculates the first order conditional independence score using the G1DBN package
#' @param data data.frame or a matrix with n rows (=time points) and p columns (=genes) containing the gene expression time series.
#' @param predPosition To be specified to reduce the set of possible predictor genes to a subset of \eqn{d<p} genes: an array included in \code{[1,p]} defining the position of the d predictor genes in the data matrix (n * p), default=NULL.
#' @return A matrix of first order conditional independence scores for each edge in the network
#' @example
#' sim <- SimulateData(genes = 50, timepoints = 20, seed = 1, prop = 0.25)
#' score_mat <- ExecuteG1DBNS1(sim$data)
#' @export
ExecuteG1DBNS1 <- function(data, predPosition = NULL)
{
  S1 <- G1DBN::DBNScoreStep1(data, method='ls', predPosition = predPosition)
  return(S1)
}

#' Simulate Data
#' This function first generates a random network and then simulates the multi variate time series data according to the following  first order Auto-Regressive process,
#' \eqn{X(t) = AX(t-1) + B + E(t)},
#' where \eqn{E(t)} follows a zero-centered multivariate gaussian distribution whose variance matrix S is diagonal.
#' It uses different functions from the G1DBN package.
#' @param genes Number of genes or variables
#' @param  timepoints Number of samples
#' @param seed Random seed,  is a number used to initialize a pseudorandom number generator.
#' @param prop The proportion of edges in the entire network
#' @param range vector with 4 elements specifying range values for the adjacency matrix generation (minimum negative value, maximum negative value, minimum positive value, maximum positive value)
#' @param errors vector of 2 elements specifying min and max value of the uniform distribution from which the error variances are drawn
#' @return A list comprising of dataset generated (data) and the network (RealNet) from which the data is generated.
#' @examples
#' sim <- SimulateData(genes = 50, timepoints = 20, seed = 1, prop = 0.25)
#' head(sim$data)
#' @export
#' @import G1DBN
#' @import stats

SimulateData <- function(genes, timepoints, seed = 1, prop = 0.05, range = c(-0.95,-0.05,0.05,0.95), errors = c(0.03, 0.08) )
{
  n = timepoints
  p = genes
  set.seed(seed)
  prop = prop
  MyNet <- G1DBN::SimulNetworkAdjMatrix(p,prop, range)

  ## initializing the B vector
  ## initializing the variance of the noise
  sigmaEps <- runif(p,errors[1], errors[2])
  B1 <- runif(p,range[1], range[2])
  B2 <- runif(p,range[3], range[4])
  B = unlist(lapply(1:p, function (x) sample(c(B1[x], B2[x]), 1)))
  ## initializing the process Xt
  X0 <- B +  runif(p,errors[1], errors[2])
  Xn <- G1DBN::SimulGeneExpressionAR1(MyNet$A,B,X0,sigmaEps,n)

  ml = list("data" = Xn, "RealNet"  = MyNet)
  return (ml)
}

#' Puts the data in a time frame with values of variables in the current time slice and their immediate prior value from the previous time slice.
#' @param data a data.frame, columns represent the variables and rows represent the temporal data
#' @return a data.frame with twice the number of columns.
#' @keywords internal

ShiftData <- function(data)
{

  n = nrow(data)
  p = ncol(data)
  data2 = as.data.frame(matrix(NA, ncol = 2*p, nrow = n-1))
  data2[1:n-1, 1:p] = data[1:n-1,]
  new_names = vector(length = p)
  colnames(data) = paste0("V", 1:p)
  for(i in 1:ncol(data))
  {
    new_names[i] = paste0(colnames(data)[i] , "_1")


    data2[,i+ncol(data)] = data[2:n, i]
  }
  colnames(data2) =  c(colnames(data), new_names)
  data2
}



#' This function applies L1 regularised Lasso regression on the dataset. The regularisation parameter is selected by cross validation.
#' It uses function from package "lars".
#' @param data a data.frame or a matrix with n rows (=time points) and p columns (=genes) containing the gene expression time series.
#' @param pred_pos indices of the potential predictor variables in the dataset. If the value is NULL, all variables are considered as potential predictors.
#' @return a matrix of Lasso regression coefficients.
#' @example
#' sim <- SimulateData(genes = 50, timepoints = 30, seed = 1, prop = 0.25)
#' LASSO <- ApplyLasso(sim$data, c(1,2,3,5))
#' @export


ApplyLasso <- function(data, pred_pos = NULL)
{
  set.seed(1)
  if(is.null(pred_pos))
  {
    #print("No predictors")
    pred_pos = 1:ncol(data)
  }
  predictors  = as.matrix(data[1:(nrow(data)-1), pred_pos])
  response = data[2:nrow(data),]
  nR = ncol(response)
  arcs = lapply(1:nR, function(y)
  {
    #print(y)
    Y = response[,y]
    arcs = LASSO(X = predictors, Y = Y)

  })
  scores = abs(do.call(rbind, arcs))
  return(scores)
}

#' Internal function. Applies Lasso regression on the dataset using the package lars
#' @param X Regressors
#' @param Y Response variable
#' @return Lasso regression coefficients for each outcome variable Y
#' @keywords internal
#' @import lars
LASSO <- function(X,Y, use.Gram = TRUE)
{
cv <- lars::cv.lars(x = X, y = Y, plot.it = FALSE,use.Gram = use.Gram)
ideal_l1_ratio <- cv$index[which.max(cv$cv - cv$cv.error <= min(cv$cv))]
obj <- lars::lars(x = X, y = Y, type = "lasso", use.Gram = use.Gram)
scaled_coefs <- scale(obj$beta, FALSE, 1 / obj$normx)
l1 <- apply(X = scaled_coefs, MARGIN = 1, FUN = function(x) sum(abs(x)))
c1 = coef(obj)[which.max(l1 / tail(l1, 1) > ideal_l1_ratio),]
return(c1)
}


#' Normalise the data values between 0 and 1.
#' @param x vector or a matrix of unscaled values
#' @return a vector or a matrix with values normalized within 0-1 range
#' @example
#' sim <- SimulateData(genes = 50, timepoints = 20, seed = 1, prop = 0.25)
#' lars_coeffs <- Normalisation_01(ApplyLasso(sim$data, c(1,2,3,4,5)))
#' @export
Normalisation_01<- function(x)
{
  (x-min(x))/(max(x)-min(x))
}


#' This function applied greedy hill climbing search to find the parents of each node in the network from the previous time slice.
#' The score for all arcs absent in the network is calculated as :
#' @param data a data.frame where each column represents the variables and rows contain the temporal data
#' @param score_mat a matrix of scores of edges to be added to the BIC score. This can be the matrix of lasso regression coefficience
#' or the p-value of the first order conditional independence tests.
#' @param maxP maximum number of allowed parents per node
#' @param gamma value of the gamma constant added to the final score, to move the very small edge score values away from zero
#' @param pred_pos position of the putative regulators in the dataset. If this value is NULL, all variables are considered as putative regulators.
#' @param score This refers to the type of score to be calculated. The default is "Score_LOPC".
#' Other possible values are Score_LASSO and BIC.
#' @return bnlearn object
#' @export
#' @import bnlearn
#' @example
#' sim <- SimulateData(genes = 30, timepoints = 10, seed = 1, prop = 0.25)
#' score_mat <- ExecuteG1DBNS1(sim$data)
#' dbn <- GreedySearchParentsWithPriors(sim$data, score_mat$S1ls, score = "Score_LOPC", gamma = 0.2, maxP = 2)

GreedySearchParentsWithPriors <- function(data,  score_mat = NULL, maxP = Inf, gamma = 0, pred_pos = NULL, score = "Score_BIC")
{
  p <- ncol(data)
  n <- nrow(data)
  dataS <- ShiftData(data)
  if(is.null(score_mat))
  {
    score_mat <- matrix(0,p,p)
  }

  if(is.null(pred_pos))
  {
    all_parents <- names(dataS)[1:p]#[!grepl("_1",  names(dataS))]
  }

  else
  {
    all_parents = names(dataS)[pred_pos]
  }

  ## Check score arg ##

  if(tolower(score) == "score_lasso")
  {
    print("Chosen score: Score_Lasso")
    score_mat <- Normalisation_01(score_mat) #Normalise between 0 and 1
    score_mat <- 1 - score_mat
  }

  else if(tolower(score) == "score_lopc")
  {
    print("Chosen score: Score_Lopc")
  }
  else if(tolower(score) == "bic")
  {
    print("Chosen score: BIC")
  }
  else
  {
    stop("In valid score argument. Possible inputs are 'Score_LASSO', 'BIC' or 'Score_LOPC'")
  }
  targets <- 1:p

  local_bns_arcs <- lapply(targets, function(x)

    GreedySearchWithPriors_OneTarget(score_mat = score_mat, targetIdx = x, all_parents = all_parents, dataS = dataS, maxP = maxP,  gamma = gamma, score = score))

  fullBN <- do.call(what = rbind, args = local_bns_arcs)
  nodes <- names(dataS)
  bn<- bnlearn::empty.graph(nodes)
  arcs(bn) <- fullBN

  return (bn)


}
#' Internal function
#' Applied greedy hill-climbing search to each local target node in the netork
#' \eqn{E_ij = log(score_mat + gamma)}
#'  \eqn{¬E_ij =  log(1-score_mat+gamma)}
#' @return set of parents arcs for each local target node
#' @keywords internal
#' @import bnlearn
GreedySearchWithPriors_OneTarget<- function(score_mat, targetIdx, all_parents, dataS, maxP = maxP, gamma, score)
{
  keep_adding = 1
  targetNode <- paste0("V", targetIdx , "_1")
 # all_parents_names <- paste0("V", all_parents)
  all_parents_names <- all_parents
  nodes <- c(all_parents_names, targetNode)
 # print(nodes)
  print(paste0("Currently finding parents for: ", targetNode))
  data_subset <- dataS[,nodes]

  if(tolower(score) == "score_lasso")
  {
    weights_absent <- log(score_mat + gamma)
    weights_present <- log(1 - score_mat+gamma)
  }

  else if(tolower(score) == "score_lopc")
  {
  #  print("LOPC")
    weights_absent <- log(score_mat + gamma)
    weights_present <- log(1 - score_mat+gamma)
  }
  else if(tolower(score) == "bic")
  {
    #score_mat <- matrix(0, p, p)
    weights_absent = score_mat * 0
    weights_present = score_mat * 0
  }

  weight_vector_absent <-  weights_absent[targetIdx, ]
  weight_vector_present <- weights_present[targetIdx,]
  names(weight_vector_absent) <- all_parents_names
  names(weight_vector_present) <- all_parents_names
  bn <- bnlearn::empty.graph(nodes)

  old_bic <-  bnlearn::score(x  =bn, data = data_subset, type = "bic-g",by.node = TRUE )[targetNode]
  empty_g1dbn <- sum(weight_vector_absent)
  old_bic <- old_bic+empty_g1dbn
  TotalParents <- 0
  while (keep_adding > 0 && TotalParents < maxP)
  {

    new_bics <- unlist(lapply(all_parents, function(x)
    {

      if(!is.na(x))
      {

        bn_1 <- bnlearn::set.arc(bn, from = x, to = targetNode)
        score <- bnlearn::score(x  = bn_1, data = data_subset, type = "bic-g", by.node = TRUE )[targetNode]

        weights_present <- weight_vector_present[x]
        weights_absent <- sum(weight_vector_absent[which(names(weight_vector_absent) != x)])
      }
      else
      {
        score <- NA
      }


      sum(score, weights_present, weights_absent)
    }))


    bic_and_weight <- new_bics #+ weight_vector
  # print(sum(new_bics))
    difference <- bic_and_weight - old_bic
    diff <- 0
    keep_adding <- length(which(difference > diff)) > 0

    if(keep_adding)
    {
    #  print(difference)
      max_idx <- which.max(difference)
      best_parent <- all_parents_names[max_idx]
      bn <- bnlearn::set.arc(bn, from = best_parent, to = targetNode)
      TotalParents <- TotalParents + 1
      print(paste0("Adding Parent: ", best_parent ))
      all_parents[max_idx] <- NA
      old_bic <- bic_and_weight[max_idx]
    }

  }

  return (bnlearn::arcs(bn))

}

#' Take a learned Bayesian network structure (bnlearn object) and
#' returns a data frame in which each row represents a pair of nodes conected in the network.
#' This function is used to retrieve the arcs from the network.
#' @param bn Bnlearn object Bayesian network
#' @param data data.frame from which the Bayesian network structure is learned.
#' @return arc_set a data.frame with two columns. First column shows the regulators, second column shows the targets.
#' @keywords internal
#'
#' @example
#' sim <- SimulateData(genes = 20, timepoints = 10, seed = 1, prop = 0.25)
#' score_mat <- ExecuteG1DBNS1(sim$data)
#' dbn <- GreedySearchParentsWithPriors(sim$data, score_mat$S1ls, score = "Score_LOPC", gamma = 0.2, maxP = 2)
#' arcs <- MarkArcSet(dbn)
#' @export


MakeArcSet <- function(bn, data)
{
  arc_set <- bnlearn::arcs(bn)
  froms <- unlist(lapply(arc_set[,1], function(from) as.numeric(gsub("V", from, replacement = ""))))
  tos <- unlist(lapply(arc_set[,2], function(to) {

    to <- gsub("V", to, replacement = "")
    to <- as.numeric(gsub("_1",to, replacement = "" ))}))

  froms <- names(data)[froms]
  tos <- names(data)[tos]
  arc_set <- as.data.frame(cbind("from" = froms, "to" = tos ))
  return(arc_set)

}

#' Exhaustively looks for best combination of parents, up to max 3 parents per node.
#' Currently limits max number of parents to 3. Will be changed in later uppdates.
#' Please note that this function currently allows you to limited the maximum number of parents to 3
#' It does not allow you to define a set of possible parents as in the GreedySearchWithPriors function.
#' @param data a data.frame/matrix with n rows (=time points) and p columns (=genes) containing the gene expression time series.
#' @param type A scoring objective. Possible values are "Score_BIC", "Score_LASSO" and "Score_LOPC".
#' @param gamma a constant parameter to be added to the score_mat, to shift the values away from zero. Required for Score_LASSO and Score_LOPC
#' @param score_mat a weight matrix to be added to the BIC score, to calculate "Score_LOPC" and "Score_LASSO"
#' @export
#' @example
#' sim_data <-  SimulateData(genes = 50, timepoints = 20, seed = 1, prop = 0.05 )
#' lasso_mat = ApplyLasso(sim_data$data)
#' dbn_score_lasso <- ExhaustiveSearchForBestParents(sim_data$data, type = "Score_LASSO",
#' gamma = 0.0, score_mat = lasso_mat)
#' @example
#' sim_data <-  SimulateData(genes = 50, timepoints = 20, seed = 1, prop = 0.05, range = c(-0.95,-0.05,0.05,0.95), errors = c(0.03, 0.08) )
#' dbn_score_bic <-  ExhaustiveSearchForBestParents(sim_data$data, type = "Score_BIC")
#' @example
#' sim <-  SimulateData(genes = 10, timepoints = 8, seed = 1, prop = 0.25, range = c(-0.95,-0.05,0.05,0.95), errors = c(0.03, 0.08) )
#' lopc <- ExecuteG1DBNS1(sim$data)
#' dbn_score_lopc <- ExhaustiveSearchForBestParents(sim$data, type = "Score_LOPC", gamma = 0.1, score_mat = lopc$S1ls)


ExhaustiveSearchForBestParents <- function(data, type = "Score_BIC", gamma = 0.0, score_mat = NULL)
{
  #print(score_mat[1:4,1:4])
  p <- ncol(data)
  n <- nrow(data)
  dataS <- ShiftData(data)
  allNodes <- names(dataS)
  All_possible_parents <-  names(dataS)[1:p]
 # print(2)
  Allcombinations = MakeAllPossibleCombination(All_possible_parents)
  All_possible_targets = allNodes[grepl("_1", allNodes)]
  #print(3)
  if(type == "Score_LOPC")
  {
  #  G1_mat <- ExecuteG1DBNS1(data)
  #  score_mat <- G1_mat$S1ls
    priors <-  lapply(All_possible_targets, function(x) CalculatePriors(allNodes, x, Allcombinations, score_mat, gamma))
    priors <- do.call(rbind,priors)

  }
  else if(type == "Score_LASSO")
  {
 #   lasso <- ApplyLars(data)
    lasso <- score_mat
    lasso <- Normalisation_01(lasso) #Normalise between 0 and 1
    lasso <- 1-lasso #Important to invert
    score_mat <- lasso
    priors <- lapply(All_possible_targets, function(x) CalculatePriors(allNodes, x, Allcombinations, score_mat, gamma))
    priors <- do.call(rbind, priors)
  }
  else
  {
    #set priors to zero, Only bic scores to be calculated
    priors <- matrix(data = 0, ncol = length(Allcombinations), nrow = length(All_possible_targets) )
    gamma <- 0

  }
  bics <-  lapply(All_possible_targets, function(x) ScoreEntirepossibleSetsOfParents(dataS, x, Allcombinations))
  bics <- do.call(rbind, bics)
  bics <- CombineBICandPriors(bics, priors)
  BN <- ConstructBNUsingMaxScoringParents(bics, Allcombinations, allNodes)
  return(BN)
}

#' Makes all possible combination of sets of parents, maximum number of parents set to 3
#' @keywords internal
#' @import utils
MakeAllPossibleCombination <- function(all_parents)
{
  combinations0 = combn(all_parents, 0)
  combinations1 = combn(all_parents, 1)
  combinations2 = combn(all_parents,2)
  combinations3 = combn(all_parents, 3)

  parents_0 <- c(NA)
  parents_1 <- lapply(1:ncol(combinations1), function(x) return(c(combinations1[x])))
  parents_2 <- lapply(1:ncol(combinations2), function(x) return(c(combinations2[1,x], combinations2[2,x])))
  parents_3 <- lapply(1:ncol(combinations3), function(x) return(c(combinations3[1,x], combinations3[2,x], combinations3[3,x])))

  all_parents <- c(parents_0, parents_1, parents_2, parents_3)
  return (all_parents)
}


#' This function takes in all scores for all parents combinations (columns) of all target nodes (rows)
#' and returns the maximum scoring network.
#' @param BICandWeigh Score_LASSO or Score_LOPC for each combination of parents for each target node
#' @param all_parent_combinations All possible combinations of parents , upto 3 parents max
#' @param allNodes all the target nodes in the network
#' @keywords internal
#' @import bnlearn
ConstructBNUsingMaxScoringParents <- function(BICandWeight,  all_parent_combinations, allNodes)
{
  max_ids <- apply(BICandWeight, 1, which.max)
  #allNodes = names(dataS)
  All_possible_targets = allNodes[grepl("_1", allNodes)]
  max_parents <- lapply(max_ids, function(x) unlist(all_parent_combinations[x]))
  bn <- bnlearn::empty.graph(allNodes)
  for (i in 1:length(All_possible_targets))

  {
    #print(i)
    target <- All_possible_targets[[i]]
    parents <- max_parents[[i]]
    print(parents)
    if(!is.na(parents))
    {
      #print(length(parents))

      for (par in parents)
      {
        bn <- bnlearn::set.arc(bn, from = par, to = target)
      }
    }
  }
  return (bn)
}


#' This methods calculates BIC score for each target node, for all the possible combinations of parents
#' @keywords internal
#' @export
ScoreEntirepossibleSetsOfParents <- function(dataS, tarNode, all_combinations)
{

  tarIndex = sub("V", "", tarNode)
  print(tarNode)
  tarIndex = as.numeric(sub("_1", "", tarIndex))
  allScores <- array(data = NA, dim = length(all_combinations))
  #  beta = beta[which(beta$to == tarNode),]
  allNodes = names(dataS)
  All_possible_parents = allNodes[!grepl("_1", allNodes)]

  bn = bnlearn::empty.graph(nodes = c(All_possible_parents, tarNode))
  ScoreOfGraph = bnlearn::score(bn, dataS[,c(bnlearn::nodes(bn))], type = "bic-g", by.node = TRUE)[tarNode]

  allScores[1] = ScoreOfGraph
  print(length(all_combinations))
  for (i in 2:length(all_combinations))
  {
   # print(i)
    bn = bnlearn::empty.graph(nodes = c(All_possible_parents, tarNode))

    parents = all_combinations[[i]]
    for (j in 1:length(parents))
    {
      if(!is.na(parents[j]))
      {
        bn <- bnlearn::set.arc(bn, from = parents[j], to = tarNode)
      }

    }

    newScore = bnlearn::score(bn, dataS[,c(bnlearn::nodes(bn))], type = "bic-g", by.node = TRUE)[tarNode]
    # if(newScore == -Inf)
    #   print("!blah")
    # else
    #  print("not blah!")

    allScores[i] = newScore

  }
  return (allScores)
}


#' This function calculates the LOPC or the LASSO parents of Score_LOPC and Score_LASSO
#' @keywords internal

CalculatePriors <- function(allNodes, tarNode, Allcombinations, score_mat, gamma)
{

  beta <- score_mat
  tarIndex = sub("V", "", tarNode)
  tarIndex = as.numeric(sub("_1", "", tarIndex))
  allScores <- array(data = NA, dim = length(Allcombinations))

  All_possible_parents = allNodes[!grepl("_1", allNodes)]
  all_targets <- allNodes[grepl("_1", allNodes)]

  rownames(beta) <- all_targets
  colnames(beta) <- All_possible_parents

  beta_target <- beta[tarNode,]
  prior_empty <- sum(log(beta_target+gamma))


  allScores[1] <- prior_empty
  for (i in 2:length(Allcombinations))
  {
    parents = Allcombinations[[i]]
    condition <- names(beta_target) %in% parents
    prior_present <- sum(log((1-beta_target[condition]+gamma)))
    prior_absent <- sum(log(beta_target[!condition]+gamma))

    prior <- sum(prior_absent, prior_present)
    allScores[i] = prior
  }
  return (allScores)

}

#' Internal function.
#' For full exhaustive search, sum up the BIC score of each possible parent
#' combination and the corresponsding Eij and ¬Eij scores for the edges.
#' @param all_target_weights the E and ¬E calculated for each possible parents combination
#' @param bic_scores BIC scores for each target nodes for each possible parent combination
#' @param score This refers to the type of score to be calculated. The default is "Score_LOPC".
#' Other possible values are Score_LASSO and BIC.
#' @param gamma value of the gamma constant added to the final score, to move the very small edge score values away from zero
#' @param all_parent_combinations List of all the combination of parent nodes, with maximum number of parents limited to 3
#' @keywords internal

CombineBICandPriors <- function(bic_scores, priors)
{

  bic_and_prior = priors + bic_scores
  return(bic_and_prior)
}


#' Applies linear interpolation for missing (NA) values in the data. This function uses the package "ImputeTS"
#' @param data data.frame with missing values
#' @import imputeTS
#' @export
#' @return data.frame with missing values imputed with linear interpolation
ImputeData <-function(data)
{
  imp <- imputeTS::na_interpolation(data, option = "linear")
  return (as.data.frame(imp))
}

#' This function calculates the number and percentage of true positive arcs.
#' @param  arc_set data.frame with two columns, "from" and "to", that represents the arcs of the learned DBN.
#' @param validated_set data.frame with two columns, "from" and "to" consisting of pairs of nodes present in the "true" network.
#' @return Number and percentage of TP arcs present in the learned network
#' @example
#'
#' yeast_data_imputed <- imputeTS::na_interpolation(yeast_data)
#' GS_Score_BIC <- GreedySearchParentsWithPriors(data = yeast_data_imputed,
#' score_mat = NULL, maxP = 1, gamma = 0.0,score = "BIC", pred_pos = pred_pos_35)
#' Score_BIC_arcs <- MakeArcSet(GS_Score_BIC, yeast_data)
#' results <- PercentageValidatedFromRecovered(Score_BIC_arcs,
#' yeasts_arcs_validated)
#' @export
PercentageValidatedFromRecovered <- function(arc_set, validated_set)
{
  common <- dplyr::inner_join(arc_set, validated_set)
  return(c("number" = nrow(common), "percentage" = nrow(common)/nrow(arc_set) * 100))

}

#' @keywords internal
ConvertToBN <-function(MyNet)
{
  adjMat = MyNet$AdjMatrix
  p = ncol(adjMat)
  newAdjMat = matrix(0, p*2,p*2)
  names = paste0("V", 1:p)
  names = c(names, paste0(names,"_1"))
  colnames(newAdjMat) = names
  rownames(newAdjMat) =  colnames(newAdjMat)
  newAdjMat[1:p, 1:p+p] = t(MyNet$AdjMatrix)
  new_bn = bnlearn::empty.graph(colnames(newAdjMat))
  amat(new_bn) = newAdjMat
  new_bn
}

#' This function is for internal use only.
#' @keywords internal

GetAMAT <- function(mat)
{
  genes = nrow(mat) / 2
  newMat = mat[1:genes, (genes+1):ncol(mat)]

}

#' This function is used to compute the results for a simulation study. It computes individual precision, recall, hamming distance and balanced accuracy
#' of each learned DBN structure and then calculates the mean and sd of all the individual results.
#' @param DBNList a list of learned DBN structures
#' @param datasets a list of datasets simulated using the \code{SimulateData()} function. It holds the true DBN network structure (Adjaceny matrix)
#' @return a list of different parameters:"precision_mean","precision_sd","recall_mean", "recall_sd",
#' "hamming_mean", "hamming_sd", "TN_mean", "TN_sd","FN_mean" , "FN_sd", "FP_mean", "FP_sd", "TP_mean",
#' "TP_sd", "Accuracy" , "Accuracy_sd".
#' @import G1DBN
#' @import caret
#' @export
CalculatePrecisionAndRecallForMultiple <- function(DBNList, datasets)
{
  dbns <- lapply(DBNList, function(x) GetAMAT(amat(x)))
  reals <- lapply(datasets, function(x) t(x$RealNet$AdjMatrix))
  realdbns <- lapply(datasets, function (x) ConvertToBN(x$RealNet))
  confs_mats <- mapply(function(dbn, real)
  {
    confusionMatrix( reference = factor(real), data = factor(dbn), positive = "1")

  }, dbns, reals, SIMPLIFY = FALSE)
  precisions = unlist(lapply(confs_mats, function(x) x$byClass[3]))
  recalls = unlist(lapply(confs_mats, function(x) x$byClass[1]))
  balanced_Acc <- unlist(lapply(confs_mats, function(x) x$byClass[11]))
  hammings = unlist(mapply(function(x,y) hamming(learned = x, true = y), DBNList, realdbns, SIMPLIFY = FALSE))
  TN <- unlist(lapply(confs_mats, function(x) x$table[1,1]))
  FN <- unlist(lapply(confs_mats, function(x) x$table[1,2]))
  FP <- unlist(lapply(confs_mats, function(x) x$table[2,1]))
  TP <- unlist(lapply(confs_mats, function(x) x$table[2,2]))

  precision_mean <- mean(precisions, na.rm = TRUE)
  precision_sd <- sd(precisions,  na.rm = TRUE)
  recall_mean <- mean(recalls, na.rm = TRUE)
  recall_sd <- sd(recalls, na.rm = TRUE)
  hamming_mean <- mean(hammings, na.rm = TRUE)
  hamming_sd <- sd(hammings, na.rm = TRUE)
  TN_mean <- mean(TN, na.rm = TRUE)
  TN_sd <- sd(TN, na.rm = TRUE)
  FN_mean <- mean(FN, na.rm = TRUE)
  FN_sd <- sd(FN, na.rm = TRUE)
  FP_mean <- mean(FP, na.rm = TRUE)
  FP_sd <- sd(FP, na.rm = TRUE)
  TP_mean <- mean(TP, na.rm = TRUE)
  TP_sd <- sd(TP, na.rm = TRUE)
  acc_mean <-  mean(balanced_Acc, na.rm = TRUE)
  acc_sd <- sd(balanced_Acc, na.rm = TRUE)

  return (list("precision_mean" = precision_mean, "precision_sd" = precision_sd,
               "recall_mean" = recall_mean, "recall_sd" = recall_sd,
               "hamming_mean" = hamming_mean, "hamming_sd" = hamming_sd,
               "TN_mean" = TN_mean, "TN_sd" = TN_sd,
               "FN_mean" = FN_mean, "FN_sd" = FN_sd,
               "FP_mean" = FP_mean, "FP_sd" = FP_sd,
               "TP_mean" = TP_mean, "TP_sd" = TP_sd,
               "Accuracy" = acc_mean, "Accuracy_sd" = acc_sd))


}
