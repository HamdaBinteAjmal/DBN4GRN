#' This function calculates the first order conditional independence score using the G1DBN package
#' @param data data.frame or a matrix with n rows (=time points) and p columns (=genes) containing the gene expression time series.
#' @param predPosition To be specified to reduce the set of possible predictor genes to a subset of d<p genes: an array included in [1,p] defining the position of the d predictor genes in the data matrix (n * p), default=NULL.
#' @return A matrix of first order conditional independence scores for each edge in the network
#' @examples
#' score_mat <- ExecuteG1DBNS1(data)
#' score_mat <- ExecuteG1DBNS1(yeast_data, pred_pos_35)
ExecuteG1DBNS1 <- function(data, predPosition = NULL)
{
  S1 <- G1DBN::DBNScoreStep1(data, method='ls', predPosition = predPosition)
  return(S1)
}

#' Simulate Data
#' This function is used to simulate the data from an AR(1) process. It makes a random network, with random values
#' in the adjaceny matrix and simulates data from that network.
#' @param genes Number of genes/variables
#' @param  timepoints Number of samples
#' @param seed Random seed,  is a number used to initialize a pseudorandom number generator.
#' @param prop The proportion of edges in the entire network
#' @param range vector with 4 elements specifying range values for the adjacency matrix generation (minimum negative value, maximum negative value, minimum positive value, maximum positive value)
#' @param errors vector of 2 elements specifying min and max value of the uniform distribution from which the error variances are drawn
#' @return A list comprising of dataset generated (data) and the network (RealNet) from which the data is generated.
#' @examples
#' SimulateData(genes = 50, timepoints = 20, seed = 1, prop = 0.05, range = c(-0.95,-0.05,0.05,0.95), errors = c(0.03, 0.08) )
#'
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
  Xn <- SimulGeneExpressionAR1(MyNet$A,B,X0,sigmaEps,n)

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



#' Applies L1-regularised Lasso regression on the dataset.
#' @param data a data.frame, columns represent the variables and rows represent the temporal data
#' @param pred_pos index of the potential predictor variables in the dataset. If the value is NULL, all variables are considered as potential predictors.
#' @return a matrix of Lasso regression coefficients.
#' @example LASSO <- ApplyLasso(yeast_data, c(1,2,3,5))

ApplyLasso <- function(data, pred_pos = NULL)
{
  if(is.null(pred_pos))
  {
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


#' Normalise the data values between 0-1
#' @param x a vector or a matrix of unscaled values
#' @return a vector or a matrix with values normalized within 0-1 range
#' @example lars_coeffs <- Normalisation_01(ApplyLasso(data))
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
#' @example
#' score_mat <- ExecuteG1DBNS1(data)
#' dbn <- GreedySearchParentsWithPriors(data, score_mat$S1ls, score = "Score_LOPC", gamma = 0.2)

GreedySearchParentsWithPriors <- function(data,  score_mat = NULL, maxP = Inf, gamma = 0, pred_pos = NULL, score = "Score_BIC")
{

  if(is.null(pred_pos))
  {
    all_parents <-  names(dataS)[1:p]
  }
  else
  {
    all_parents = names(dataS)[pred_pos]
  }
  p <- ncol(data)
  n <- nrow(data)
  dataS <- ShiftData(data)
  targets <- 1:p

  local_bns_arcs <- lapply(targets, function(x)

    GreedySearchWithPriors_OneTarget(score_mat = score_mat, targetIdx = x, all_parents = all_parents, dataS = dataS, maxP = maxP,  val = gamma, score = score))

  fullBN <- do.call(what = rbind, args = local_bns_arcs)
  nodes <- names(dataS)
  bn<- bnlearn::empty.graph(nodes)
  arcs(bn) <- fullBN

  return (bn)


}
#' Internal function
#' Applied greedy hill-climbing search to each local target node in the netork
#' E_ij = log(score_mat + gamma)
#' ¬E_ij =  log(1-score_mat+gamma
#' @return set of parents arcs for each local target node
#' @keywords internal
GreedySearchWithPriors_OneTarget<- function(score_mat, targetIdx, all_parents, dataS, maxP = maxP, gamma, score)
{
  keep_adding = 1
  targetNode <- paste0("V", targetIdx , "_1")
  all_parents_names <- paste0("V", all_parents)

  nodes <- c(all_parents_names, targetNode)

  print(paste0("Currently finding parents for: ", targetNode))
  data_subset <- dataS[,nodes]

  if(score == "Score_LASSO")
  {
    weights_absent <- log(1- score_mat + gamma)
    weights_present <- log(score_mat+gamma)
  }

  else if(score == "Score_LOPC")
  {
    weights_absent <- log(score_mat + gamma)
    weights_present <- log(1 - score_mat+gamma)
  }
  else
  {
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

    new_bics <- unlist(lapply(all_parents_names, function(x)
    {

      if(!is.na(x))
      {

        bn_1 <- bnlearn::set.arc(bn, from = x, to = targetNode)
        score <- bnlearn::score(x  = bn_1, data = data_subset, type = scoreType, by.node = TRUE )[targetNode]

        weights_present <- weight_vector_present[x]
        weights_absent <- sum(weight_vector_absent[which(names(weight_vector_absent) != x)])
      }
      else
      {
        score <- NA
      }

      sum(score, weights_present, weights_absent)
    }))
    #weight_vector <- sum(weights_present, weights_absent)
    bic_and_weight <- new_bics #+ weight_vector
    difference <- bic_and_weight - old_bic
    keep_adding <- length(which(difference > diff)) > 0

    if(keep_adding)
    {
      print(difference)
      max_idx <- which.max(difference)
      best_parent <- all_parents_names[max_idx]
      bn <- bnlearn::set.arc(bn, from = best_parent, to = targetNode)
      TotalParents <- TotalParents + 1
      print(paste0("Adding Parent: ", best_parent ))
      all_parents[max_idx] <- NA
      old_bic <- bic_and_weight[max_idx]
    }

  }

  return (arcs(bn))

}

#' Take a learned bnlearn Bayesian network structure and
#' returns a data frame in which each row represents a pair of nodes conected in the network. This function is used to retrieve the arcs from the network.
#' @param bn Bnlearn object Bayesian network
#' @param data data.frame from which the Bayesian network structure is learned.
#' @return arc_set a data.frame with two columns. First column shows the regulators, second column shows the targets.
#' @keywords internal
#'
#' @example
#' score_mat <- ExecuteG1DBNS1(data)
#' dbn <- GreedySearchParentsWithPriors(data, score_mat$S1ls, score = "Score_LOPC", gamma = 0.2)
#' arcs <- MarkArcSet(dbn)

MakeArcSet <- function(bn, data)
{
  arc_set <- arcs(bn)
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
#' @param data data.frame where columns represent the variables and rows contain the temporal data
ExhaustiveSearchForBestParents <- function(data, type = "Score_BIC", gamma = NULL, score_mat = NULL)
{
  p <- ncol(data)
  n <- nrow(data)
  dataS <- ShiftData(data)
  allNodes <- names(dataS)
  All_possible_parents <-  names(dataS)[1:p]
  Allcombinations = MakeAllPossibleCombination(All_possible_parents)
  All_possible_targets = allNodes[grepl("_1", allNodes)]

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
  BICandWeight <- CombineBICandPriors(bics, priors)
  BN <- ConstructBNUsingMaxScoringParents(full_score, Allcombinations, allNodes)
  return(BN)
}

#' This function takes in all scores for all parents combinations (columns) of all target nodes (rows)
#' and returns the maximum scoring network.
#' @param BICandWeigh Score_LASSO or Score_LOPC for each combination of parents for each target node
#' @param all_parent_combinations All possible combinations of parents , upto 3 parents max
#' @param allNodes all the target nodes in the network
#' @keywords internal
ConstructBNUsingMaxScoringParents <- function(BICandWeight,  all_parent_combinations, allNodes)
{
  max_ids <- apply(BICandWeight, 1, which.max)
  #allNodes = names(dataS)
  All_possible_targets = allNodes[grepl("_1", allNodes)]
  max_parents <- lapply(max_ids, function(x) unlist(all_parent_combinations[x]))
  bn <- empty.graph(allNodes)
  for (i in 1:length(All_possible_targets))

  {
    print(i)
    target <- All_possible_targets[[i]]
    parents <- max_parents[[i]]
    print(parents)
    if(!is.na(parents))
    {
      print(length(parents))

      for (par in parents)
      {
        bn <- set.arc(bn, from = par, to = target)
      }
    }
  }
  return (bn)
}


#' This methods calculates BIC score for each target node, for all the possible combinations of parents
#' @keywords internal
ScoreEntirepossibleSetsOfParents <- function(dataS, tarNode, all_combinations)
{

  tarIndex = sub("V", "", tarNode)
  tarIndex = as.numeric(sub("_1", "", tarIndex))
  allScores <- array(data = NA, dim = length(all_combinations))
  #  beta = beta[which(beta$to == tarNode),]
  allNodes = names(dataS)
  All_possible_parents = allNodes[!grepl("_1", allNodes)]

  bn = empty.graph(nodes = c(All_possible_parents, tarNode))
  ScoreOfGraph = bnlearn::score(bn, dataS[,c(nodes(bn))], type = "bic-g", by.node = TRUE)[tarNode]

  allScores[1] = ScoreOfGraph

  for (i in 2:length(all_parents))
  {
    bn = bnlearn::empty.graph(nodes = c(All_possible_parents, tarNode))

    parents = all_combinations[[i]]
    for (j in 1:length(parents))
    {
      if(!is.na(parents[j]))
      {
        bn <- set.arc(bn, from = parents[j], to = tarNode)
      }

    }

    newScore = bnlearn::score(bn, dataS[,c(nodes(bn))], type = "bic-g", by.node = TRUE)[tarNode]
    if(newScore == -Inf)
      print("!blah")

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
    prior_present <- sum(log((1-beta_target[condition]+val)))
    prior_absent <- sum(log(beta_target[!condition]+val))

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

CombineBICandWeight <- function(all_targets_weights, bic_scores, gamma = 0, score = "Score_LOPC", all_parent_combinations)
{
  lengths <- unlist(lapply(all_parent_combinations,length))
  adjusted <- lapply(all_targets_weights, function(x)
  {
    mapply(function(wt,lts)
    {

      multiplier * (wt + (adjustment * lts))
    }, x,lengths, SIMPLIFY = FALSE)
  })
  adjusted <- lapply(adjusted, unlist)


  BICandWeight <- t(mapply(function(x,y) {
    x +y
  }
  , bic_scores, adjusted, SIMPLIFY = TRUE))
  BICandWeight
}

#' Applies linear interpolation for missing (NA) values in the data
#' @param data data.frame with missing values
#' @return data.frame with missing values imputed with linear interpolation
ImputeData <-function(data)
{
  # imp <- imputeTS::na_interpolation(data, option = "linear")
  imp <- imputeTS::na.locf(data)
  return (as.data.frame(imp))
}

#' Calculates the number and percentage of true positive arcs
#' @param  arc_set data.frame with two columns, "from" and "to"
#' @param validated_set data.frame with two columns, "from" and "to" consisting of pairs of nodes present in the "true" network.
#' @return number and percentage of TP arcs present in the learned network
#' @example
#' GS_Score_BIC <- GreedySearchParentsWithPriors(data = yeast_data, score_mat = NULL, maxP = Inf, gamma = 0.0,score = "Score_BIC", pred_pos = pred_pos_35)
#' Score_BIC_arcs <- MakeArcSet(GS_Score_BIC, yeast_data)
#' results <- PercentageValidatedFromRecovered(Score_BIC_arcs, yeast_data)
PercentageValidatedFromRecovered <- function(arc_set, validated_set)
{
  common <- inner_join(arc_set, validated_set)
  return(c("number" = nrow(common), "percentage" = nrow(common)/nrow(arc_set) * 100))

}
