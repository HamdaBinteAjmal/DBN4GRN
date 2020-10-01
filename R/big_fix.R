# Fixing bug #
GreedySearchParentsWithPriors_fix <- function(data,  score_mat = NULL, maxP = Inf, gamma = 0, pred_pos = NULL, score = "Score_BIC")
{
  p <- ncol(data)
  n <- nrow(data)
  dataS <- ShiftData(data)
  if(is.null(score_mat))
  {
   #break()
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

    GreedySearchWithPriors_OneTarget_fix(score_mat = score_mat, targetIdx = x, all_parents = all_parents, dataS = dataS, maxP = maxP,  gamma = gamma, score = score))

  fullBN <- do.call(what = rbind, args = local_bns_arcs)
  nodes <- names(dataS)
  bn <- bnlearn::empty.graph(nodes)
  bnlearn::arcs(bn) <- fullBN

  return (bn)


}
GreedySearchWithPriors_OneTarget_fix <- function(score_mat, targetIdx, all_parents, dataS, maxP = maxP, gamma, score)
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
#  print(paste0("old: ", old_bic))
  iterations <- 0
  while (keep_adding > 0 )#&& TotalParents < maxP)
  {
    parents <- bnlearn::parents(x = bn, node = targetNode)
    absents <- all_parents_names[which(!all_parents_names %in% parents)]
         new_bics <- unlist(lapply(all_parents_names, function(x)
    {


      if(x %in% absents)
      {
        if(TotalParents < maxP)
        {
         bn_1 <- bnlearn::set.arc(bn, from = x, to = targetNode)
        }
        else
        {
          bn_1 <- bn
        }

      }
      else
      {
        bn_1 <- bnlearn::drop.arc(bn, from = x, to = targetNode)

      }

      score <- bnlearn::score(x  = bn_1, data = data_subset, type = "bic-g", by.node = TRUE )[targetNode]
       current_parents <- bnlearn::parents(x = bn_1, node = targetNode)
      current_absents <- all_parents_names[which(!all_parents_names %in% current_parents)]
      weights_present <- sum(weight_vector_present[current_parents])
      weights_absent <- sum(weight_vector_absent[current_absents])
       if(score == -Inf)
      {
        print("ALERT -Inf ")
       # break()
      }
        sum(score, weights_present, weights_absent)

    }))


    bic_and_weight <- new_bics #+ weight_vector
    # print(sum(new_bics))
    difference <- bic_and_weight - old_bic
    #print(difference)
    diff <- 0
    keep_adding <- length(which(difference > diff)) > 0

    if(keep_adding)
    {

   #   print(difference)
      max_idx <- which.max(difference)
      if(all_parents[max_idx] %in% absents)
      {
        best_parent <- all_parents_names[max_idx]
        bn <- bnlearn::set.arc(bn, from = best_parent, to = targetNode)
        #TotalParents <- TotalParents + 1
        print(paste0("Adding Parent: ", best_parent ))
      }
      else
      {
        worst_parent <- all_parents_names[max_idx]
        bn <- bnlearn::drop.arc(bn, from = worst_parent, to = targetNode)
        #TotalParents <- TotalParents - 1
        print(paste0("Removing Parent: ", worst_parent ))
      }
      TotalParents <- nrow(bnlearn::arcs(bn))
      #all_parents[max_idx] <- NA
      old_bic <- bic_and_weight[max_idx]
    }


  }

  return (bnlearn::arcs(bn))

}


CreateBlackList <- function(data, targetNode)
{
  Slice1names = paste0("V", 1:ncol(data), "_1")
  Slice0names = paste0("V", 1:ncol(data))

  bl1 = expand.grid(Slice1names, Slice0names)
  bl2 = expand.grid(Slice0names, Slice0names)
  bl3 = expand.grid(Slice1names, Slice1names)

  bl4 = expand.grid(Slice0names, Slice1names)
 # colnames(bl4) = c("from", "to")
 # print(bl4)
 # print(Slice1names)
  nonTargets <- Slice1names[Slice1names!=targetNode]
# print(nonTargets)
  #bl4$to %in% nonTargets
# bl4$to
# bl4$to %in% nonTargets
  bl4 = bl4[bl4[,2] %in% nonTargets,]
 # print(bl4)
#  break()
  blacklist = rbind(bl1,bl2,bl3, bl4)
  blacklist <- as.matrix(blacklist)
  colnames(blacklist) = c("from", "to")
  blacklist
}


## Test greeddy search

for (j in 1:100)
{
  print(paste0("dataset: ", j))
for (i in 1:50)
{
  targetNode <- paste0("V",i,"_1")
#  print(targetNode)
  bl <- CreateBlackList(datasets[[j]]$data,targetNode )
  bnlearn_gs <- bnlearn::hc(x = ShiftData(datasets[[j]]$data), blacklist = bl, score = "bic-g" , maxp = 3, debug = FALSE)
  a = as.data.frame(arcs(bnlearn_gs))
  a = a[which(a$to == targetNode ),]
  #print(a)

  b = as.data.frame(arcs(BNs_BIC_greedy_3_fix[[j]]))
  b = b[which(b$to == targetNode ),]
 # print(b)
  if(!setequal(a$from, b$from))
  {
   #print(a)
   # print(b)
    print(paste0("ERROR - MISMATCH in dataset ",j, " node= ", i ))
   # break()
  }


}
}

for (j in 1:100)
{
  print(paste0("dataset: ", j))
for (i in 1:50)
{
  targetNode <- paste0("V",i,"_1")
 # print(targetNode)
  bl <- CreateBlackList(datasets[[j]]$data,targetNode )
  bnlearn_gs <- bnlearn::hc(x = ShiftData(datasets[[j]]$data), blacklist = bl, score = "bic-g" , maxp = Inf, debug = FALSE)
  a = as.data.frame(arcs(bnlearn_gs))
 # a = a[which(a$to == targetNode ),]
#  print(a)

  b = as.data.frame(arcs(BNs_BIC_greedy_fix[[j]]))
  b = b[which(b$to == targetNode ),]
 #  print(b)
  if(!setequal(a$from, b$from))
  {
  #  print(a)
 #   print(b)
    print(paste0("ERROR - MISMATCH in dataset ",j, " node= ", i ))
   # break()
  }


}
}
