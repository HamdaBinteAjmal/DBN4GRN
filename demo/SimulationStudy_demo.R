
# This documents contains a demonstration to learn Dynamic Bayesian Networks from simulated timeseries data,
# using BIC score, Score_LOPC and Score_LASSO. The name of our R package is 'DBN4GRN'.
# For further details of algorithms and methods, please see our paper "Dynamic Bayesian network learning to infer sparse models from time series data",
# which is currently submitted to the \textit{Bioinformatics} journal for review, or email us at \href{h.ajmal1@nuigalway.ie} or \href{michael.madden@nuigalway.ie}.

# For this demo, we simulate 100 random datasets, with number of genes set to 50 and number of timepoints set to 20.
# You can change this number as needed.
# In our paper, we have shown the results on 100 simulated time series datasets. This is an
# expensive process. We executed the code on Kay, which is ICHEC's primary supercomputer
# and Ireland's national supercomputer for academic researchers.
# We used mc.cores = 40 in our apply functions using R package "parallel".
# If using on a multicore system, please change the mc.cores and apply function to mcapply and import
# package "parallel"

library(DBN4GRN)
## DATA GENERATION ##
sim_num <- 100 #Change this number to the desired number of simulations
seeds <- 1:sim_num
datasets <- lapply(seeds, function(x) SimulateData(genes = 50, timepoints = 20, seed = x, prop = 0.05))


## EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = BIC

BNs_BIC <- lapply(datasets, function(dataset) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_BIC"))
results = CalculatePrecisionAndRecallForMultiple(BNs_BIC, datasets)
print(results)


##  EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = Score_LOPC, gamma = 0.3
LOPCs <- lapply(datasets, function(dataset) ExecuteG1DBNS1(dataset$data))
BNs_Score_LOPC_0p3 <- mapply(function(dataset, lopc) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_LOPC", gamma = 0.3, score_mat = lopc$S1ls), datasets, LOPCs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LOPC_0p3, datasets)
print(results)


##  EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = Score_LOPC, gamma = 0.4
LOPCs <- lapply(datasets, function(dataset) ExecuteG1DBNS1(dataset$data))
BNs_Score_LOPC_0p4 <- mapply(function(dataset, lopc) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_LOPC", gamma = 0.4, score_mat = lopc$S1ls), datasets, LOPCs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LOPC_0p4, datasets)
print(results)


##  EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = Score_LOPC, gamma = 0.6
LOPCs <- lapply(datasets, function(dataset) ExecuteG1DBNS1(dataset$data))
BNs_Score_LOPC_0p6 <- mapply(function(dataset, lopc) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_LOPC", gamma = 0.6, score_mat = lopc$S1ls), datasets, LOPCs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LOPC_0p6, datasets)
print(results)

###  EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = Score_LASSO, gamma = 0.2


LASSOs <- lapply(datasets, function(dataset) ApplyLasso(dataset$data))
BNs_Score_LASSO <- mapply(function(dataset, lasso) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_LASSO", gamma = 0.2, score_mat = lasso), datasets, LASSOs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LASSO, datasets)
print(results)

###  EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = Score_LASSO, gamma = 0.3


#LASSOs <- lapply(datasets, function(dataset) ApplyLasso(dataset$data))
BNs_Score_LASSO_0p3 <- mapply(function(dataset, lasso) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_LASSO", gamma = 0.3, score_mat = lasso), datasets, LASSOs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LASSO_0p3, datasets)
print(results)

###  EXHAUSTIVE SEARCH with maximum number of parents set to 3. Score = Score_LASSO, gamma = 0.5


#LASSOs <- lapply(datasets, function(dataset) ApplyLasso(dataset$data))
BNs_Score_LASSO_0p5 <- mapply(function(dataset, lasso) ExhaustiveSearchForBestParents(data = dataset$data, type = "Score_LASSO", gamma = 0.5, score_mat = lasso), datasets, LASSOs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LASSO_0p5, datasets)
print(results)

### Greedy hill-climbing search, max number of parents set to Inf. Score = BIC
BNs_BIC_greedy <- lapply(datasets, function(dataset) GreedySearchParentsWithPriors(data = dataset$data, score = "BIC", maxP = Inf , gamma = 0.6))
results = CalculatePrecisionAndRecallForMultiple(BNs_BIC_greedy, datasets)
print(results)


### Greedy hill-climbing search, max number of parents set to 3. Score = BIC
BNs_BIC_greedy_3 <- lapply(datasets, function(dataset) GreedySearchParentsWithPriors(data = dataset$data, score = "BIC", maxP = 3 , gamma = 0.6))
results = CalculatePrecisionAndRecallForMultiple(BNs_BIC_greedy_3, datasets)
print(results)



### Greedy hill-climbing search, max number of parents set to 3. Score = Score_LOPC, gamma = 0.6
#LOPCs <- lapply(datasets, function(dataset) ExecuteG1DBNS1(dataset$data))
BNs_Score_LOPC_greedy_3 <- mapply(function(dataset, lopc) GreedySearchParentsWithPriors(data = dataset$data, score = "Score_LOPC", maxP = 3, gamma = 0.6, score_mat = lopc$S1ls ), datasets, LOPCs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LOPC_greedy_3, datasets)
print(results)



### Greedy hill-climbing search, max number of parents set to Inf. Score = Score_LOPC, gamma = 0.6
#LOPCs <- lapply(datasets, function(dataset) ExecuteG1DBNS1(dataset$data))

BNs_Score_LOPC_greedy <- mapply(function(dataset, lopc) GreedySearchParentsWithPriors(data = dataset$data, score = "Score_LOPC", maxP = Inf, gamma = 0.6, score_mat = lopc$S1ls ), datasets, LOPCs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LOPC_greedy, datasets)
print(results)



### Greedy hill-climbing search, max number of parents set to 3. Score = Score_LASSO, gamma = 0.6
BNs_Score_LASSO_greedy_3 <- mapply(function(dataset, lasso) GreedySearchParentsWithPriors(data = dataset$data, score = "Score_LASSO", maxP = 3, gamma = 0.6, score_mat = lasso ), datasets, LASSOs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LASSO_greedy_3, datasets)
print(results)





### Greedy hill-climbing search, no limit on max number of parents. Score = Score_LASSO, gamma = 0.6
BNs_Score_LASSO_greedy <- mapply(function(dataset, lasso) GreedySearchParentsWithPriors(data = dataset$data, score = "Score_Lasso", maxP = Inf, gamma = 0.6, score_mat = lasso ), datasets, LASSOs, SIMPLIFY = FALSE)
results = CalculatePrecisionAndRecallForMultiple(BNs_Score_LASSO_greedy, datasets)
print(results)

