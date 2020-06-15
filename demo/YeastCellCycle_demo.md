knit: (function(input\_file, encoding) { out\_dir &lt;- 'docs'; rmarkdown::render(input\_file, encoding=encoding, output\_file=file.path(dirname(input\_file), out\_dir, 'index.html'))})

Learning sparse DBN structures from Yeast Cell Cycle gene expression data.
--------------------------------------------------------------------------

This document contains a demonstration to reverse engineer gene regulatory networks using DBN learning with BIC, *S**c**o**r**e*<sub>*L**O**P**C*</sub> and *S**c**o**r**e*<sub>*L**A**S**S**O*</sub> from Spellman's Yeast Cell Cycle data (Spellman et al. (1998) Mol Biol Cell. 1998 Dec; 9(12): 3273–3297. ( The name of our R package is 'DBN4GRN'. For further details of algorithms and methods, please see our paper "Dynamic Bayesian network learning to infer sparsemodels from time series data", which is currently submitted to the journal for review, or email us at r

The first step is to load the Spellmans Yeast Cell Cycle data, that is provided in the package.

``` r
data("yeast_data")
```

Now we impute the data using linear interpolation function from the "ImputeTS" package.

``` r
yeast_data_imputed <- ImputeData(yeast_data)
```

Now, using the data frame of validated arcs downloaded from YEASTRACT, we filter a list of all the possible regulator genes documented in YEASTRACT and their column indices in the data frame.

``` r
data("yeasts_arcs_validated")
potential_parents <- c("ACE2","ASH1","CHA4","FKH1","FKH2","GAT3","GCR1","HCM1","KAR4","MBP1", "MCM1","MET28","MIG2","NDD1","PHD1","PLM2","RAP1","RDS2","RME1","RPI1","SRD1","STB5","STP2","STP4","SUT1","SWI4","SWI5","SWI6","TBF1","TEC1","TOS4","WTM1","WTM2","YHP1","YOX1")

pred_pos_35 <- which(names(yeast_data) %in% potential_parents)
```

To learn a structure of a DBN using the BIC score, we use the following code.

``` r
GS_Score_BIC <- GreedySearchParentsWithPriors(data = yeast_data_imputed, score_mat = NULL, maxP = Inf, gamma = 0.0,score = "BIC", pred_pos = pred_pos_35)
#Get a list of pair of connected nodes in the network
Score_BIC_arcs <- MakeArcSet(GS_Score_BIC, yeast_data)
# To calculate the percentages of arcs learned in the DBN, we use the following function
results <- PercentageValidatedFromRecovered(Score_BIC_arcs, yeasts_arcs_validated)
xlsx::write.xlsx(Score_BIC_arcs, file = "GreedySearch_BIC.xlsx", col.names = TRUE, row.names = FALSE)
```

``` r
print(paste0("Percentage of validated arcs = ", results[2]))
```

    ## [1] "Percentage of validated arcs = 19.1770756796473"

``` r
print(paste0("Total number of arcs = ", nrow(Score_BIC_arcs) ))
```

    ## [1] "Total number of arcs = 8166"

To learn a DBN from yeast data using *S**c**o**r**e*<sub>*L**O**P**C*</sub> with *γ* = 0.0, we use the following code.

``` r
LOPC  = ExecuteG1DBNS1(yeast_data_imputed, predPosition = pred_pos_35)
GS_Score_LOPC_0p0 <- GreedySearchParentsWithPriors(data = yeast_data_imputed, score_mat = LOPC$S1ls, maxP = Inf, gamma = 0.0,score = "Score_LOPC", pred_pos = pred_pos_35)
#Get a list of pair of connected nodes in the network
Score_LOPC_0p0_arcs <- MakeArcSet(GS_Score_LOPC_0p0, yeast_data)
# To calculate the percentages of arcs learned in the DBN, we use the following function
results <- PercentageValidatedFromRecovered(Score_LOPC_0p0_arcs, yeasts_arcs_validated)
```

``` r
print(paste0("Percentage of validated arcs = ", results[2]))
```

    ## [1] "Percentage of validated arcs = 18.6193793540215"

``` r
print(paste0("Total number of arcs = ", nrow(Score_LOPC_0p0_arcs) ))
```

    ## [1] "Total number of arcs = 1579"

``` r
xlsx::write.xlsx(Score_BIC_arcs, file = "GreedySearch_LOPC.xlsx", col.names = TRUE, row.names = FALSE)
```

To learn a DBN from yeast data using *S**c**o**r**e*<sub>*L**O**P**C*</sub> with *γ* = 0.2, we use the following code.

``` r
#LOPC  = ExecuteG1DBNS1(yeast_data_imputed)
GS_Score_LOPC_0p2 <- GreedySearchParentsWithPriors(data = yeast_data_imputed, score_mat = LOPC$S1ls, maxP = Inf, gamma = 0.2,score = "Score_LOPC", pred_pos = pred_pos_35)
#Get a list of pair of connected nodes in the network
Score_LOPC_0p2_arcs <- MakeArcSet(GS_Score_LOPC_0p2, yeast_data)
# To calculate the percentages of arcs learned in the DBN, we use the following function
results <- PercentageValidatedFromRecovered(Score_LOPC_0p2_arcs, yeasts_arcs_validated)
xlsx::write.xlsx(Score_LOPC_0p2_arcs, file = "GreedySearch_LOPC_0p2.xlsx", col.names = TRUE, row.names = FALSE)
```

``` r
print(paste0("Percentage of validated arcs = ", results[2]))
```

    ## [1] "Percentage of validated arcs = 20.3909673070442"

``` r
print(paste0("Total number of arcs = ", nrow(Score_LOPC_0p2_arcs) ))
```

    ## [1] "Total number of arcs = 2967"

To learn a DBN from yeast data using *S**c**o**r**e*<sub>*L**A**S**S**O*</sub> with *γ* = 0.0, we use the following code.

``` r
LASSO  = ApplyLasso(yeast_data_imputed, pred_pos = pred_pos_35)
GS_Score_LASSO_0p0 <- GreedySearchParentsWithPriors(data = yeast_data_imputed, score_mat = LASSO, maxP = Inf, gamma = 0.0,score = "Score_LASSO", pred_pos = pred_pos_35)
Score_LASSO_0p0_arcs <- MakeArcSet(GS_Score_LASSO_0p0, yeast_data_imputed)
results <- PercentageValidatedFromRecovered(Score_LASSO_0p0_arcs, yeasts_arcs_validated)
xlsx::write.xlsx(Score_LASSO_0p0_arcs, file = "GreedySearch_LASSO.xlsx", col.names = TRUE, row.names = FALSE)
```

``` r
print(paste0("Percentage of validated arcs = ", results[2]))
```

    ## [1] "Percentage of validated arcs = 18.9097103918228"

``` r
print(paste0("Total number of arcs = ", nrow(Score_LASSO_0p0_arcs) ))
```

    ## [1] "Total number of arcs = 1174"

To learn a DBN from yeast data using *S**c**o**r**e*<sub>*L**A**S**S**O*</sub> with *γ* = 0.2, we use the following code.

``` r
#LASSO  = ApplyLasso(yeast_data_imputed, pred_pos = pred_pos_35)
GS_Score_LASSO_0p2 <- GreedySearchParentsWithPriors(data = yeast_data_imputed, score_mat = LASSO, maxP = Inf, gamma = 0.2,score = "Score_LASSO", pred_pos = pred_pos_35)
#Get a list of pair of connected nodes in the network
Score_LASSO_0p2_arcs <- MakeArcSet(GS_Score_LASSO_0p2, yeast_data_imputed)
# To calculate the percentages of arcs learned in the DBN, we use the following function
results <- PercentageValidatedFromRecovered(Score_LASSO_0p2_arcs, yeasts_arcs_validated)
xlsx::write.xlsx(Score_LASSO_0p2_arcs, file = "GreedySearch_LASSO_0p2.xlsx", col.names = TRUE, row.names = FALSE)
```

``` r
print(paste0("Percentage of validated arcs = ", results[2]))
```

    ## [1] "Percentage of validated arcs = 19.3619200936631"

``` r
print(paste0("Total number of arcs = ", nrow(Score_LASSO_0p2_arcs) ))
```

    ## [1] "Total number of arcs = 6833"
