# Code to test if Lebres old code and new code performs same.
G1s_old <- lapply(datasets[1:10], function(x) inferG1(x$data))

for( i in 1:10)
{
 k <- identical(G1s[[i]]$mat$S1ls, G1s_old[[i]]$ls)
  print(k)
}

G2s_new <- mapply(function(G1, dataset) G1DBN::DBNScoreStep2(S1 = G1$mat$S1ls, data = dataset$data, alpha1 = 0.7), G1s[1:10], datasets[1:10], SIMPLIFY = FALSE)
G2s_old <- mapply(function(g1, dataset) GfromG1(g1$ls, dataset$data, alpha1 = 0.7), G1s_old, datasets[1:10], SIMPLIFY = FALSE)
for(i in 1:10)
{
  G2s_old[[i]]$S2[which(G2s_old[[i]]$S2 == 1)] <- NA

}

  for( i in 1:10)
  {
    k <- identical(G2s_new[[i]], G2s_old[[i]]$S2)
    print(k)
  }
