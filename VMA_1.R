rm(list=ls())  

library(R.matlab)
library(rstan)

setwd("/Users/soyeonkim/Desktop/AMT_group/AMT_analysis/VBA_subfunctions")

# read the data file
rm(list=ls())
data = read.table("data.txt", header=T, sep="\t")

subjList <- unique(data[,"subjID"])
numSubjs <- length(subjList)
T = table(data$subjID)[1]
Tsubj <- as.vector(rep(0, numSubjs)) # number of trials for each subject
for (sIdx in 1:numSubjs)  {
  curSubj     <- subjList[sIdx]
  Tsubj[sIdx] <- sum(data$subjID == curSubj)  # Tsubj[N]
}

POI     <- c("mu_kappa", "mu_eta", "mu_beta", "mu_bias", "log_lik",
             "kappa", "eta", "beta", "bias")

r <- array(0, c(numSubjs, T))
a <- array(0, c(numSubjs, T))
s <- array(0, c(numSubjs, T))
Avoid <- array(-1, c(numSubjs, T))

for (subjIdx in 1:numSubjs)   {
  #number of trials for each subj.
  useTrials                      <- Tsubj[subjIdx]
  currID                         <- subjList[subjIdx]
  data_curSubj                <- subset(data, data$subjID == currID)
  r[subjIdx, 1:useTrials] <- data_curSubj[, "r"]
  a[subjIdx, 1:useTrials] <- data_curSubj[, "a"]
  s[subjIdx, 1:useTrials] <- data_curSubj[, "s"]
  
  for (tIdx in 1:useTrials) {
    Y_t   <- data_curSubj[tIdx, "a"] # chosen Y on trial "t"
    Avoid[subjIdx , tIdx] <- Y_t
  }
}

dataList <- list(
  N      = numSubjs,
  T      = T,
  Tsubj  = Tsubj,
  r      = r,
  a      = a,
  s      = s, 
  Avoid = Avoid
)

output = stan("VMA_1.stan", data = dataList, pars = POI,
              iter = 4000, warmup=2000, chains=4, cores=4)

save(output, file = "VMA_1.Rdata")

traceplot(output)
plot(output, pars = "bias", show_density = T)
print(output)
output_log_param <- rstan::extract(output_log)
output_param <- rstan::extract(output)
output_param$bias
