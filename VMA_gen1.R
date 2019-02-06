library(R.matlab)
library(rstan)

setwd("/home/soyeon/")

# read the data file
rm(list=ls())
data = read.table("data.txt", header=T, sep="\t")
data = data[1:9500, ]

subjList <- unique(data[,"subjID"])
numSubjs <- length(subjList)
T = table(data$subjID)[1]
Tsubj <- as.vector(rep(0, numSubjs)) # number of trials for each subject
for (sIdx in 1:numSubjs)  {
  curSubj     <- subjList[sIdx]
  Tsubj[sIdx] <- sum(data$subjID == curSubj)  # Tsubj[N]
}

POI     <- c("mu_kappa", "mu_eta", "mu_theta", "mu_beta", "mu_bias", "log_lik",
             "kappa", "eta", "theta", "beta", "bias")

r <- array(0, c(numSubjs, T))
s <- array(0, c(numSubjs, T))
a <- array(0, c(numSubjs, T))
Avoid <- array(-1, c(numSubjs, T))

for (subjIdx in 1:numSubjs)   {
  #number of trials for each subj.
  useTrials                      <- Tsubj[subjIdx]
  currID                         <- subjList[subjIdx]
  data_curSubj                <- subset(data, data$subjID == currID)
  r[subjIdx, 1:useTrials] <- data_curSubj[, "r"]
  s[subjIdx, 1:useTrials] <- data_curSubj[, "s"]
  a[subjIdx, 1:useTrials] <- data_curSubj[, "a"]
  
  for (tIdx in 1:useTrials) {
    Y_t   <- data_curSubj[tIdx, "a"] # chosen Y on trial "t"
    Avoid[subjIdx , tIdx] <- Y_t
  }
}

dataList <- list(
  N      = numSubjs,
  T      = T,
  Tsubj  = Tsubj ,
  r      = r,
  s      = s,
  a      = a,
  Avoid  = Avoid
)

output2 = stan("VMA_gen1.stan", data = dataList, pars = POI,
              iter = 4000, warmup=2000, chains=4, cores=4)

save(output2, flie = "VMA_gen1.Rdata")

traceplot(output2)
plot(output2, pars = "bias", show_density = T)
print(output2)
output2_param <- rstan::extract(output2)

