or = read.csv("AMT_summary_data.csv", header=T)
library('dplyr')
or <- or %>% select(sigmaN, sigmaA, eta, kappa, beta, bias)
or <- or[1:20, ]
kappa_non_mean = apply(output3_non_param$kappa, 2, mean)
eta_non_mean = apply(output3_non_param$eta, 2, mean)
sigmaN_non_mean = apply(output3_non_param$sigmaN, 2, mean)
sigmaA_non_mean = apply(output3_non_param$sigmaA, 2, mean)
beta_non_mean = apply(output3_non_param$beta, 2, mean)
bias_non_mean = apply(output3_non_param$bias, 2, mean)

kappa_mean = apply(output3_param$kappa, 2, mean)
eta_mean = apply(output3_param$eta, 2, mean)
sigmaN_mean = apply(output3_param$sigmaN, 2, mean)
sigmaA_mean = apply(output3_param$sigmaA, 2, mean)
beta_mean = apply(output3_param$beta, 2, mean)
bias_mean = apply(output3_param$bias, 2, mean)

comp_non <- data.frame(kappa_mean, kappa_non_mean,
                        eta_mean, eta_non_mean,
                        sigmaN_mean, sigmaN_non_mean,
                        sigmaA_mean, sigmaA_non_mean,
                        beta_mean, beta_non_mean,
                        bias_mean, bias_non_mean)

library(ggplot2)
ggplot(data = comp, aes(x= kappa_non_mean, y= kappa_mean)) +
  geom_point()+
  xlab("individual kappa")+
  ylab("hierarchical kappa")+
  xlim(0, 1)+
  ylim(0, 1)+
  geom_smooth(method=lm)

ggplot(data = comp, aes(x= eta_non_mean, y= eta_mean)) +
  geom_point()+
  xlab("individual eta")+
  ylab("hierarchical eta")+
  geom_smooth(method=lm)

ggplot(data = comp, aes(x= sigmaN_non_mean, y= sigmaN_mean)) +
  geom_point()+
  xlab("individual sigmaN")+
  ylab("hierarchical sigmaN")+
  ylim(0.025, 1)+
  geom_smooth(method=lm)

ggplot(data = comp, aes(x= sigmaA_non_mean, y= sigmaA_mean)) +
  geom_point()+
  xlab("individual sigmaA")+
  ylab("hierarchical sigmaA")+
  geom_smooth(method=lm)

ggplot(data = comp, aes(x= beta_non_mean, y= beta_mean)) +
  geom_point()+
  xlab("individual beta")+
  ylab("hierarchical beta")+
  geom_smooth(method=lm)

ggplot(data = comp, aes(x= bias_non_mean, y= bias_mean)) +
  geom_point()+
  xlab("individual bias")+
  ylab("hierarchical bias")+
  geom_smooth(method=lm)
