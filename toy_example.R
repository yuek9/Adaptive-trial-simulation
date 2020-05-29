totaln=720
library(sandwich)
event.rates = c(0.3, 0.3)
stages=2
ntimes=500
pick=NULL
drift=2

see = function(null){
{
arms <- length(event.rates)
## nvec is the sample size allocation for each stage (the total sample size by the end of each stage)
nvec=round(seq(totaln/stages,totaln,length.out = stages))
p.assign <- rep(1/arms,arms)

## loop over each interim analysis, 
## adjusting randomization frequencies based on outcomes ######## this is the current rationale of adaption
## assume all outcomes are available each month

## store each stage events and sample size separately
events.stage = NULL
sample.stage = NULL
for (i in 1:stages) {
  
  ## simulate random binomial failures (function rbinom) according to event.rates,
  ##  use number of persons divided by number of stages for each analysis 
  
  stage_size = nvec[i]-ifelse(i==1, 0, nvec[i-1])
  n_each_arm = p.assign*stage_size 
  n_each_arm = round(n_each_arm)
  diff = stage_size - sum(n_each_arm)
  if(diff<0){n_each_arm = n_each_arm - rmultinom(n=1, size=-diff, prob=rep(1/4, arms))}
  if(diff>0){n_each_arm = n_each_arm + rmultinom(n=1, size=diff, prob=rep(1/4, arms))}  ######## this deal with rounding issue, by rounding and assigning the difference to random selected several arms
  
  if ((sum(n_each_arm) != stage_size)) stop('Allocation of arms mismatch the sample size')
  
  events.observed <- mapply(rbinom, size=n_each_arm, prob=event.rates, n=1)
  
  ## keep record of total number of allocated samples to each arm
  sample.stage = cbind(sample.stage, n_each_arm)
  events.stage = cbind(events.stage, events.observed)
  cum.events = rowSums(events.stage)
  cum.sample = rowSums(sample.stage)
  
  ## now update adaptive randomization probabilities based on observed events numbers
  ## can weight each arm by proportion of events occurring in that arm out of total events
  
  ####### naive adaption allocation: (but then zero allocaiton proportion will happen; may not want this due to randomly low events)
  # p.assign <- cum.events / cum.sample
  
  ####### another more stable allocation method: mentioned in Korn2011 paper "Outcome-Adaptive Randomization: Is It Useful?" 
  ####### the method by Thall and Wathen 
  # a = nvec[i]/(2*tail(nvec, 1))
  # p = cum.events / cum.sample
  # p.assign <- p^a/(p^a+(1-p)^a)  ############ still not able to deal with zero counts
  
  
  ####### another assigning method: based on Bayesian posterior mean of event rate, and assuming more event desired, assigning prob proportional to posterior mean
  ####### uniform prior used
  # a = nvec[i]/(2*tail(nvec, 1))
  # pos_mean = (cum.events^a+1)/(cum.sample^a+2)
  pos_mean = (cum.events+1)/(cum.sample+2)
  p.assign <- pos_mean/sum(pos_mean)
  
  
  ## the drift happens at the end of the stage (and do not further drift at the end of the trail)
  ## currently we assume multiplicative drifting effect
  if(i!=stages){
    event.rates = event.rates*drift
    event.rates[event.rates>1] <-1
  }
  
}  ## end loop over stages

cum.events = rowSums(events.stage)
cum.sample = rowSums(sample.stage)

## perform statistical testing of differences in event rates
## currently we only do statistical analysis once at the very end, do not drop arm/early stopping at interim analysis

return_arm_con = NULL
p_arm = NULL
events.stage.vec = NULL
sample.stage.vec = NULL
for(i in 1:stages){
  events.stage.vec = c(events.stage.vec, events.stage[,i])
  sample.stage.vec = c(sample.stage.vec, sample.stage[,i])
}

arm.short = as.factor(rep(1:arms, stages))
stage.short = as.factor(rep(1:stages, each=arms))


#transform variables into long formate (binary response)
binary.response.vec = do.call(c, lapply(1:length(events.stage.vec), function(i) c(rep(1, events.stage.vec[i]), rep(0, sample.stage.vec[i]-events.stage.vec[i]))))
arm.binomial = rep(1:arms, stages)
stage.binomial = rep(1:stages, each=arms)
arm.binary = as.factor(do.call(c, lapply(1:length(events.stage.vec), function(i)rep(arm.binomial[i], sample.stage.vec[i]))))
stage.binary = as.factor(do.call(c, lapply(1:length(events.stage.vec), function(i)rep(stage.binomial[i], sample.stage.vec[i]))))

}

tmp = glm(cbind(events.stage.vec, sample.stage.vec - events.stage.vec)~arm.short+stage.short, family=binomial)
#summary(tmp)
sqrt(diag(vcovHC(tmp, 'HC')))

data.model = data.frame(y = binary.response.vec,
                        arm = arm.binary, stage = stage.binary)

model_l = glm(y~arm*stage, family=quasibinomial, data=data.model)
#summary(model_l)
coefs_l = summary(model_l)$coefficient[1:arms,]
coef_cov = vcovHC(model_l, type='HC')[1:arms,1:arms] # sandwich estimator is used
coefs_l[,4] = pt(-abs(coefs_l[,1])/sqrt(diag(coef_cov)), df=model_l$df.residual)
coefs_l[,2] = sqrt(diag(coef_cov))[1:arms]
round(coefs_l[2,],4)


# tmp_index = lower.tri(matrix(0, arms, arms))
# index = cbind(matrix(rep(1:arms, arms), byrow=T)[tmp_index], matrix(rep(1:arms, each=arms), byrow=T)[tmp_index])
# p_value_con = rep(0, dim(index)[1])
# trt_ratio_est_con = rep(0, dim(index)[1])
# logit_p_arm = coefs[,1]+c(0, rep(coefs[1,1], arms-1))
# p_arm = exp(logit_p_arm)/(1+exp(logit_p_arm))



model = glm(y~arm*stage, family = quasipoisson, data=data.model)
# model = glm(y~arm, family = poisson, data=data.model)  
#summary(model)

coefs = summary(model)$coefficient[1:arms,]
coef_cov = vcovHC(model, type='HC')[1:arms, 1:arms]
coefs[,4] = pt(-abs(coefs[,1])/sqrt(diag(coef_cov)), df=model$df.residual)
coefs[,2] = sqrt(diag(coef_cov))[1:arms]
round(coefs[2,], 4)

# tmp_index = lower.tri(matrix(0, arms, arms))
# index = cbind(matrix(rep(1:arms, arms), byrow=T)[tmp_index], matrix(rep(1:arms, each=arms), byrow=T)[tmp_index])
# p_value_con_poi = rep(0, dim(index)[1])
# trt_ratio_est_con_poi = rep(0, dim(index)[1])
# log_p_arm = coefs[,1]+c(0, rep(coefs[1,1], arms-1))
# p_arm_poi = exp(log_p_arm)


# modify the p value by robust sd; notice the df=1, should use t distribution
return(c(coefs_l[2,4], coefs[2,4]))
}

res = replicate(500, see(1))
mean(res[1,]>0.05, na.rm=T)
mean(res[2,]>0.05)
