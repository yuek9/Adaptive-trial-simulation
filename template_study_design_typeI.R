
## author = Amalia Magaret, Kun Yue
## date started = 24 Jan 2019
## objective = provide skeleton for consideration of
#
## run this function using 'template_run_study.R' where
##   the parameters can be selected and the output manipulated

############# inputs #####################
## nvec = number of persons enrolled, a vector including all sizes of interest
## event.rates = vector of events rates within the arms
## stages = number of times you will perform interim analyses and
##   update randomization frequencies
## ntimes = number of times to run simulation

## for drift part, apply to event rate and at each end of stage
## if we look at the treatment effect in relative sense (look at the ratio, with drift happen), then consistent effect require multiplicative drift change
## so the bias will be comparede to true rate ratio (does not drift)

#####some check
# 1. if we do analysis conditional on stage: should be fine with power but lose power: question: how much power lost
# (information: conditional on stage, include this as a random effect, (logistic regerssion model/) or better with poisson model to get risk ratio and drift is multiplicative to the baseline rate and does not affect risk ratio)
# 2. if we do not conditional on stage, then able to update allocation ratio continuously: question: how high drift to be to induce large bias


###### updates of the codes:
# 03/19: changed re-allocation adaptive methods from  pos_mean = (cum.events^a+1)/(cum.sample^a+2) to pos_mean = (cum.events+1)/(cum.sample+2):
#        the re-allocation will be more extreme. a=cumulative_n/total_n, down-weighting the early results

# for nocondition version, used regression models not conditional on stage (only fixed effect of arm)
# for typeI version, under null hupothesis should pick all arms

library(lme4)
library(rlist)
library(sandwich)

totaln=720

event.rates = c(0.2, 0.2)
stages=2
ntimes=500
pick=NULL
drift=4

get.effect.power <- 
  function(totaln = 2000, # total sample size to enroll for all arms and stages together
           event.rates= seq(.05,.2,.05), # event rates for each arm
           stages= 4 ,
           ntimes= 500 , #repetition times for simulation of power/estimate
           pick = NULL, # specify the arm you want o pick for power analysis; pick=c(3, 4)
           drift = 1.3 # effect drift, multiplicative
           ) {


## arms is the number of tx arms to randomize to  
arms <- length(event.rates)
## nvec is the sample size allocation for each stage (the total sample size by the end of each stage)
nvec=round(seq(totaln/stages,totaln,length.out = stages))

## loop over number of simulations
one_simulation = function(seed){
  # set.seed(seed)
  
  ## initially randomize evenly
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
  
  ########### version one: do not conditional on stage; treatment effect estimated by sample mean; 
  if(T){
    return_arm = NULL
    event.rates_est = NULL
    message_raw = tryCatch({
      event.rates_est = cum.events/cum.sample
      event.rates_est
      
      
      # to pick out the one with the highest effect: might be able to do a two step testing
      # 1. do chi-square test on proportions 
      # 2. do post hoc 
      
      tmp_index = lower.tri(matrix(0, arms, arms))
      index = cbind(matrix(rep(1:arms, arms), byrow=T)[tmp_index], matrix(rep(1:arms, each=arms), byrow=T)[tmp_index])
      p_value = rep(0, dim(index)[1])
      trt_ratio_est = rep(0, dim(index)[1])
      trt_ratio_true = rep(0, dim(index)[1])
      options(warn=-1)
      for (k in 1:dim(index)[1]){
        invisible(capture.output(test_res <- prop.test(x=cum.events[index[k,]], n=cum.sample[index[k,]])))
        p_value[k] = test_res$p.value
        tmp = (cum.events/cum.sample)[index[k,]] #
        trt_ratio_est[k] =  tmp[1]/tmp[2]
        tmp = event.rates[index[k,]]
        trt_ratio_true[k] = tmp[1]/tmp[2]
      }
      options(warn=1)
      ## adjust p values to control multiple comparison familywise error
      p_value = p.adjust(p_value, method = 'holm')
      # print(round(p_value, 4))
      
      
      ## pick out the 'best' arm if there is significant results
      if(sum(p_value<0.05)==0){ # none is significant, nothing picked
        return_arm = 0
      } else { # at least one is significant
        
        # pick out significant pairs, and order to put the higher one in the first column
        sig_pair = index[p_value<0.05,, drop=F] 
        trt_ratio_sig_pair = trt_ratio_est[p_value<0.05]
        for(i in 1:dim(sig_pair)[1]) sig_pair[i,]<- sig_pair[i,c(-1/2, 1/2)*sign(trt_ratio_sig_pair[i]-1)+3/2]
        
        # pick out the highest mean group, and check significancy with the rest; start with the second highest
        ranking = order(event.rates_est, decreasing = T)
        arm_index= 2
        a1 = ranking[1]; a2 = a1; 
        return_arm = NULL
        while(sum((sig_pair[,1]==a1) * (sig_pair[,2]==a2))==0 && arm_index<=arms){ # if the current comparison pair not significant, the highest is same as this arm
          return_arm = c(return_arm, a2)
          a2 = ranking[arm_index]
          arm_index = arm_index+1
        }
      }
      
    }, error = function(e) e)
}
  ########## version 2: conditional on stage, logistic model
  if(T){
    return_arm_con = NULL
    p_arm = NULL
    events.stage.vec = NULL
    sample.stage.vec = NULL
    for(i in 1:stages){
      events.stage.vec = c(events.stage.vec, events.stage[,i])
      sample.stage.vec = c(sample.stage.vec, sample.stage[,i])
    }
    
    #transform variables into long formate (binary response)
    binary.response.vec = do.call(c, lapply(1:length(events.stage.vec), function(i) c(rep(1, events.stage.vec[i]), rep(0, sample.stage.vec[i]-events.stage.vec[i]))))
    arm.binomial = rep(1:arms, stages)
    stage.binomial = rep(1:stages, each=arms)
    arm.binary = as.factor(do.call(c, lapply(1:length(events.stage.vec), function(i)rep(arm.binomial[i], sample.stage.vec[i]))))
    stage.binary = as.factor(do.call(c, lapply(1:length(events.stage.vec), function(i)rep(stage.binomial[i], sample.stage.vec[i]))))
    
    data.model = data.frame(y = binary.response.vec,
                            arm = arm.binary, stage = stage.binary)
    
    
    message_logistic = tryCatch(error = function(e)e, warning = function(w) w,
                                {
                                  # model = glmer(cbind(events, sample-events)~arm+(1|stage), family = binomial, data=data.model)
                                  
                                  # attention: it seems robust sd is very different from model based sd, we try to include stage as dummy variable to be adjusted
                                  model = glm(y~arm+stage, family=quasibinomial, data=data.model)
                                  
         
        # fixed effect only, no condition on stage
        # model = glm(cbind(events, sample-events)~arm, family = binomial, data=data.model)
        
  
            # plot(as.numeric(data.model$arm), summary(model)$residual)
            # plot(as.numeric(data.model$stage), summary(model)$residual)
            
            
            #testing based on coefficients of arm, taken into consideration of covariance matrix, normal approximation
            ## notice: testing is based on logit(p) scale, which is equivalent to test the equality of each arm
            #treatment effect (ratio) estimation based on estimated (\hat p) by transformation from logit(p)
            
            coefs = summary(model)$coefficient[1:arms,]
            coef_cov = vcovHC(model, type='HC')[1:arms,1:arms] # sandwich estimator is used
            
            tmp_index = lower.tri(matrix(0, arms, arms))
            index = cbind(matrix(rep(1:arms, arms), byrow=T)[tmp_index], matrix(rep(1:arms, each=arms), byrow=T)[tmp_index])
            p_value_con = rep(0, dim(index)[1])
            trt_ratio_est_con = rep(0, dim(index)[1])
            logit_p_arm = coefs[,1]+c(0, rep(coefs[1,1], arms-1))
            p_arm = exp(logit_p_arm)/(1+exp(logit_p_arm))
            
            # modify the p value by robust sd; notice the df=1, should use t distribution
            coefs[,4] = pt(-abs(coefs[,1])/sqrt(diag(coef_cov)), df=model$df.residual)*2
            
            
            for (k in 1:dim(index)[1]){
              # reference arm is arm 1(always)
              if(1 %in% index[k,]){ #if one arm is arm 1, then the difference is just the coefficient
                other_arm = ifelse(index[k,1]==1, index[k,2], index[k,1])
                p_value_con[k] = as.numeric(coefs[other_arm, 4])
              }else{ # if both not arm 1, then difference is the contrast of coefficients
                t_value = as.numeric(c(1, -1)%*%coefs[index[k,],1]/sqrt(c(1, 1)%*%coef_cov[index[k,], index[k,]]%*%c(1, 1)))
                p_value_con[k] = (1-pt(abs(t_value), df=model$df.residual))*2
              }
              
              trt_ratio_est_con[k] =  p_arm[index[k,1]]/p_arm[index[k,2]]
            }
            
            ## adjust p values to control multiple comparison familywise error
            p_value_con = p.adjust(p_value_con, method = 'holm')
            
            
            return_arm_con = NULL
            ## pick out the 'best' arm if there is significant results
            if(sum(p_value_con<0.05)==0){ # none is significant, nothing picked
              return_arm_con = 0
            } else { # at least one is significant
              
              # pick out significant pairs, and order to put the higher one in the first column
              sig_pair_con = index[p_value_con<0.05,, drop=F] 
              trt_ratio_sig_pair_con = trt_ratio_est_con[p_value_con<0.05]
              
              for(i in 1:dim(sig_pair_con)[1]) sig_pair_con[i,]<- sig_pair_con[i,c(-1/2, 1/2)*sign(trt_ratio_sig_pair_con[i]-1)+3/2]
              
              # pick out the highest mean group, and check significancy with the rest; start with the second highest
              ranking = order(p_arm, decreasing = T)
              arm_index= 2
              a1 = ranking[1]; a2 = a1; 
              return_arm_con = NULL
              while(sum((sig_pair_con[,1]==a1) * (sig_pair_con[,2]==a2))==0 && arm_index<=arms){ # if the current comparison pair not significant, the highest is same as this arm
                return_arm_con = c(return_arm_con, a2)
                a2 = ranking[arm_index]
                arm_index = arm_index+1
              }
            }
            
  })
  
  
}
  ################## version 3: use poisson model (quasipoisson can not be used with mixed effect glmer???)
  if(T){
    return_arm_con_poi = NULL
    p_arm_poi = NULL
  events.stage.vec = NULL
  sample.stage.vec = NULL
  for(i in 1:stages){
    events.stage.vec = c(events.stage.vec, events.stage[,i])
    sample.stage.vec = c(sample.stage.vec, sample.stage[,i])
  }
  
  #transform variables into long formate (binary response)
  binary.response.vec = do.call(c, lapply(1:length(events.stage.vec), function(i) c(rep(1, events.stage.vec[i]), rep(0, sample.stage.vec[i]-events.stage.vec[i]))))
  arm.binomial = rep(1:arms, stages)
  stage.binomial = rep(1:stages, each=arms)
  arm.binary = as.factor(do.call(c, lapply(1:length(events.stage.vec), function(i)rep(arm.binomial[i], sample.stage.vec[i]))))
  stage.binary = as.factor(do.call(c, lapply(1:length(events.stage.vec), function(i)rep(stage.binomial[i], sample.stage.vec[i]))))
  
  data.model = data.frame(y = binary.response.vec,
                          arm = arm.binary, stage = stage.binary)
  
  message_poisson = tryCatch(error = function(e)e, warning = function(w) w, 
      {
             
        # model = glmer(y~(1|stage)+arm, family = quasipoisson, data=data.model)
        model = glm(y~arm+stage, family = quasipoisson, data=data.model)
        # model = glm(y~arm, family = poisson, data=data.model)  
           
  
        see = summary(model)
            # plot(as.numeric(arm.binary), see$residuals)
            # plot(as.numeric(stage.binary), see$residuals)
            
            #testing based on coefficients of arm, taken into consideration of covariance matrix, normal approximation
            ## notice: testing is based on log(p) scale, which is equivalent to test the equality of each arm
            #treatment effect (ratio) estimation based on estimated (\hat p) by transformation from log(p)
            
            coefs = summary(model)$coefficient[1:arms,]
            coef_cov = vcovHC(model, type='HC')[1:arms, 1:arms]
            
            tmp_index = lower.tri(matrix(0, arms, arms))
            index = cbind(matrix(rep(1:arms, arms), byrow=T)[tmp_index], matrix(rep(1:arms, each=arms), byrow=T)[tmp_index])
            p_value_con_poi = rep(0, dim(index)[1])
            trt_ratio_est_con_poi = rep(0, dim(index)[1])
            log_p_arm = coefs[,1]+c(0, rep(coefs[1,1], arms-1))
            p_arm_poi = exp(log_p_arm)
            
        
            # modify the p value by robust sd; notice the df=1, should use t distribution
            coefs[,4] = pt(-abs(coefs[,1])/sqrt(diag(coef_cov)), df=model$df.residual)*2
            
            
            for (k in 1:dim(index)[1]){
              # reference arm is arm 1(always)
              if(1 %in% index[k,]){ #if one arm is arm 1, then the difference is just the coefficient
                other_arm = ifelse(index[k,1]==1, index[k,2], index[k,1])
                p_value_con_poi[k] = as.numeric(coefs[other_arm, 4])
              }else{ # if both not arm 1, then difference is the contrast of coefficients
                t_value = as.numeric(c(1, -1)%*%coefs[index[k,],1]/sqrt(c(1, 1)%*%coef_cov[index[k,], index[k,]]%*%c(1, 1)))
                p_value_con_poi[k] = (1-pt(abs(t_value), df=model$df.residual))*2
              }
              
              trt_ratio_est_con_poi[k] =  p_arm_poi[index[k,1]]/p_arm_poi[index[k,2]]
            }
            
            
            ## adjust p values to control multiple comparison familywise error
            p_value_con_poi = p.adjust(p_value_con_poi, method = 'holm')

            return_arm_con_poi = NULL
            ## pick out the 'best' arm if there is significant results
            if(sum(p_value_con_poi<0.05)==0){ # none is significant, nothing picked
              return_arm_con_poi = 0
            } else { # at least one is significant
              
              # pick out significant pairs, and order to put the higher one in the first column
              sig_pair_con_poi = index[p_value_con_poi<0.05,, drop=F] 
              trt_ratio_sig_pair_con_poi = trt_ratio_est_con_poi[p_value_con_poi<0.05]
              
              for(i in 1:dim(sig_pair_con_poi)[1]) sig_pair_con_poi[i,]<- sig_pair_con_poi[i,c(-1/2, 1/2)*sign(trt_ratio_sig_pair_con_poi[i]-1)+3/2]
              
              # pick out the highest mean group, and check significancy with the rest; start with the second highest
              ranking = order(p_arm_poi, decreasing = T)
              arm_index= 2
              a1 = ranking[1]; a2 = a1; 
              return_arm_con_poi = NULL
              while(sum((sig_pair_con_poi[,1]==a1) * (sig_pair_con_poi[,2]==a2))==0 && arm_index<=arms){ # if the current comparison pair not significant, the highest is same as this arm
                return_arm_con_poi = c(return_arm_con_poi, a2)
                a2 = ranking[arm_index]
                arm_index = arm_index+1
              }
            }
           })
  
  }
  
  
  
  return(list('arm_picked' = return_arm, # the picked highest mean groups, based on multiple comparison corrected p values
              'all_treatment_effect' = event.rates_est, #based on empirical event ratio
              'arm_picked_conditional_logistic' = return_arm_con,
              'all_treatment_effect_conditional_logistic' = p_arm, # based on conditional regression removing stage effect
              'arm_picked_conditional_poisson' = return_arm_con_poi,
              'all_treatment_effect_conditional_poisson' = p_arm_poi, # based on conditional regression removing stage effect
              'message_raw' = message_raw,
              'message_logistic' = message_logistic,  #record the message and check for nonconvergence
              'message_poisson' = message_poisson)) 
} 
# bias of results are focused on ratio od effects
  
results = lapply(1:ntimes, one_simulation)

bad_estimate_logistic = 0
bad_estimate_poisson = 0
error_raw = 0
results_raw = list(NULL)
results_logistic = list(NULL)
results_poisson = list(NULL)

### it seems to have two types of bad estimation: Model failed to converge and Model is nearly unidentifiable
for(i in 1:ntimes){ #check if output contains warning messages; only keep analysis with that method converge (can have different number of analysis for each method)
  # it also happened that some unexplanined error occurred for raw analysis, remove those as well
  
  # the attributes of warning message is not null
  if(is.null(attributes(results[[i]]$message_logistic)) & is.null(attributes(results[[i]]$message_poisson))
     & is.null(attributes(results[[i]]$message_raw))){
    results_raw = list.append(results_raw, results[[i]])
    results_logistic = list.append(results_logistic, results[[i]])
    results_poisson = list.append(results_poisson, results[[i]])
  }else{
    cat('warning/error occurred:', c(i, event.rates, totaln, drift), '\n')
    error_raw = error_raw+1
    bad_estimate_logistic =  bad_estimate_logistic+1
    bad_estimate_poisson = bad_estimate_poisson+1
    
    if(is.null(attributes(results[[i]]$message_raw))){ #logistic no problem, keep it
      error_raw =  error_raw-1
      results_raw = list.append(results_raw, results[[i]])
      
    }
    
    
    if(is.null(attributes(results[[i]]$message_logistic))){ #logistic no problem, keep it
      bad_estimate_logistic =  bad_estimate_logistic-1
      results_logistic = list.append(results_logistic, results[[i]])
      
    }
    
    if(is.null(attributes(results[[i]]$message_poisson))){ #poisson no problem, keep it
      bad_estimate_poisson = bad_estimate_poisson - 1
      results_poisson = list.append(results_poisson, results[[i]])
    }
    
  }
} 
results_raw[[1]]<-NULL
results_logistic[[1]]<- NULL
results_poisson[[1]]<-NULL

# check = one_simulation(178)
# attributes(check$message_logistic)
# check$message_logistic$message

# Power to pick the best arm (we should specify the target arm/arms to pick for power analysis
# under null hypothesis, should pick every arm

arm_pick = sapply(results_raw, `[`, 1)
arm_pick_con = sapply(results_logistic, `[`, 3)
arm_pick_con_poi = sapply(results_poisson, `[`, 5)

# best_arm = (1:arms)[order(event.rates, decreasing=T)][1]
best_arm = 0

if(is.null(pick)){  # if not specified then pick only the best arm (need to remove null outputs)
  power_pick_exact_arm  = 1/(ntimes-error_raw)*sum(sapply(1:(ntimes-error_raw), function(i) {
    (length(arm_pick[[i]])==1) & (arm_pick[[i]][1]==best_arm)}))
  power_pick_exact_arm_con  = 1/(ntimes-bad_estimate_logistic) *sum(sapply(1:(ntimes-bad_estimate_logistic), function(i) {
    (length(arm_pick_con[[i]])==1) & (arm_pick_con[[i]][1]==best_arm)}))
  power_pick_exact_arm_con_poi  = 1/(ntimes-bad_estimate_poisson) * sum(sapply(1:(ntimes-bad_estimate_poisson), function(i) {
    (length(arm_pick_con_poi[[i]])==1) & (arm_pick_con_poi[[i]][1]==best_arm)}))
  
} else{ # if specified, pick out all arms specified or part of the arms specified
  power_pick_exact_arm = mean(sapply(1:(ntimes-error_raw), function(i) 
    floor(sum(arm_pick[[i]] %in% pick)/length(arm_pick[[i]]))), na.rm=T)
  power_pick_exact_arm_con = mean(sapply(1:(ntimes-bad_estimate_logistic), function(i) 
    floor(sum(arm_pick_con[[i]] %in% pick)/length(arm_pick_con[[i]]))),  na.rm=T)
  power_pick_exact_arm_con_poi = mean(sapply(1:(ntimes-bad_estimate_poisson), function(i) 
    floor(sum(arm_pick_con_poi[[i]] %in% pick)/length(arm_pick_con_poi[[i]]))), na.rm=T)
  
}
  
  
  

# mean of treatment ratios, for all pairwise
arm_effect = lapply(1:arms, function(i) sapply(sapply(results_raw, `[`, 2), `[`, i))
arm_effect_con = lapply(1:arms, function(i) sapply(sapply(results_logistic, `[`, 4), `[`, i))
arm_effect_con_poi = lapply(1:arms, function(i) sapply(sapply(results_poisson, `[`, 6), `[`, i))

tmp_index = lower.tri(matrix(0, arms, arms))
index = cbind(matrix(rep(1:arms, arms), byrow=T)[tmp_index], matrix(rep(1:arms, each=arms), byrow=T)[tmp_index])

summary_table = matrix(0,ncol=6, nrow=dim(index)[1], dimnames = list(NULL, c('arm_a', 'arm_b', 'estimated_ratio', 'estimated_ratio_logistic', 'estimated_ratio_poisson','true_ratio')))
summary_table[,1:2] <- index

for (k in 1:dim(index)[1]){
  summary_table[k,3] = mean(arm_effect[[index[k,1]]]/arm_effect[[index[k,2]]]) # mean of ratio, not ratio of mean here
  summary_table[k,4] = mean(arm_effect_con[[index[k,1]]]/arm_effect_con[[index[k,2]]])
  summary_table[k,5] = mean(arm_effect_con_poi[[index[k,1]]]/arm_effect_con_poi[[index[k,2]]])
  summary_table[k,6] = event.rates[index[k,1]]/event.rates[index[k,2]]
}

summary_table
  
power = c(power_pick_exact_arm, power_pick_exact_arm_con, power_pick_exact_arm_con_poi)
names(power) <-c('pairwise prop test power', 'logistic power', 'poisson power')

effect_table = matrix(c(sapply(1:arms, function(i)mean(arm_effect[[i]])),
sapply(1:arms, function(i)mean(arm_effect_con[[i]])),
sapply(1:arms, function(i)mean(arm_effect_con_poi[[i]])), event.rates), byrow=F, ncol=4, dimnames = list(c(1:arms), c('effect empirical', 'effect logistic', 'effect poisson', 'effect true')))

bad_estimation = c(error_raw/ntimes, bad_estimate_logistic/ntimes, bad_estimate_poisson/ntimes)
names(bad_estimation)<-c('error raw', 'logistic', 'poisson')
return(list('Effect estimate (arm ratio)' = summary_table,'Effect estimate (arm event rate)' = effect_table, 'power'= power,
            'Bad estimation case proportion' = bad_estimation))
} 


## we can run a total of N settings 
run.multiple.designs <- function(design.n, # length N vector
                                 design.rate, # N by a matrix, a is number of arms
                                 design.stages, # length N vector
                                 ntimes=500, 
                                 design.pick, # length N list 
                                 design.drift # length N vector
                                 ){
  N = length(design.n)
  res = lapply(1:N, function(i){
    totaln = design.n[i]
    event.rates= design.rate[i,]
    stages= design.stages[i]
    ntimes= ntimes
    pick = design.pick[[i]]
    drift = design.drift[i]
    settings = list(totaln = design.n[i],
                    event.rates= design.rate[i,],
                    stages= design.stages[i],
                    ntimes= ntimes,
                    pick = design.pick[[i]],
                    drift = design.drift[i])
    return(list(Settings = settings, Results = get.effect.power(totaln,event.rates,stages,ntimes,pick,drift)))
    })
  return(res)
  }

# see1 = get.effect.power(totaln=5000,event.rates=c(0.2, 0.03),stages=2,ntimes=200,pick=1,drift=4)
# see2 = get.effect.power(totaln=720,event.rates=c(0.2, 0.1),stages=4,ntimes=200,pick=1,drift=1.2)
# see3 = get.effect.power(totaln=720,event.rates=c(0.2, 0.1),stages=8,ntimes=200,pick=1,drift=1.2)



# some cases where the estimation could be very off (seems to happen with small sample)
# test_off_1 = get.effect.power(totaln=200,event.rates=c(0.2, 0.15, 0.1, 0.05),stages=4,ntimes=500,pick=1,drift=1.2)
# test_off_2 = get.effect.power(totaln=200,event.rates=c(0.2, 0.15, 0.1, 0.05),stages=4,ntimes=500,pick=1,drift=0.8)
# test_off_22 = get.effect.power(totaln=4000,event.rates=c(0.2, 0.15, 0.1, 0.05),stages=4,ntimes=500,pick=1,drift=0.6)

###############################################
# systematically run the simulations
################################################


# filepath = '/Users/Kun/Desktop/Dropbox/study/ind_amalia/'
filepath = '~/Desktop/amalia_results/'


# filepath = "//fs2-vip-nfs.nfs.biost.priv/students/yuek/Desktop/amalia_results/"
data = read.csv(paste0(filepath, 'null_2_arm.csv'))
# three groups, 2 arms/3 arms/4 arms
arm=2
group = 1:nrow(data)
set.seed(2)
res_2arm = run.multiple.designs(design.n = data$totalN.less[group],
                                design.rate = as.matrix(data[group, 2:(2+arm-1)]),
                                design.stages = rep(2, length(group)),
                                design.pick = NULL,
                                design.drift = data$drift[group],
                                ntimes=500
)
save.image(paste0(filepath, 'GLM_sandwich_extreme_2_arm_typeI.RData'))

# test = run.multiple.designs(design.n = c(720, 82, 720), # length N vector
#                      design.rate = matrix(c(0.2, 0.1, 
#                                             0.8, 0.4, 
#                                             0.2, 0.1), ncol=2, byrow=T), # N by a matrix, a is number of arms
#                      design.stages = rep(4, 3), # length N vector
#                      ntimes=100, 
#                      design.pick = list(1, 1, 1), # length N list 
#                      design.drift = c(1,1,1.2))
# 
# test = run.multiple.designs(design.n = c(720, 720, 720, 325, 325), # length N vector
#                             design.rate = matrix(c(0.2,0.1,
#                                                    0.2,0.1,
#                                                    0.2,0.1,
#                                                    0.6, 0.4,
#                                                    0.6, 0.4), 
#                                                  ncol=2, byrow=T), # N by a matrix, a is number of arms
#                             design.stages = rep(4, 5), # length N vector
#                             ntimes=500, 
#                             design.pick = list(1, 1, 1, 1, 1), # length N list 
#                             design.drift = c(1,
#                                              0.8,
#                                              1.2,
#                                              1,
#                                              1.1))
# 
# test = run.multiple.designs(design.n = c(4228,4228,4228,450,450), # length N vector
#                             design.rate = matrix(c(0.2,	0.15,	0.1,
#                                                    0.2,	0.15,	0.1,
#                                                    0.2,	0.15,	0.1,
#                                                    0.5,0.3,0.2,
#                                                    0.5,0.3,0.2), 
#                                                  ncol=3, byrow=T), # N by a matrix, a is number of arms
#                             design.stages = rep(4, 5), # length N vector
#                             ntimes=500, 
#                             design.pick = list(1, 1, 1, 1, 1), # length N list 
#                             design.drift = c(1,
#                                              0.8,
#                                              1.2,
#                                              1,
#                                              1.1))
# 
# 
# test = run.multiple.designs(design.n = c(2182, 2182, 2182, 578, 578), # length N vector
#                             design.rate = matrix(c(0.4,0.3,0.2,0.1,
#                                                    0.4,0.3,0.2,0.1,
#                                                    0.4,0.3,0.2,0.1,
#                                                    0.5,0.3,0.2,0.1,
#                                                    0.5,0.3,0.2,0.1), 
#                                                    ncol=4, byrow=T), # N by a matrix, a is number of arms
#                             design.stages = rep(4, 5), # length N vector
#                             ntimes=500, 
#                             design.pick = list(1, 1, 1, 1, 1), # length N list 
#                             design.drift = c(1,
#                                              0.8,
#                                              1.2,
#                                              1,
#                                              1.1))
# check_res = run.multiple.designs(design.n = c(720,
#                                               325,
#                                               271,
#                                               1156),
#                                  design.rate = matrix(c(0.2,	0.1,
#                                                         0.6,	0.4,
#                                                         0.8,	0.6,
#                                                         0.4,	0.3), byrow=T, ncol=2),
#                                  design.stages = rep(2, 4),
#                                  design.pick = NULL,
#                                  design.drift = rep(1, 4),
#                                  ntimes=500
# )



