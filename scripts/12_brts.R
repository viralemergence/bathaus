# BRT analysis for bat viruses
# babetke@utexas.edu
# updated 07/24/24

# clear environment and graphics
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
# library(fastDummies) # dummy coding family this is done in data cleaning script now
library(gbm)
library(rsample)
library(ROCR)
library(caret) 
library(InformationValue)

## load data with logged predictors
setwd("/Users/brianabetke/Desktop/bathaus")
data <- readRDS("flat files/Log cleaned dataset 30 cutoff.rds")

# how many NAs are there - 234
length(data$Synurbic[is.na(data$Synurbic)])

# # zoonotic virus proportions
# data <- data %>% 
#   mutate(zoo_prop = zvirus/virus) %>%
#   mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop))
# 
# # reservoir status
# data <- data %>%
#   mutate(dum_virus = if_else(virus <= 0, 0, 1),
#          dum_zvirus = if_else(zvirus <= 0, 0, 1)) 

# # transform zoonotic proportion to gaussian space here so it is included 
# # in both datasets
# c <- data[c("species","zoo_prop")]
# c$xi <- c$zoo_prop * 100
# c$ni <- 100
# 
# # Freeman tukey double arcsine transformation
# t <- metafor::escalc(measure = "PFT", xi=xi, ni=ni, data=c)
# # t$back <- metafor::transf.ipft(t$yi, t$ni) # back transformation
# 
# # subset 
# tsub <- t[c("species","yi")]
# 
# # rename yi to zoo_pft
# names(tsub) <- c("species","zoo_pft")
# 
# # merge to data by species
# data <- merge(data, tsub, by = "species")
# 
# # remove
# rm(c,t,tsub)

# remove redundant diet variable
data$det_inv <- NULL

# save before removing species (need this for predictions)
fdata <- data

# remove species
data$species <- NULL

# run tuning grid
# ifelse for running grid search function
gsrun = "no" 
if(gsrun == "yes"){# run grid search 
  
  # Set up BRT tuning via search grid
  # function to make grids?
  ## hyperparameter grid, maybe allow the number of seeds to change?
  makegrid <- function(seed,trees) {
    
    # create grid
    tgrid <- expand.grid(n.trees = trees,
                         interaction.depth = c(2, 3, 4),
                         shrinkage = c(0.01, 0.001, 0.0005),
                         #n.minobsinnode = c(5, 10, 15), 
                         #bag.fraction = c(0.5, 0.8, 1),
                         seed = seq(1,seed,by = 1))
    
    # adjust combinations with low n.trees and high shrinkage
    tgrid <- mutate(tgrid, n.trees = ifelse(shrinkage < 0.01 & n.trees == 5000, 15000, n.trees))
    
    ## trees, depth, shrink, min, prop
    tgrid$id <- with(tgrid,paste(n.trees,interaction.depth,shrinkage))
    
    ## sort by id then seed
    tgrid <- tgrid[order(tgrid$id,tgrid$seed),]
    
    # Number rows
    tgrid$row <- 1:nrow(tgrid)
    tgrid$id2=factor(as.numeric(factor(tgrid$id)))
    
    return(tgrid)
  }
  
  #### Function to assess each hyperparameter combination
  # this function takes a search grid and runs them through gbms for a given split of data
  # then assesses performance of each combination of parameters to get the optimal parameters for
  # final gbms.
  
  # Takes:
  # row <- row in the hgrid
  # data_df <- dataframe without species
  # hgrid <- change grid based on response
  # response <- indicate which response, determines the error distribution used
  # cv <- null to avoid warnings when running poisson or gaussian dist. Set to TRUE for class stratify.
  grid_search <- function(row, data_df, hgrid, response, cv = NULL){
    
    # new data
    ndata <- data_df
    
    # correct the response and remove raw response variables
    ndata <- ndata %>% 
      mutate(response = !!sym(response)) %>%
      select(-c("Virus", "dum_virus", "dum_zvirus", "zoo_sp", "vfam", "zfam"))
    
    # assign name of hgrid
    hgrid <- hgrid
    
    ## use rsample to split (allow the proportion to be changed incase)
    set.seed(hgrid$seed[row])
    split <- initial_split(ndata,prop=0.8,strata="response")
    
    ## test and train
    dataTrain <- training(split)
    dataTest <- testing(split)
    
    ## yTest and yTrain
    yTrain <- dataTrain$response
    yTest <- dataTest$response
    
    # ifelse
    dist <- ifelse(response == "vfam", "poisson", 
                   ifelse(response == "zfam", "poisson", 
                          ifelse(response == "dum_virus", "bernoulli", 
                                 ifelse(response == "dum_zvirus", "bernoulli", "error"))))
    
    ## BRT
    set.seed(1)
    gbmOut=gbm(response ~ . ,data=dataTrain,
               n.trees=hgrid$n.trees[row],
               distribution=dist,
               shrinkage=hgrid$shrinkage[row],
               interaction.depth=hgrid$interaction.depth[row],
               n.minobsinnode=10,
               cv.folds=10,
               class.stratify.cv=cv, # will throw warning for poisson and gaussian unless = NULL
               bag.fraction=0.5,
               train.fraction=1,
               n.cores=1,
               verbose=F)
    
    ## performance
    par(mfrow=c(1,1),mar=c(4,4,1,1))
    best.iter=gbm.perf(gbmOut,method="cv")
    
    ## predict with test data
    preds <- predict(gbmOut,dataTest,n.trees=best.iter,type="response")
    
    ## known
    result <- dataTest$response
    
    if(dist == "bernoulli"){
      
      ## sensitiviy and specificity
      sen <- InformationValue::sensitivity(result,preds)
      spec <- InformationValue::specificity(result,preds)
      
      ## AUC on train
      auc_train <- gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
      
      ## AUC on test
      auc_test <- gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
      
      ## print
      print(paste("hpar row ",row," done; test AUC is ",auc_test, " for ", response, sep=""))
      
      ## save outputs
      return(list(best=best.iter,
                  trainAUC=auc_train,
                  testAUC=auc_test,
                  spec=spec,
                  sen=sen,
                  wrow=row))
      
    } else {
      
      # another if else statement or try to transform rmse after? 
      # since it will be in the unit of the response.
      
      # calculate RMSE
      testMSE <- mean((preds - result)^2)
      testRMSE <- sqrt(testMSE)
      
      # training RMSE
      preds <- predict(gbmOut,dataTrain,n.trees=best.iter,type="response")
      trainMSE <- mean((preds - yTrain)^2)
      trainRMSE <- sqrt(trainMSE)
      
      # print
      print(paste("hpar row ",row," done; Test RMSE is ", testRMSE, " for ", response, sep=""))
      
      # save outputs
      return(list(best = best.iter,
                  trainMSE = trainMSE,
                  trainRMSE = trainRMSE,
                  testMSE = testMSE,
                  testRMSE = testRMSE,
                  wrow = row))
    }
    
  }
  
  # grids for each model type
  pgrid <- makegrid(10, 5000)  # smaller grid for poisson models (seem to have low no of trees)
  bgrid <- makegrid(10, c(5000, 25000)) # binary models
  
  # ### trim just for testing
  # pgrid <- pgrid[71, ]
  # pgrid$n.trees <- 5000
  # pgrid$shrinkage <- 0.001
  # pgrid$row <- 1
  
  ### run for all types of responses
  # Virus Richness and transformed zoonotic proportion
  vpars <- lapply(1:nrow(pgrid),function(x) grid_search(x, data_df = data, hgrid = pgrid, response="vfam"))
  zpars <- lapply(1:nrow(pgrid),function(x) grid_search(x, data_df = data, hgrid = pgrid, response="zfam"))
  
  # trim just for testing
  # bgrid <- bgrid[84, ]
  
  # Reservoir status models
  vbinpars <- lapply(1:nrow(bgrid),function(x) grid_search(x, data_df = data, hgrid = bgrid, response="dum_virus", cv = TRUE))
  zbinpars <- lapply(1:nrow(bgrid),function(x) grid_search(x, data_df = data, hgrid = bgrid, response="dum_zvirus", cv = TRUE))
  
  ## get virus results
  vresults <- data.frame(sapply(vpars,function(x) x$trainRMSE),
                         sapply(vpars,function(x) x$testRMSE),
                         sapply(vpars,function(x) x$wrow),
                         sapply(vpars,function(x) x$best))
  names(vresults) <- c("trainRMSE","testRMSE","row","best")
  
  # Merge with hgrid
  vcomplete <- merge(vresults, pgrid, by = "row")
  
  # get zoonotic results
  zresults <- data.frame(sapply(zpars,function(x) x$trainRMSE),
                          sapply(zpars,function(x) x$testRMSE),
                          sapply(zpars,function(x) x$wrow),
                          sapply(zpars,function(x) x$best))
  names(zresults) <- c("trainRMSE","testRMSE","row","best")
  
  # Merge with hgrid
  zcomplete <- merge(zresults, pgrid, by = "row")
  
  # Overall virus binary models
  vbresults <- data.frame(sapply(vbinpars,function(x) x$trainAUC),
                          sapply(vbinpars,function(x) x$testAUC),
                          sapply(vbinpars,function(x) x$spec),
                          sapply(vbinpars,function(x) x$sen),
                          sapply(vbinpars,function(x) x$wrow),
                          sapply(vbinpars,function(x) x$best))
  names(vbresults) <- c("trainAUC","testAUC","spec","sen","row","best")
  
  # Merge with hgrid
  vbcomplete <- merge(vbresults, bgrid, by = "row")
  
  # zoonotic virus reservoir models
  ## get virus results
  zbresults <- data.frame(sapply(zbinpars,function(x) x$trainAUC),
                          sapply(zbinpars,function(x) x$testAUC),
                          sapply(zbinpars,function(x) x$spec),
                          sapply(zbinpars,function(x) x$sen),
                          sapply(zbinpars,function(x) x$wrow),
                          sapply(zbinpars,function(x) x$best))
  names(zbresults) <- c("trainAUC","testAUC","spec","sen","row","best")
  
  # Merge with hgrid
  zbcomplete <- merge(zbresults, bgrid, by = "row")
  
  # write as csv
  setwd("/Users/brianabetke/Desktop/bathaus")
  write_csv(vcomplete, "flat files/virus family grid search.csv")
  write_csv(zcomplete, "flat files/zoonotic family grid search.csv")
  write_csv(vbcomplete, "flat files/virus binary grid search.csv")
  write.csv(zbcomplete, "flat files/zoonotic virus binary grid search.csv")
  
} else {# read in grid search results
  
  # overall_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/virus grid search.csv")
  # zoo_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/pft zoonotic grid search.csv")
  # vres_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/virus reservoir grid search.csv")
  # zres_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/zoonotic virus reservoir grid search.csv")
 
  # lab comp paths
  overall_search <- read_csv("flat files/virus family grid search.csv")
  zoo_search <- read_csv("flat files/zoonotic family grid search.csv")
  vbin_search <- read_csv("flat files/virus binary grid search.csv")
  zbin_search <- read_csv("flat files/zoonotic virus binary grid search.csv")
}

### plots of parameters - figure S2
library(patchwork)

#r <- overall_search %>% filter(shrinkage < 0.1)

# plots for family richness models
virus_gg <- ggplot(overall_search, aes(x = factor(shrinkage), y = testRMSE)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = "Test RMSE", fill = "Interaction Depth", title = "Virus Family Richness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

zvirus_gg <- ggplot(zoo_search, aes(x = factor(shrinkage), y = testRMSE)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = NULL, fill = "Interaction Depth", title = "Zoonotic Family Richness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

# overall virus reservoir
vresAUC_gg <- ggplot(vbin_search, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = "AUC", fill = "Interaction Depth", title = "Virus Hosting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

vresSen_gg <- ggplot(vbin_search, aes(x = factor(shrinkage), y = sen)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = "Sensitiviy", fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

vresSpec_gg <- ggplot(vbin_search, aes(x = factor(shrinkage), y = spec)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning rate", y = "Specificity", fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

# zoonotic reservoir
zresAUC_gg <- ggplot(zbin_search, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Interaction Depth", title = "Zoonotic Hosting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

zresSen_gg <- ggplot(zbin_search, aes(x = factor(shrinkage), y = sen)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

zresSpec_gg <- ggplot(zbin_search, aes(x = factor(shrinkage), y = spec)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning rate", y = NULL, fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

png("figs/figure S2.png", width=5, height=7,units="in",res=300)
guide_area () + 
  (virus_gg + zvirus_gg) / (vresAUC_gg + zresAUC_gg) / 
  (vresSen_gg + zresSen_gg) / (vresSpec_gg + zresSpec_gg) +
  plot_layout(guides = "collect", heights = c(1, 15)) & theme(legend.position = "top", legend.text = element_text(size = 8), legend.title = element_text(size = 8))
dev.off()

# remove plots from envrionment
rm(virus_gg, zvirus_gg, vresAUC_gg, vresSen_gg, vresSpec_gg, zresAUC_gg, zresSen_gg, zresSpec_gg)

# unload patchwork
detach("package:patchwork", unload = TRUE)

# look at summary of ids across seeds, arranged by lowest testRMSE
overall_search %>% 
  group_by(id) %>% 
  summarise(med = median(testRMSE)) %>% 
  arrange(med)

# look at summary of ids across seeds, arranged by lowest testRMSE
zoo_search %>% 
  group_by(id) %>% 
  summarise(med = median(testRMSE)) %>% 
  arrange(med)

# summarize
vbin_search %>% 
  group_by(id) %>% 
  summarise(medAUC = median(testAUC),
            medSen = median(sen),
            medSpec = median(spec)) %>% 
  arrange(desc(medAUC))

# summarize
zbin_search %>% 
  group_by(id) %>% 
  summarise(medAUC = median(testAUC),
            medSen = median(sen),
            medSpec = median(spec)) %>% 
  arrange(desc(medAUC))

# look at best iteration of trees
# filter by best perf
overall_search %>% filter(id == "5000 4 0.01") # doesn't exceed 2000
zoo_search %>% filter(id == "5000 4 0.01") # doesn't exceed 1000
vbin_search %>% filter(id == "5000 2 0.01") # doesn't exceed 1500
zbin_search %>% filter(id == "5000 3 0.01") # close to maxing out at 1000 

# clean
rm(overall_search, vbin_search, zbin_search, zoo_search)

#### Define brt function 
# seed - the seeds to be used for the unique splits
# response - indicating what the response variable is
# nt - number of trees based off search grid
# shr - optimal shrinkage value based on search grid
# int.d - interaction depth value based on search grid
# syn - to indicate if it is a model without synurbic as a predictor. If == "yes" then keep it, remove if anything else
# cv - change class stratify to NULL or TRUE depending on response variable (NULL for poisson and gaussian, TRUE for bernoulli)
# this model has been expanded to run BRTs for all response types.
brts <- function(seed, response, nt, shr, int.d, syn, cv = NULL){
  
  # new data
  ndata <- data
  
  if(syn == "yes"){
    # correct the response and remove raw response variables
    ndata <- ndata %>% 
      mutate(response = !!sym(response)) %>%
      select(-c("Virus", "dum_virus", "dum_zvirus", "zoo_sp", "vfam", "zfam"))
    
  }else{ # remove synurbic from data
    # correct the response and remove raw response variables
    ndata <- ndata %>% 
      mutate(response = !!sym(response)) %>%
      select(-c("Virus", "dum_virus", "dum_zvirus", "zoo_sp", "vfam", "zfam", "Synurbic"))
    
    #fdata$Synurbic <- NULL # make null in fdata for prediction? Not sure if this is necessary if the formula dictates variables
  }
  
  # ifelse to determine distribution. If not any of the below options, throw error
  dist <- ifelse((response == "vfam" | response == "zfam" | response == "log_cites" | response == "log_vcites"), "poisson", 
                 ifelse((response == "dum_zvirus" | response == "dum_virus"),"bernoulli", 
                        ifelse(response == "zoo_pft","gaussian","Error")))
  
  # remove citations if predicting citations
  if(response == "log_cites" | response == "log_vcites") {
    ndata$response = as.integer(exp(ndata$response)-1) # backtransform to count
    ndata$log_cites = NULL
    ndata$log_vcites = NULL
    # } else if(response == "vcites"){
    #   ndata$vcites = NULL
  } else {
    ndata = ndata
  }
  
  ## use rsample to split
  set.seed(seed)
  split <- initial_split(ndata,prop=0.8,strata="response")
  
  ## test and train for data
  train <- training(split)
  test <- testing(split)
  
  # pull response test and train
  yTrain <- train$response
  yTest <- test$response
  
  ## BRT
  set.seed(1)
  gbmOut <- gbm(response ~ . ,data=train,
                n.trees=nt,
                distribution=dist,
                shrinkage=shr,
                interaction.depth=int.d,
                n.minobsinnode=4,
                cv.folds=10,
                class.stratify.cv=cv,
                bag.fraction=0.5,
                train.fraction=1,
                n.cores=1,
                verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter <- gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds <- predict(gbmOut,test,n.trees=best.iter,type="response")
  
  # output dependent on dist
  # AUC, sen, spec, ROC for binary models
  if(dist == "bernoulli"){ 
    
    ## known 
    result <- test$response
    
    ## sensitiviy and specificity
    sen <- InformationValue::sensitivity(result,preds)
    spec <- InformationValue::specificity(result,preds)
    
    ## AUC on train
    auc_train <- gbm.roc.area(yTrain,predict(gbmOut,train,n.trees=best.iter,type="response"))
    
    ## AUC on test
    auc_test <- gbm.roc.area(yTest,predict(gbmOut,test,n.trees=best.iter,type="response"))
    
    ## inner loop if yTest is all 0
    ## I removed the portion for cites but this portion seems important for when variance of ytest is zero?
    if(var(yTest)==0){
      
      perf=NA
      
    }else{
      
      ## ROC
      pr=prediction(preds,test$response)
      perf=performance(pr,measure="tpr",x.measure="fpr")
      perf=data.frame(perf@x.values,perf@y.values)
      names(perf)=c("fpr","tpr")
      
      ## add seed
      perf$seed=seed
      
    }
    
    ## relative importance
    bars <- summary(gbmOut,n.trees=best.iter,plotit=F)
    bars$rel.inf <- round(bars$rel.inf,2)
    
    ## predict with cites
    preds <- predict(gbmOut,fdata,n.trees=best.iter,type="response")
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","vfam", "zfam", "Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # then mean cites
    pdata <- fdata
    pdata$log_cites <- mean(pdata$log_cites)
    pdata$log_vcites <- mean(pdata$log_vcites)
    pred_data$cpred <- predict(gbmOut, pdata, n.trees=best.iter, type="response")
    
    ## print
    print(paste("BRT ", seed," done; test AUC = ", auc_test, " for ", response, sep=""))
    
    ## save outputs
    return(list(mod = gbmOut,
                best = best.iter,
                trainAUC = auc_train,
                testAUC = auc_test,
                spec = spec,
                sen = sen,
                roc = perf,
                rinf = bars,
                predict = pred_data,
                traindata = train,
                testdata = test,
                seed = seed))
    
  }else{ # pseudo R2 for poisson and gaussian models
    
    ## relative importance
    bars <- summary(gbmOut,n.trees=best.iter,plotit=F)
    bars$rel.inf <- round(bars$rel.inf,2)
    
    # calculate RMSE
    testMSE <- mean((preds - yTest)^2)
    testRMSE <- sqrt(testMSE)
    
    # training RMSE
    preds <- predict(gbmOut,train,n.trees=best.iter,type="response")
    trainMSE <- mean((preds - yTrain)^2)
    trainRMSE <- sqrt(trainMSE)
    
    # calculate train pseudo R squared
    r2train <- (1-(sum((yTrain - predict(gbmOut, newdata=train, n.trees=best.iter, type='response'))^2)/
                     sum((yTrain - mean(yTrain))^2)))
    
    # calculate test pseudo R squared
    r2test <- (1-(sum((yTest - predict(gbmOut, newdata=test, n.trees=best.iter, type='response'))^2)/
                    sum((yTest - mean(yTest))^2)))
    
    ## predict with cites
    preds <- predict(gbmOut,fdata,n.trees=best.iter,type="response")
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","vfam","zfam","Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # then mean cites
    pdata <- fdata
    pdata$log_cites <- mean(pdata$log_cites)
    pdata$log_vcites <- mean(pdata$log_vcites)
    pred_data$cpred <- predict(gbmOut, pdata, n.trees=best.iter, type="response")
    
    if(dist == "poisson"){
      
      # back transforming predictions - don't if count models
      pred_data$ni <- NULL
      pred_data$back <- NULL
      
      # sort
      pred_data <- pred_data[order(pred_data$pred,decreasing=T),]
      
    } else { # back transformations for predictions to proportion
      
      # back transforming predictions
      pred_data$ni <- 100
      pred_data$back <- metafor::transf.ipft(pred_data$pred, pred_data$ni)
      
      # back transform cpred
      pred_data$cback <- metafor::transf.ipft(pred_data$cpred, 100)
      
      # sort
      pred_data <- pred_data[order(pred_data$back,decreasing=T),]
      
    }
    
    ## print
    print(paste("BRT ",seed," done; test R2 = ", r2test, " for ", response, sep=""))
    
    ## save outputs
    return(list(mod = gbmOut,
                best = best.iter,
                testRMSE = testRMSE,
                trainRMSE = trainRMSE,
                testr2 = r2test,
                trainr2 = r2train,
                rinf = bars,
                predict = pred_data,
                traindata = train,
                testdata = test,
                seed = seed))
  }
  
}

##### apply across specified number of splits smax
smax=50
#### transformed versions
vfams_brts <- lapply(1:smax,function(x) brts(seed = x,response = "vfam", nt = 2500, shr = 0.01, int.d = 4, syn = "yes", cv = NULL))
no_vfams_brts <- lapply(1:smax,function(x) brts( seed = x,response = "vfam", nt = 2500, shr = 0.01, int.d = 4, syn = "no", cv = NULL))

# save then remove
saveRDS(vfams_brts,"flat files/virus family richness with brts.rds")
saveRDS(no_vfams_brts,"flat files/virus family richness without brts.rds")

rm(vfams_brts, no_vfams_brts)

# Proportion with and without
zfams_brts <- lapply(1:smax,function(x) brts(seed = x,response = "zfam", nt = 2500, shr = 0.01, int.d = 4, syn = "yes", cv = NULL))
no_zfams_brts <- lapply(1:smax,function(x) brts(seed = x,response = "zfam", nt = 2500, shr = 0.01, int.d = 4, syn = "no", cv = NULL))

saveRDS(zfams_brts, "flat files/zoonotic family richness with brts.rds")
saveRDS(no_zfams_brts, "flat files/zoonotic family richness without brts.rds")

rm(zfams_brts,no_zfams_brts)

# Overall virus reservoir status
vbin_brts <- lapply(1:smax,function(x) brts(seed = x,response = "dum_virus", nt = 2500, shr = 0.01, int.d = 2, syn = "yes", cv = TRUE))
no_vbin_brts <- lapply(1:smax,function(x) brts(seed = x,response ="dum_virus", nt = 2500, shr = 0.01, int.d = 2, syn = "no", cv = TRUE))

saveRDS(vbin_brts, "flat files/binary virus with brts.rds")
saveRDS(no_vbin_brts, "flat files/binary virus without brts.rds")

rm(vbin_brts,no_vbin_brts)

# Zoonotic virus reservoir
zbin_brts <- lapply(1:smax,function(x) brts(seed = x,response = "dum_zvirus",  nt = 2500, shr = 0.01, int.d = 3, syn = "yes", cv = TRUE))
no_zbin_brts <- lapply(1:smax,function(x) brts(seed = x,response ="dum_zvirus",  nt = 2500, shr = 0.01, int.d = 3, syn = "no", cv = TRUE))

saveRDS(zbin_brts, "flat files/zoonotic binary with brts.rds")
saveRDS(no_zbin_brts, "flat files/zoonotic binary without brts.rds")

rm(zbin_brts,no_zbin_brts)

# citation models
smax=20 # no. seeds for citation models
cite_brts <- lapply(1:smax,function(x) brts(seed = x,response = "log_cites", nt = 25000, shr = 0.001, int.d = 3, syn = "yes", cv = NULL))
vcite_brts <- lapply(1:smax,function(x) brts(seed = x,response = "log_vcites", nt = 25000, shr = 0.001, int.d = 3, syn = "yes", cv = NULL))

saveRDS(cite_brts, "flat files/citation brts.rds")
saveRDS(vcite_brts, "flat files/virus citation brts.rds")
