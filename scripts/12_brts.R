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
  # cv <- null to avoid warnings when running poisson or gaussian dist. Set to TRUE for class stratified cv.
  grid_search <- function(row, data_df, hgrid, response, cv = NULL){
    
    # new data
    ndata <- data_df
    
    # correct the response and remove raw response variables
    ndata <- ndata %>% 
      mutate(response = !!sym(response)) %>%
      select(-c("Virus", "dum_virus", "zVirus", "dum_zvirus"))
    
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
    dist <- ifelse(response == "Virus", "poisson", 
                   ifelse(response == "zVirus", "poisson", 
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
  vpars <- lapply(1:nrow(pgrid),function(x) grid_search(x, data_df = data, hgrid = pgrid, response="Virus"))
  zpars <- lapply(1:nrow(pgrid),function(x) grid_search(x, data_df = data, hgrid = pgrid, response="zVirus"))
  
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
  write_csv(vcomplete, "flat files/virus species grid search.csv")
  write_csv(zcomplete, "flat files/zoonotic species grid search.csv")
  write_csv(vbcomplete, "flat files/virus binary grid search.csv")
  write_csv(zbcomplete, "flat files/zoonotic virus binary grid search.csv")
  
  # rename for consistency with rest of code
  overall_search <- vcomplete
  zoo_search <- zcomplete
  vbin_search <- vbcomplete
  zbin_search <- zbcomplete
  
  # remove dfs with old names
  rm(vcomplete, zcomplete, vbcomplete, zbcomplete)
  
} else {# read in grid search results
  
  # overall_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/virus grid search.csv")
  # zoo_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/pft zoonotic grid search.csv")
  # vres_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/virus reservoir grid search.csv")
  # zres_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/zoonotic virus reservoir grid search.csv")
 
  # lab comp paths
  overall_search <- read_csv("flat files/virus species grid search.csv")
  zoo_search <- read_csv("flat files/zoonotic species grid search.csv")
  vbin_search <- read_csv("flat files/virus binary grid search.csv")
  zbin_search <- read_csv("flat files/zoonotic virus binary grid search.csv")
}

### plots of parameters - figure S2
library(patchwork)

# plots for family richness models
virus_gg <- ggplot(overall_search, aes(x = factor(shrinkage), y = testRMSE)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = "Test RMSE", fill = "Interaction Depth", title = "Virus Species Richness") +
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
  labs(x = "Learning Rate", y = NULL, fill = "Interaction Depth", title = "Zoonotic Species Richness") +
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
zoo_search %>% filter(id == "5000 3 0.01") # doesn't exceed 1000
vbin_search %>% filter(id == "5000 2 0.01") # doesn't exceed 1500
zbin_search %>% filter(id == "5000 4 0.01") # close to maxing out at 1000 

# clean
rm(overall_search, vbin_search, zbin_search, zoo_search)

#### Define brt function 
# seed - the seeds to be used for the unique splits
# response - indicating what the response variable is
# nt - number of trees based off search grid
# shr - optimal shrinkage value based on search grid
# int.d - interaction depth value based on search grid
# syn - to indicate if it is a model without synurbic as a predictor. If == "yes" then keep it, remove if anything else
# cv_strat - change class stratify to NULL or TRUE depending on response variable (NULL for poisson and gaussian, TRUE for bernoulli)
# hvar - the type of data being used ("full" for all variables, "no syn" for no roosting variable, "no hd" for no human density variables,
# and "no all" for no roosting or human density variables.)
brts <- function(seed, df, response, nt, shr, int.d, cv_strat = NULL, hvar){
  
  # new data
  ndata <- df
  
  # correct the response and remove raw response variables
  ndata <- ndata %>% 
    mutate(response = !!sym(response)) %>%
    select(-c("Virus", "zVirus", "dum_virus", "dum_zvirus"))
  
  # ifelse to determine distribution. If not any of the below options, throw error
  dist <- ifelse((response == "Virus" | response == "zVirus" | response == "log_cites" | response == "log_vcites"), "poisson", 
                 ifelse((response == "dum_zvirus" | response == "dum_virus"),"bernoulli", "Error"))
  
  # remove citations if predicting citations
  if(response == "log_cites" | response == "log_vcites") {
    ndata$response = as.integer(exp(ndata$response)-1) # backtransform to count
    ndata$log_cites = NULL
    ndata$log_vcites = NULL
    
  } else {
    
    # if not citation model, keep data the same
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
                class.stratify.cv=cv_strat,
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
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","Virus", "zVirus", "Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # then mean cites
    pdata <- fdata
    pdata$log_cites <- mean(pdata$log_cites)
    pdata$log_vcites <- mean(pdata$log_vcites)
    pred_data$cpred <- predict(gbmOut, pdata, n.trees=best.iter, type="response")
    
    ## print
    print(paste("BRT ", seed," done; test AUC = ", auc_test, " for ", response, " ", hvar, sep=""))
    
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
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","Virus","zVirus","Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # then mean cites
    pdata <- fdata
    pdata$log_cites <- mean(pdata$log_cites)
    pdata$log_vcites <- mean(pdata$log_vcites)
    pred_data$cpred <- predict(gbmOut, pdata, n.trees=best.iter, type="response")
    
    # sort
    pred_data <- pred_data[order(pred_data$pred,decreasing=T),]
    
    # if(dist == "poisson"){
    #   
    #   # back transforming predictions - don't if count models
    #   pred_data$ni <- NULL
    #   pred_data$back <- NULL
    #   
    #   # sort
    #   pred_data <- pred_data[order(pred_data$pred,decreasing=T),]
    #   
    # } else { # back transformations for predictions to proportion
    #   
    #   # back transforming predictions
    #   pred_data$ni <- 100
    #   pred_data$back <- metafor::transf.ipft(pred_data$pred, pred_data$ni)
    #   
    #   # back transform cpred
    #   pred_data$cback <- metafor::transf.ipft(pred_data$cpred, 100)
    #   
    #   # sort
    #   pred_data <- pred_data[order(pred_data$back,decreasing=T),]
    #   
    # }
    
    ## print
    print(paste("BRT ",seed," done; test R2 = ", r2test, " for ", response, " ", hvar, sep=""))
    
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

# full models
vsp_brts <- lapply(1:smax,function(x) brts(seed = x, df = data, response = "Virus", nt = 2500, shr = 0.01, int.d = 4, cv_strat = NULL, hvar = "full"))
zsp_brts <- lapply(1:smax,function(x) brts(seed = x, df = data, response = "zVirus", nt = 2500, shr = 0.01, int.d = 3, cv_strat = NULL, hvar = "full"))
vbin_brts <- lapply(1:smax,function(x) brts(seed = x, df = data, response = "dum_virus", nt = 2500, shr = 0.01, int.d = 2, cv_strat = TRUE, hvar = "full"))
zbin_brts <- lapply(1:smax,function(x) brts(seed = x, df = data, response = "dum_zvirus",  nt = 2500, shr = 0.01, int.d = 4, cv_strat = TRUE, hvar = "full"))

# save
saveRDS(vsp_brts,"flat files/virus species richness full brts.rds")
saveRDS(zsp_brts, "flat files/zoonotic species richness full brts.rds")
saveRDS(vbin_brts, "flat files/binary virus full brts.rds")
saveRDS(zbin_brts, "flat files/zoonotic binary full brts.rds")

# remove
rm(vsp_brts, zsp_brts, vbin_brts, zbin_brts)

# models without anthro
data %>% select(!Synurbic) -> data_no_syn

vsp_no_syn_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_syn, response = "Virus", nt = 2500, shr = 0.01, int.d = 4, cv_strat = NULL, hvar = "no_syn"))
zsp_no_syn_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_syn, response = "zVirus", nt = 2500, shr = 0.01, int.d = 3, cv_strat = NULL, hvar = "no_syn"))
vbin_no_syn_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_syn, response = "dum_virus", nt = 2500, shr = 0.01, int.d = 2, cv_strat = TRUE, hvar = "no_syn"))
zbin_no_syn_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_syn, response = "dum_zvirus",  nt = 2500, shr = 0.01, int.d = 4, cv_strat = TRUE, hvar = "no_syn"))

# save
saveRDS(vsp_no_syn_brts,"flat files/virus species richness sans syn brts.rds")
saveRDS(zsp_no_syn_brts, "flat files/zoonotic species richness sans syn brts.rds")
saveRDS(vbin_no_syn_brts, "flat files/binary virus sans syn brts.rds")
saveRDS(zbin_no_syn_brts, "flat files/zoonotic binary sans syn brts.rds")

# remove
rm(vsp_no_syn_brts, zsp_no_syn_brts, vbin_no_syn_brts, zbin_no_syn_brts, data_no_syn)

# models without human density variables
data %>% select(!c("log_X27.1_HuPopDen_Min_n.km2", "log_X27.2_HuPopDen_Mean_n.km2", "log_X27.3_HuPopDen_5p_n.km2","X27.4_HuPopDen_Change")) -> data_no_hd

vsp_no_hd_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_hd, response = "Virus", nt = 2500, shr = 0.01, int.d = 4, cv_strat = NULL, hvar = "no_hd"))
zsp_no_hd_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_hd, response = "zVirus", nt = 2500, shr = 0.01, int.d = 3, cv_strat = NULL, hvar = "no_hd"))
vbin_no_hd_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_hd, response = "dum_virus", nt = 2500, shr = 0.01, int.d = 2, cv_strat = TRUE, hvar = "no_hd"))
zbin_no_hd_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_hd, response = "dum_zvirus",  nt = 2500, shr = 0.01, int.d = 4, cv_strat = TRUE, hvar = "no_hd"))

# save
saveRDS(vsp_no_hd_brts,"flat files/virus species richness sans hd brts.rds")
saveRDS(zsp_no_hd_brts, "flat files/zoonotic species richness sans hd brts.rds")
saveRDS(vbin_no_hd_brts, "flat files/binary virus sans hd brts.rds")
saveRDS(zbin_no_hd_brts, "flat files/zoonotic binary sans hd brts.rds")

# remove
rm(vsp_no_hd_brts, zsp_no_hd_brts, vbin_no_hd_brts, zbin_no_hd_brts, data_no_hd)

# models without human density variables and synurbic
data %>% select(!c("Synurbic", "log_X27.1_HuPopDen_Min_n.km2", "log_X27.2_HuPopDen_Mean_n.km2", "log_X27.3_HuPopDen_5p_n.km2","X27.4_HuPopDen_Change")) -> data_no_all

vsp_no_all_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_all, response = "Virus", nt = 2500, shr = 0.01, int.d = 4, cv_strat = NULL, hvar = "no all"))
zsp_no_all_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_all, response = "zVirus", nt = 2500, shr = 0.01, int.d = 3, cv_strat = NULL, hvar = "no all"))
vbin_no_all_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_all, response = "dum_virus", nt = 2500, shr = 0.01, int.d = 2, cv_strat = TRUE, hvar = "no all"))
zbin_no_all_brts <- lapply(1:smax,function(x) brts(seed = x, df = data_no_all, response = "dum_zvirus",  nt = 2500, shr = 0.01, int.d = 4, cv_strat = TRUE, hvar = "no all"))

# save
saveRDS(vsp_no_all_brts,"flat files/virus species richness sans all brts.rds")
saveRDS(zsp_no_all_brts, "flat files/zoonotic species richness sans all brts.rds")
saveRDS(vbin_no_all_brts, "flat files/binary virus sans all brts.rds")
saveRDS(zbin_no_all_brts, "flat files/zoonotic binary sans all brts.rds")

rm(vsp_no_all_brts, zsp_no_all_brts, vbin_no_all_brts, zbin_no_all_brts, data_no_all)

# citation models
smax=20 # no. seeds for citation models
cite_brts <- lapply(1:smax,function(x) brts(seed = x, df = data, response = "log_cites", nt = 2500, shr = 0.01, int.d = 4, cv_strat = NULL, hvar = "full"))
vcite_brts <- lapply(1:smax,function(x) brts(seed = x, df = data, response = "log_vcites", nt = 2500, shr = 0.01, int.d = 4, cv_strat = NULL, hvar = "full"))

saveRDS(cite_brts, "flat files/citation brts.rds")
saveRDS(vcite_brts, "flat files/virus citation brts.rds")

rm(cite_brts, vcite_brts)
