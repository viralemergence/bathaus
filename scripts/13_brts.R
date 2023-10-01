# BRT analysis for bat viruses
# babetke@utexas.edu

# clear envrionement and graphics
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

# read in rds
# data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")

# # lab comp directory
data <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/cleaned dataset 30 cutoff.rds")

# how many NAs are there - 237
length(data$Synurbic[is.na(data$Synurbic)])

# zoonotic virus proportions
data <- data %>% 
  mutate(zoo_prop = zvirus/virus) %>%
  mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop))

# reservoir status
data <- data %>%
  mutate(dum_virus = if_else(virus <= 0, 0, 1),
         dum_zvirus = if_else(zvirus <= 0, 0, 1)) 

# transform zoonotic proportion to gaussian space here so it is included 
# in both datasets
c <- data[c("species","zoo_prop")]
c$xi <- c$zoo_prop * 100
c$ni <- 100

# Freeman tukey double arcsine transformation
t <- metafor::escalc(measure = "PFT", xi=xi, ni=ni, data=c)
# t$back <- metafor::transf.ipft(t$yi, t$ni) # back transformation

# subset 
tsub <- t[c("species","yi")]

# rename yi to zoo_pft
names(tsub) <- c("species","zoo_pft")

# merge to data by species
data <- merge(data, tsub, by = "species")

# remove
rm(c,t,tsub)

# remove duplicated diet variable and geographic realm
data <- data %>% 
  select(-c(det_inv, biogeographical_realm)) 

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
      select(-c("zvirus", "virus", "zoo_prop", "dum_zvirus", "dum_virus", "zoo_pft"))
    
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
    dist <- ifelse(response == "virus", "poisson", 
                   ifelse(response == "zvirus", "poisson", 
                          ifelse(response == "zoo_pft", "gaussian", "bernoulli")))
    
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
    bgrid <- makegrid(10, c(5000, 25000)) # binary models, tend to use much larger values
    
    # ### trim just for testing
    # pgrid <- pgrid[71, ]
    # pgrid$n.trees <- 5000
    # pgrid$shrinkage <- 0.001
    # pgrid$row <- 1
    
    ### run for all types of responses
    # Virus Richness and transformed zoonotic proportion
    vpars <- lapply(1:nrow(pgrid),function(x) grid_search(x, data_df = data, hgrid = pgrid, response="virus"))
    pftpars <- lapply(1:nrow(pgrid),function(x) grid_search(x, data_df = data, hgrid = pgrid, response="zoo_pft"))
    
    # trim just for testing
    # bgrid <- bgrid[84, ]
    
    # Reservoir status models
    vrespars <- lapply(1:nrow(bgrid),function(x) grid_search(x, data_df = data, hgrid = bgrid, response="dum_virus", cv = TRUE))
    zrespars <- lapply(1:nrow(bgrid),function(x) grid_search(x, data_df = data, hgrid = bgrid, response="dum_zvirus", cv = TRUE))
    
    # comment out for now because I already ran this
    ## get virus results
    vresults <- data.frame(sapply(vpars,function(x) x$trainRMSE),
                           sapply(vpars,function(x) x$testRMSE),
                           sapply(vpars,function(x) x$wrow),
                           sapply(vpars,function(x) x$best))
    names(vresults) <- c("trainRMSE","testRMSE","row","best")

    # Merge with hgrid
    vcomplete <- merge(vresults, pgrid, by = "row")
    
    # get zoonotic results
    zpresults <- data.frame(sapply(pftpars,function(x) x$trainRMSE),
                           sapply(pftpars,function(x) x$testRMSE),
                           sapply(pftpars,function(x) x$wrow),
                           sapply(pftpars,function(x) x$best))
    names(zpresults) <- c("trainRMSE","testRMSE","row","best")
    
    # Merge with hgrid
    zpcomplete <- merge(zpresults, pgrid, by = "row")
    
    # Overall virus reservoir models not sure if you want a separate grid for this b/c of what you ran b4
    ## get virus reservoir model results
    vrresults <- data.frame(sapply(vrespars,function(x) x$trainAUC),
                           sapply(vrespars,function(x) x$testAUC),
                           sapply(vrespars,function(x) x$spec),
                           sapply(vrespars,function(x) x$sen),
                           sapply(vrespars,function(x) x$wrow),
                           sapply(vrespars,function(x) x$best))
    names(vrresults) <- c("trainAUC","testAUC","spec","sen","row","best")
    
    # Merge with hgrid
    vrcomplete <- merge(vrresults, bgrid, by = "row")
    
    # zoonotic virus reservoir models
    ## get virus results
    zrresults <- data.frame(sapply(zrespars,function(x) x$trainAUC),
                           sapply(zrespars,function(x) x$testAUC),
                           sapply(zrespars,function(x) x$spec),
                           sapply(zrespars,function(x) x$sen),
                           sapply(zrespars,function(x) x$wrow),
                           sapply(zrespars,function(x) x$best))
    names(zrresults) <- c("trainAUC","testAUC","spec","sen","row","best")
    
    # Merge with hgrid
    zrcomplete <- merge(zrresults, bgrid, by = "row")
    
    # write as csv
    write_csv(vcomplete, "/Volumes/BETKE 2021/bathaus/flat files/virus grid search.csv")
    write_csv(zpcomplete, "/Volumes/BETKE 2021/bathaus/flat files/pft zoonotic grid search.csv")
    write_csv(vrcomplete, "/Volumes/BETKE 2021/bathaus/flat files/virus reservoir grid search.csv")
    write.csv(zrcomplete, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus reservoir grid search.csv")
    # may actually want to combine all of these for multiplots? or do it for like model outcomes.
    
} else {# read in grid search results

  # overall_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/virus grid search.csv")
  # zoo_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/pft zoonotic grid search.csv")
  # vres_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/virus reservoir grid search.csv")
  # zres_search <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/zoonotic virus reservoir grid search.csv")
  # 
  # lab comp paths
  overall_search <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/virus grid search.csv")
  zoo_search <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/pft zoonotic grid search.csv")
  vres_search <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/virus reservoir grid search.csv")
  zres_search <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus reservoir grid search.csv")
  
}

### plots of parameters - figure S1
library(patchwork)

#r <- overall_search %>% filter(shrinkage < 0.1)

# plots for richness models
virus_gg <- ggplot(overall_search, aes(x = factor(shrinkage), y = testRMSE)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = "Test RMSE", fill = "Interaction Depth", title = "Virus Richness") +
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
  labs(x = "Learning Rate", y = NULL, fill = "Interaction Depth", title = "Arcsine Zoonotic Proportion") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

# overall virus reservoir
vresAUC_gg <- ggplot(vres_search, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = "AUC", fill = "Interaction Depth", title = "Virus Reservoir") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

vresSen_gg <- ggplot(vres_search, aes(x = factor(shrinkage), y = sen)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = "Sensitiviy", fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

vresSpec_gg <- ggplot(vres_search, aes(x = factor(shrinkage), y = spec)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning rate", y = "Specificity", fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

# zoonotic reservoir
zresAUC_gg <- ggplot(zres_search, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Interaction Depth", title = "Zoonotic Reservoir") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

zresSen_gg <- ggplot(zres_search, aes(x = factor(shrinkage), y = sen)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

zresSpec_gg <- ggplot(zres_search, aes(x = factor(shrinkage), y = spec)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning rate", y = NULL, fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_fill_brewer(palette="Accent")

png("/Volumes/BETKE 2021/bathaus/figs/figure S1.png", width=5, height=7,units="in",res=300)
guide_area () + 
(virus_gg + zvirus_gg) / (vresAUC_gg + zresAUC_gg) / 
  (vresSen_gg + zresSen_gg) / (vresSpec_gg + zresSpec_gg) +
  plot_layout(guides = "collect", heights = c(1, 15)) & theme(legend.position = "top", legend.text = element_text(size = 8), legend.title = element_text(size = 8))
dev.off()

# remove plots from envrionment
rm(virus_gg, zvirus_gg, vresAUC_gg, vresSen_gg, vresAUC_gg, vresSpec_gg, zresAUC_gg, zresSen_gg, zresSpec_gg)

# unload patchwork
detach("package:patchwork", unload = TRUE)

# Sort output to view top model combinations (you want the lowest value for rmse)
sort <- overall_search %>% 
  arrange(testRMSE)

# look at summary of ids across seeds, arranged by lowest testRMSE
overall_search %>% 
  group_by(id) %>% 
  summarise(med = median(testRMSE)) %>% 
  arrange(med)

# Sort output to view top model combinations
psort <- zoo_search %>% 
  arrange(testRMSE)

# look at summary of ids across seeds, arranged by lowest testRMSE
zoo_search %>% 
  group_by(id) %>% 
  summarise(med = median(testRMSE)) %>% 
  arrange(med)

# sort
vrsort <- vres_search %>% 
  arrange(desc(testAUC))

# summarize
vres_search %>% 
  group_by(id) %>% 
  summarise(medAUC = median(testAUC),
            medSen = median(sen),
            medSpec = median(spec)) %>% 
  arrange(desc(medAUC))

# sort
zrsort <- zres_search %>% 
  arrange(desc(testAUC))

# summarize
zres_search %>% 
  group_by(id) %>% 
  summarise(medAUC = median(testAUC),
            medSen = median(sen),
            medSpec = median(spec)) %>% 
  arrange(desc(medAUC))

# remove sorts
rm(sort, psort, vrsort, zrsort)

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
      select(-c("virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop", "zoo_pft"))
    
  }else{ # remove synurbic from data
    # correct the response and remove raw response variables
    ndata <- ndata %>% 
      mutate(response = !!sym(response)) %>%
      select(-c("virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop", "zoo_pft", "Synurbic"))
    
    #fdata$Synurbic <- NULL # make null in fdata for prediction? Not sure if this is necessary if the formula dictates variables
  }
  
  # ifelse to determine distribution. If not any of the below options, throw error
  dist <- ifelse((response == "virus" | response == "zvirus" | response == "cites" | response == "vcites"), "poisson", 
                 ifelse((response == "dum_zvirus" | response == "dum_virus"),"bernoulli", 
                        ifelse(response == "zoo_pft","gaussian","Error")))
  
  # remove citations if predicting citations
  if(response == "cites" | response == "vcites") {
    ndata$cites = NULL
    ndata$vcites = NULL
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
    
    ## known - isn't this the same as ytest?
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
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","virus", "zvirus", "dum_zvirus", "dum_virus","zoo_prop","Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # then mean cites
    pdata <- fdata
    pdata$cites <- mean(pdata$cites)
    pdata$vcites <- mean(pdata$vcites) #just incase we also mean to do this as well
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
    
    # calculate train pseudo R squared
    r2train <- (1-(sum((yTrain - predict(gbmOut, newdata=train, n.trees=best.iter, type='response'))^2)/
                     sum((yTrain - mean(yTrain))^2)))
    
    # calculate test pseudo R squared
    r2test <- (1-(sum((yTest - predict(gbmOut, newdata=test, n.trees=best.iter, type='response'))^2)/
                    sum((yTest - mean(yTest))^2)))
    
    ## predict with cites
    preds <- predict(gbmOut,fdata,n.trees=best.iter,type="response")
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","virus","zvirus","zoo_prop","zoo_pft","Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # then mean cites
    pdata <- fdata
    pdata$cites <- mean(pdata$cites)
    pdata$vcites <- mean(pdata$vcites) #just incase we also mean to do this as well
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
      
      # backtransform cpred? not sure but I should probably ask dan
      pred_data$cback <- metafor::transf.ipft(pred_data$cpred, 100)
      
      # sort
      pred_data <- pred_data[order(pred_data$back,decreasing=T),]
      
    }
    
    ## print
    print(paste("BRT ",seed," done; test R2 = ", r2test, " for ", response, sep=""))
    
    ## save outputs
    return(list(mod = gbmOut,
                best = best.iter,
                testr2 = r2test,
                trainr2 = r2train,
                rinf = bars,
                predict = pred_data,
                traindata = train,
                testdata = test,
                seed = seed))
    
  }
  
}


#### apply across specified number of splits smax
# Filtered dataset
smax=100 # 1 just for testing on personal comp

# Richness models with and without synurbic
vrichness_brts <- lapply(1:smax,function(x) brts(seed = x,response = "virus", nt = 15000, shr = 0.001, int.d = 4, syn = "yes", cv = NULL))
no_vrichness_brts <- lapply(1:smax,function(x) brts( seed = x,response = "virus", nt = 15000, shr = 0.001, int.d = 4, syn = "no", cv = NULL))

# save then remove
saveRDS(vrichness_brts,"/Volumes/BETKE 2021/bathaus/flat files/virus with brts.rds")
saveRDS(no_vrichness_brts,"/Volumes/BETKE 2021/bathaus/flat files/virus without brts.rds")

rm(vrichness_brts, no_vrichness_brts)

# Proportion with and without
zoo_pft_brts <- lapply(1:smax,function(x) brts(seed = x,response = "zoo_pft", nt = 15000, shr = 0.001, int.d = 4, syn = "yes", cv = NULL))
no_zoo_pft_brts <- lapply(1:smax,function(x) brts(seed = x,response = "zoo_pft", nt = 15000, shr = 0.001, int.d = 4, syn = "no", cv = NULL))

saveRDS(zoo_pft_brts, "/Volumes/BETKE 2021/bathaus/flat files/zoo_prop with brts.rds")
saveRDS(no_zoo_pft_brts, "/Volumes/BETKE 2021/bathaus/flat files/zoo_prop without brts.rds")

rm(zoo_pft_brts,no_zoo_pft_brts)

# Overall virus reservoir status
vbinary_brts <- lapply(1:smax,function(x) brts(seed = x,response = "dum_virus", nt = 5000, shr = 0.01, int.d = 3, syn = "yes", cv = NULL))
no_vbinary_brts <- lapply(1:smax,function(x) brts(seed = x,response ="dum_virus", nt = 5000, shr = 0.01, int.d = 3, syn = "no", cv = NULL))

saveRDS(vbinary_brts, "/Volumes/BETKE 2021/bathaus/flat files/dum_virus with brts.rds")
saveRDS(no_vbinary_brts, "/Volumes/BETKE 2021/bathaus/flat files/dum_virus without brts.rds")

rm(vbinary_brts,no_vbinary_brts)

# Zoonotic virus reservoir
zbinary_brts <- lapply(1:smax,function(x) brts(seed = x,response = "dum_zvirus",  nt = 15000, shr = 0.0005, int.d = 3, syn = "yes", cv = NULL))
no_zbinary_brts <- lapply(1:smax,function(x) brts(seed = x,response ="dum_zvirus",  nt = 15000, shr = 0.0005, int.d = 3, syn = "no", cv = NULL))

saveRDS(zbinary_brts, "/Volumes/BETKE 2021/bathaus/flat files/dum_zvirus with brts.rds")
saveRDS(no_zbinary_brts, "/Volumes/BETKE 2021/bathaus/flat files/dum_zvirus without brts.rds")

rm(zbinary_brts,no_zbinary_brts)

# add functions for predicting citations (maybe I don't need to do the with/without data situation for these)
smax=20 # fewer iterations for citation models
cite_brts <- lapply(1:smax,function(x) brts(seed = x,response = "cites", nt = 15000, shr = 0.001, int.d = 4, syn = "yes", cv = NULL))
vcite_brts <- lapply(1:smax,function(x) brts(seed = x,response = "vcites", nt = 15000, shr = 0.001, int.d = 4, syn = "yes", cv = NULL))

saveRDS(cite_brts, "/Volumes/BETKE 2021/bathaus/flat files/citation brts.rds")
saveRDS(vcite_brts, "/Volumes/BETKE 2021/bathaus/flat files/virus citation brts.rds")
