# BRT analysis for bat viruses
# babetke@utexas.edu

# clear envrionement and graphics
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
# library(fastDummies) # dummy coding family
library(gbm)
library(rsample)
library(ROCR)
library(caret) 
library(InformationValue)

# read in rds
data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")

# # lab comp directory
# data <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/cleaned dataset 30 cutoff.rds")

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

# remove duplicated diet variable and geographic realm
data <- data %>% 
  select(-c(det_inv, biogeographical_realm)) 

# save before removing species
fdata <- data

# remove species
data$species <- NULL

# run tuning grid
# ifelse for running grid search function
gsrun = "yes" 
if(gsrun == "yes"){# run grid search 

  # Set up BRT tuning via search grid
  # function to make grids?
  ## hyperparameter grid, maybe allow the number of seeds to change?
  makegrid <- function(seed,trees) {
    
    # create grid
    tgrid <- expand.grid(n.trees = trees,
                         interaction.depth = c(2, 3, 4),
                         shrinkage = c(0.1, 0.01, 0.001, 0.0005),
                         #n.minobsinnode = c(5, 10, 15), 
                         #bag.fraction = c(0.5, 0.8, 1),
                         seed = seq(1,seed,by = 1))
    
    # adjust combinations with low n.trees and high shrinkage
    tgrid <- mutate(tgrid, n.trees = ifelse(shrinkage < 0.01 & n.trees == 5000, 15000, n.trees))
    
    ## trees, depth, shrink, min, prop
    tgrid$id <- with(tgrid,paste(n.trees,interaction.depth,shrinkage))
    
    ## sort by id then seed
    tgrid <- tgrid[order(tgrid$id,tgrid$seed),]
    
    return(tgrid)
  }
  
  #### Function to assess each hyperparameter combination
  # this function takes a search grid and runs them through gbms for a given split of data
  # then assess performance of each combination of parameters to get the optimal parameters for
  # final gbms.

  # row <- row in the hgrid
  # set <- the subset of data I want to use (two full datasets and two trimmed to complete)
  # response <- indicate which response I want to use. Either virus or zoo_prop
  # folds <- indicating either 10 or 5 folds depending on the dataset (5 for subset of 544 and 10 for full)
  grid_search <- function(row, data_df, response, folds, cv = NULL){
    
    # new data
    ndata <- data_df
    
    # correct the response and remove raw response variables
    ndata <- ndata %>% 
      mutate(response = !!sym(response)) %>%
      select(-c("zvirus", "virus", "zoo_prop", "dum_zvirus", "dum_virus"))
    
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
                          ifelse(response == "zoo_prop_arcine", "gaussian", "bernoulli")))
    
    ## BRT
    set.seed(1)
    gbmOut=gbm(response ~ . ,data=dataTrain,
               n.trees=hgrid$n.trees[row],
               distribution=dist,
               shrinkage=hgrid$shrinkage[row],
               interaction.depth=hgrid$interaction.depth[row],
               n.minobsinnode=10,
               cv.folds=folds,
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
      print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))
      
      ## save outputs
      return(list(best=best.iter,
                  trainAUC=auc_train,
                  testAUC=auc_test,
                  spec=spec,
                  sen=sen,
                  wrow=row))
      
    } else {
      
      # calculate RMSE
      mse <- mean((preds - result)^2)
      rmse <- sqrt(mse)
      
      # print
      print(paste("hpar row ",row," done; RMSE is ",rmse,sep=""))
      
      # save outputs
      return(list(best = best.iter,
                  mse = mse,
                  rmse = rmse,
                  wrow = row))
    }
    
  }
  
    # hgrid <- makegrid(1, c(10000, 20000))
    # vgrid <- makegrid(1,c(500, ))
    hgrid <- makegrid(10, 5000) # just use 5000 and 15000, best iteration for 25000 on very small int. depth doesn't max out
    
    # # Removing combinations that constantly max out
    # hgrid %>% 
    #   filter(!(n.trees == 5000 & shrinkage < 0.01)) %>%
    #   filter(!(n.trees == 10000 & shrinkage < 0.001)) -> hgrid
    
    # renumber the rows for matching 
    hgrid$row <- 1:nrow(hgrid)
    hgrid$id2=factor(as.numeric(factor(hgrid$id)))
    
    ### trim for testing
    # hgrid <- hgrid[48, ]
    # hgrid$n.trees <- 15000
    # hgrid$shrinkage <- 0.0005
    
    # run for the two types of data?
    vpars <- lapply(1:nrow(hgrid),function(x) grid_search(x, data_df = data, response="virus", folds = 10))
    zpars <- lapply(1:nrow(hgrid),function(x) grid_search(x, data_df = data, response="zvirus", folds = 10))
    
    ## get virus results
    vresults <- data.frame(sapply(vpars,function(x) x$mse),
                           sapply(vpars,function(x) x$rmse),
                           sapply(vpars,function(x) x$wrow),
                           sapply(vpars,function(x) x$best))
    names(vresults) <- c("MSE","RMSE","row","best")
    
    # Merge with hgrid
    vcomplete <- merge(vresults, hgrid, by = "row")
    
    # get zoonotic results
    zresults <- data.frame(sapply(zpars,function(x) x$mse),
                           sapply(zpars,function(x) x$rmse),
                           sapply(zpars,function(x) x$wrow),
                           sapply(zpars,function(x) x$best))
    names(zresults) <- c("MSE","RMSE","row","best")
    
    # Merge with hgrid
    zcomplete <- merge(zresults, hgrid, by = "row")
    
    # write as csv
    write_csv(vcomplete, "/Volumes/BETKE 2021/bathaus/flat files/virus grid search.csv")
    write_csv(zcomplete, "/Volumes/BETKE 2021/bathaus/flat files/zvirus grid search.csv")
    
} else {# read in grid search results
  
  overall_search <- read_csv("~/Desktop/Bats and Viurses/bathaus/flat files/virus grid search.csv")
  zoo_search <- read_csv("~/Desktop/Bats and Viurses/bathaus/flat files/zvirus grid search.csv")
  
}

### plots of parameters - these are from synurbat. need to update to this proj.
auc_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = "AUC", fill = "Interaction Depth", title = "Initial Model") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

tree_gg <- ggplot(na.complete, aes(x = factor(n.trees), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "No.Trees", y = "AUC", fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

pauc_gg <- ggplot(p.complete, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = NULL, fill = "Interaction Depth", title = "Pseudoabsence Model") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

ptree_gg <- ggplot(p.complete, aes(x = factor(n.trees), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "No.Trees", y = NULL, fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

library(patchwork) # multiplot and save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure S4.png", width=6, height=6,units="in",res=300)
guide_area () + 
  (auc_gg + pauc_gg) / (tree_gg + ptree_gg) + 
  plot_layout(guides = "collect", heights = c(1, 15)) & theme(legend.position = "top")
dev.off()

# remove
rm(auc_gg, tree_gg, pauc_gg, ptree_gg)

# unload patchwork
detach("package:patchwork", unload = TRUE)

# Sort output to view top model combinations
sort <- na.complete %>% 
  arrange(desc(testAUC))

# Sort output to view top model combinations
psort <- p.complete %>% 
  arrange(desc(testAUC))

# remove sorts
rm(sort, psort)

#### Define brt function 
# modify this to make it appropriate for the bat data
# Original takes:
# seed - the seeds to be used for the unique splits, this also dictates the number of splits
# response - indicating what the response variable is

# what I am adding:
# data - control the dataset I plan to use because I have full and complete versions
# parameter_df - subsets for that response and then filters by best AUC, use it to pull values to ref for brt
# full - cv folds differ between full and filtered - maybe I could just do a ifelse statement that takes yes/no
brts <- function(data_df, seed, response, nt, shr, int.d, cv = NULL){
  
  # new data
  ndata <- data_df
  
  # correct the response and remove raw response variables
  ndata <- ndata %>% 
    mutate(response = !!sym(response)) %>%
    select(-c("virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop"))
  
  # # parameter dataset filter by response and select the best parameters by AUC, resulting in a single row of data
  # # reference it for the gbm parameters
  # params <- parameter_df %>%
  #   filter(type == response) %>%
  #   filter(testAUC == max(testAUC))
  
  # ifelse
  dist <- ifelse(response == "virus", "poisson", 
                 ifelse(response == "zoo_prop", "gaussian", "bernoulli"))
  
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
  
  ## known - isnt this the same as test?
  result <- test$response
  
  # output dependent on dist
  # AUC, sen, spec, ROC for binary models
  if(dist == "bernoulli"){ 
  
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
  pred_data <- fdata[c("species","dum_virus","dum_zvirus","virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop","Synurbic")]
  pred_data$pred <- preds
  pred_data$type <- response
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ", auc_test,sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=train,
              testdata=test,
              seed=seed))
  
  }else{ # pseudo R2 for poisson and gaussian models
    
    ## relative importance
    bars <- summary(gbmOut,n.trees=best.iter,plotit=F)
    bars$rel.inf <- round(bars$rel.inf,2)
    
    if(dist == "poisson"){
    ## predict with cites
    preds <- predict(gbmOut,fdata,n.trees=best.iter,type="response")
    pred_data <- fdata[c("species","dum_virus","dum_zvirus","virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop","Synurbic")]
    pred_data$pred <- preds
    pred_data$type <- response
    
    # ## predict with mean cites
    # pdata <- data_df
    # pdata$cites <- mean(pdata$cites)
    # pred_data$cpred <- predict(gbmOut,pdata,n.trees=best.iter,type="response")
    
    } else {
      
      ## predict with cites
      preds <- predict(gbmOut,fdata,n.trees=best.iter,type="response")
      pred_data <- fdata[c("species","dum_virus","dum_zvirus","virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop","Synurbic")]
      pred_data$pred <- preds
      pred_data$type <- response
      
      # back transforming predictions
      pred_data$ni <- 100
      pred_data$back <- metafor::transf.ipft(pred_data$pred, pred_data$ni)
      
    }
    
    ## sort
    pred_data <- pred_data[order(pred_data$pred,decreasing=T),]
    
    # calculate train pseudo R squared
    r2train <- (1-(sum((yTrain - predict(gbmOut, newdata=dataTrain, n.trees=best.iter, type='response'))^2)/
                    sum((yTrain - mean(yTrain))^2)))
    
    # calculate test pseudo R squared
    r2test <- (1-(sum((yTest - predict(gbmOut, newdata=dataTest, n.trees=best.iter, type='response'))^2)/
                   sum((yTest - mean(yTest))^2)))
    
    ## print
    print(paste("BRT ",seed," done; test R2 = ", r2test, sep=""))
    
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
smax=100

# Richness models with and without synurbic
vrichness_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response = "virus", parameter_df = search_full, full = "yes"))
no_vrichness_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response="virus", parameter_df = search_full, full = "yes"))

# Proportion with and without
zoo_prop_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response = "zoo_prop", parameter_df = search_full, full = "yes"))
no_zoo_prop_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response="zoo_prop", parameter_df = search_full, full = "yes"))

# Overall virus reservoir status
vbinary_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response = "dum_virus", parameter_df = search_full, full = "yes"))
no_vbinary_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response="dum_zvirus", parameter_df = search_full, full = "yes"))

# Zoonotic virus reservoir
zbinary_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response = "dum_virus", parameter_df = search_full, full = "yes"))
no_zbinary_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response="dum_zvirus", parameter_df = search_full, full = "yes"))

# write to files
#setwd("~/Desktop/hantaro/data/clean files")
saveRDS(virus_brts,"virus brts.rds")
saveRDS(zvirus_brts,"zvirus brts.rds")



####################### TEST CODE
# new data
ndata <- data

# correct the response and remove raw response variables
ndata <- ndata %>% 
  mutate(response = virus) %>%
  select(-c("zvirus", "virus", "zoo_prop"))

## use rsample to split (allow the proportion to be changed incase)
#set.seed(hgrid$seed[row])
split=initial_split(ndata,prop=0.8,strata="response")

## test and train
dataTrain=training(split)
dataTest=testing(split)

## yTest and yTrain
yTrain=dataTrain$response
yTest=dataTest$response

## BRT
set.seed(1)
gbmOut=gbm(response ~ . ,data=dataTrain,
           n.trees=15000,
           distribution="poisson",
           shrinkage=0.0001,
           interaction.depth=4,
           n.minobsinnode=5,
           cv.folds=10,
           bag.fraction=0.5,
           train.fraction=1,
           n.cores=1,
           verbose=F)

## performance
par(mfrow=c(1,1),mar=c(4,4,1,1))
best.iter=gbm.perf(gbmOut,method="cv")

## predict with test data
preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")

## known
result=dataTest$response

# ## sensitiviy and specificity
# sen=InformationValue::sensitivity(result,preds)
# spec=InformationValue::specificity(result,preds)
# 
# ## AUC on train
# auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
# 
# ## AUC on test
# auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
# 
# ## print
# print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))

# calculate RMSE
mse = mean((preds - result)^2)
rmse = sqrt(mse)


#Perfomring caculation to calculate th pseduo R2 of the training dataset
R2.train<-(1-(sum((train.response - predict(gbm.obj, newdata=data, n.trees=best.ntrees, type='response'))^2)/
                sum((train.response - mean(train.response))^2)))

#Performing calculation to calcualted the pseduo R2 of the test dataset
R2.test<-(1-(sum((test.response - predict(gbm.obj, newdata=data.test, n.trees=best.ntrees, type='response'))^2)/
               sum((test.response - mean(test.response))^2)))



# zoonotic virus proportions
data <- data %>% 
  mutate(zoo_prop = zvirus/virus) %>%
  mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop))

# new data
ndata <- data

# correct the response and remove raw response variables
ndata <- ndata %>% 
  mutate(response = zoo_prop) %>%
  select(-c("zvirus", "virus", "zoo_prop"))

## use rsample to split (allow the proportion to be changed incase)
set.seed(345)
split=initial_split(ndata,prop=0.8,strata="response")

## test and train
dataTrain=training(split)
dataTest=testing(split)

## yTest and yTrain
yTrain=dataTrain$response
yTest=dataTest$response

## BRT
set.seed(1)
gbmOut=gbm(response ~ . ,data=dataTrain,
           n.trees=5000,
           distribution="gaussian",
           shrinkage=0.001,
           interaction.depth=4,
           n.minobsinnode=5,
           cv.folds=10,
           bag.fraction=0.5,
           train.fraction=1,
           n.cores=1,
           verbose=F)

## performance
par(mfrow=c(1,1),mar=c(4,4,1,1))
best.iter=gbm.perf(gbmOut,method="cv")

## predict with test data
preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")

## known
result=dataTest$response

# calculate RMSE
mse = mean((preds - result)^2)
rmse = sqrt(mse)





tgrid <- makegrid(10, c(5000, 20000))

tgrid$n.trees=ifelse(tgrid$shrinkage<0.001,tgrid$n.trees*3,tgrid$n.trees)




seq(500, 1000, 50)

tgrid <- expand.grid(n.trees = seq(500, 1000, 250),
                     interaction.depth = c(2, 3, 4),
                     shrinkage = c(0.1, 0.01, 0.001))
                     #n.minobsinnode = 10) 
                     #bag.fraction = c(0.5, 0.8, 1)


library(car)
c <- data
c$logit_boot <- boot::logit(c$zoo_prop)
c$logitBack <- inv.logit(c$logitzoo)
c$logitBackPlus <- c$logitBack + 0.025
c$logitStack <- inv.logit(c$logitzoo,a=0.025)*100

# stack exchange code - https://stackoverflow.com/questions/23845283/logit-transformation-backwards
inv.logit <- function(f,a) {
  a <- (1-2*a)
 (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
  }

format(c$logitStack, scientific = FALSE)

asin(sqrt(c$zoo_prop))


##### testing transformations with metafor package
install.packages("metafor")
library(metafor) #no mention of masking but some problems with %>%. may want to use metafor::escalc()

c <- fdata[c("species","zoo_prop")]
c$xi <- c$zoo_prop * 100
c$ni <- 100

g <- metafor::escalc(measure = "PLO", xi=xi, ni=ni, data=c)
g$back <- metafor::transf.ilogit(g$yi)

t <- metafor::escalc(measure = "PFT", xi=xi, ni=ni, data=c)
t$back <- metafor::transf.ipft(t$yi, t$ni)

# subset 
g[c("species","yi")]
tsub <- t[c("species","yi")]

zdata <- merge(fdata, tsub, by = "species")

########
# new data
ndata <- zdata

# correct the response and remove raw response variables
ndata <- ndata %>% 
  mutate(response = yi) %>%
  select(-c("virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop", "yi", "species"))

## use rsample to split (allow the proportion to be changed incase)
#set.seed(hgrid$seed[row])
split=initial_split(ndata,prop=0.8,strata="response")

## test and train
dataTrain=training(split)
dataTest=testing(split)

## yTest and yTrain
yTrain=dataTrain$response
yTest=dataTest$response

## BRT
set.seed(1)
gbmOut=gbm(response ~ . ,data=dataTrain,
           n.trees=500,
           distribution="poisson",
           shrinkage=0.01,
           interaction.depth=4,
           n.minobsinnode=5,
           cv.folds=10,
           bag.fraction=0.5,
           train.fraction=1,
           n.cores=1,
           verbose=F)

## performance
par(mfrow=c(1,1),mar=c(4,4,1,1))
best.iter=gbm.perf(gbmOut,method="cv")

## predict with test data
preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")

## known
result=dataTest$response

## relative importance
bars <- summary(gbmOut,n.trees=best.iter,plotit=F)
bars$rel.inf <- round(bars$rel.inf,2)

## predict with cites
preds <- predict(gbmOut,fdata,n.trees=best.iter,type="response")
pred_data <- fdata[c("species","dum_virus","dum_zvirus","virus", "zvirus", "dum_zvirus", "dum_virus", "zoo_prop","Synurbic")]
pred_data$pred <- preds
pred_data$type <- "zoo_prop"

# ## predict with mean cites
# pdata <- data_df
# pdata$cites <- mean(pdata$cites)
# pred_data$cpred <- predict(gbmOut,pdata,n.trees=best.iter,type="response")

## sort
pred_data <- pred_data[order(pred_data$pred,decreasing=T),]

# back transforming
pred_data$ni <- 100
pred_data$back <- transf.ipft(pred_data$pred, pred_data$ni)

# calculate train pseudo R squared
r2train <- (1-(sum((yTrain - predict(gbmOut, newdata=dataTrain, n.trees=best.iter, type='response'))^2)/
                 sum((yTrain - mean(yTrain))^2)))

# calculate test pseudo R squared
r2test <- (1-(sum((yTest - predict(gbmOut, newdata=dataTest, n.trees=best.iter, type='response'))^2)/
                sum((yTest - mean(yTest))^2)))

