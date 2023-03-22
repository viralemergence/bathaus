# 10_data cleaning - processing data for use in brt models
# Adapted cleaning code from hantaro brt

# Script for data cleaning
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
library(fastDummies) # dummy coding family
library(gbm)
library(rsample)
library(ROCR)
library(caret) 
library(InformationValue)

# Read in the merged dataset
setwd("/Users/brianabetke/Desktop/Bats and Viruses/bathaus")
data <- read_csv("flat files/zoonotic virus to merged PanTHERIA.csv")

# glimpse and see variables
glimpse(data)

## make binary columns for genus
dums <- dummy_cols(data["fam"])

## unique, does not appear to have duplicates?
dums <- dums[!duplicated(dums$fam),]

# factor all variables
# the for loop turned everything to NA and gave 19 warnings
dums <- dums %>% 
  mutate(across(where(is.numeric), factor))

## merge - variables dont stay as factors when you merge
data <- merge(data,dums,by="fam",all.x=T)
rm(dums)

# Factor COMBINE and remove unnecessary columns
data <- data %>% # Synurbic and variables that are factors according to COMBINE
  mutate(across(c("Synurbic","hibernation_torpor","fossoriality","trophic_level",
                  "foraging_stratum","activity_cycle", "freshwater", 
                  "marine","terrestrial_non-volant", "terrestrial_volant","island_endemicity",
                  "disected_by_mountains", "glaciation", "biogeographical_realm"), 
                factor)) %>% 
  select(-c("MSW3_sciName_matched","vfilter","filter"))

colnames(data) # 96 columns

# Variation
## mode function
mode.prop <- function(x) {
  ux <- unique(x[is.na(x)==FALSE])
  tab <- tabulate(match(na.omit(x), ux))
  max(tab)/length(x[is.na(x)==FALSE])
}

ux <- unique(data$weaning_mass_g[is.na(data$weaning_mass_g)==FALSE])
tab <- tabulate(match(na.omit(data$weaning_mass_g), ux))
max(tab)/length(data$weaning_mass_g[is.na(data$weaning_mass_g)==FALSE])

## assess variation across columns
vars=data.frame(apply(data,2,function(x) mode.prop(x)),
                apply(data,2,function(x) length(unique(x))))

## get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")
vars$var=round(vars$var,2)

## if homogenous (100%)
vars$keep=ifelse(vars$var<1,"keep","cut")
#vars$keep=ifelse(vars$column%in%c('hPCR','competence',"fam"),'keep',vars$keep)
vars=vars[order(vars$keep),]

table(vars$keep)
# cut keep 
# 17   77 

## trim
keeps=vars[-which(vars$keep=="cut"),]$column

## drop if no variation
data=data[keeps]
rm(keeps,vars)

## assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(data)))

## get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")
mval$comp=round(mval$comp,2)

# ggplot of coverage
setwd("/Users/brianabetke/Desktop/Bats and Viruses")
png("ESA 2022/figs/trait_coverage.png", width=9.5,height=5.5,units="in",res=600)
mval %>% 
  filter(!column %in% c("fam","cnames", "tip", "gen", "cites", "vcites","Complete")) %>%
  ggplot(aes(comp)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept=0.30,linetype=2,size=0.5) +
  #geom_vline(xintercept=0.25,linetype=2,size=0.5) +
  theme_bw() +
  labs(y="frequency",
      x="trait coverage across bat species") +
  scale_x_continuous(labels = scales::percent)
dev.off()
  
mval$keep=ifelse(mval$comp>=0.30,"keep","cut")
table(mval$keep)
# cut keep 
# 18   63

## order
mval=mval[order(mval$comp),]
keeps=mval[-which(mval$keep=="cut"),]$column

# cleaned up table of keeps
coverage_table <- mval %>% 
  filter(keep == "keep") %>%
  rename(Variable = column,
         Coverage = comp) %>%
  filter(!Variable %in% c("clade", "gen", "tip", "species", "fam", "virus", "zvirus", "Complete")) %>%
  select(-keep) %>%
  relocate(Coverage, .after = Variable)

rownames(coverage_table) <- NULL

# save as csv
write.csv(coverage_table, "ESA Poster/figs/coverage_table.csv", row.names = FALSE)

# Write to .txt
write.table(coverage_table, file = "ESA Poster/figs/coverage_table.txt", sep = ",", quote = FALSE, row.names = F)

## drop if not well represented
data=data[keeps]
rm(mval,keeps,coverage_table)

## Clean out remaining variables?
data <- data %>% 
  select(-c("tip","gen","fam","clade"))

colnames(data) # resulting in 58 variables total, 56 covariates

# create binary variable for viruses and zviruses
data <- data %>%
  mutate(dum_virus = if_else(virus <= 0, 0, 1),
         dum_zvirus = if_else(zvirus <= 0, 0, 1)) 
glimpse(data)

# Without species and complete 
full <- data %>%
  select(-c("species", "Complete"))

# Write as csv
write.csv(full, "flat files/cleaned dataset.csv", row.names = FALSE)

# virus data trimmed by complete
filtered <- data %>%
  filter(Complete == "yes") %>%
  select(-c("species", "Complete"))

# you need to check the variation of variables again after subsetting because you got like
# 50 errors
vars=data.frame(apply(filtered,2,function(x) mode.prop(x)),
                apply(filtered,2,function(x) length(unique(x))))

## get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")
vars$var=round(vars$var,2) # bumps three additional variables to 1 

# which ones are = 1?
remove <- vars$column[vars$var == 1]
filtered <- filtered %>%
  select(-remove)
rm(remove)

# Write as csv
write.csv(filtered, "flat files/cleaned dataset filtered.csv", row.names = FALSE)

#### Adapted parameter tuning for all subsets of data
## hyperparameter grid
hgrid=expand.grid(n.trees=5000,
                  interaction.depth=c(2,3,4),
                  shrinkage=c(0.01,0.001,0.0005),
                  n.minobsinnode=4,
                  seed=50)

## fix trees
hgrid$n.trees=ifelse(hgrid$shrinkage<0.001,hgrid$n.trees*3,hgrid$n.trees)

## trees, depth, shrink, min, prop
hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))

## sort by id then seed
hgrid=hgrid[order(hgrid$id,hgrid$seed),]

## now add rows
hgrid$row=1:nrow(hgrid)

## factor id
hgrid$id2=factor(as.numeric(factor(hgrid$id)))

#### Function to assess each hyperparameter combination
# this function takes a search grid and runs them through gbms for a given split of data
# then assess performance of each combination of parameters to get the optimal parameters for
# final gbms.


# I added set and split_prop so could change datasets and proportion of splits
# row <- row in the hgrid (I assume)
# set <- the subset of data I want to use (two full datasets and two trimmed to complete)
# response <- indicate which response I want to use. Either dum_virus or dum_zvirus
# folds <- indicating either 10 or 5 folds depending on the dataset (5 for subset of 544 and 10 for full)
hfit=function(row, set, response, folds){
  
  ## make new data
  ndata=set
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata <- ndata %>% 
    select(-c("virus", "zvirus", "dum_zvirus", "dum_virus"))
  
  ## use rsample to split (allow the propotion to be changed incase)
  set.seed(hgrid$seed[row])
  split=initial_split(ndata,prop=0.7,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=hgrid$n.trees[row],
             distribution="bernoulli",
             shrinkage=hgrid$shrinkage[row],
             interaction.depth=hgrid$interaction.depth[row],
             n.minobsinnode=hgrid$n.minobsinnode[row],
             cv.folds=folds,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## print
  print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))
  
  ## save outputs
  return(list(best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              wrow=row))
}


vpars=lapply(1:nrow(hgrid),function(x) hfit(x, set = filtered, response="dum_virus", folds = 5))

## get results
vresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                    sapply(hpars,function(x) x$testAUC),
                    sapply(hpars,function(x) x$spec),
                    sapply(hpars,function(x) x$sen),
                    sapply(hpars,function(x) x$wrow),
                    sapply(hpars,function(x) x$best))
names(vresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
vsearch_complete = merge(vresults,hgrid,by="row")

## save
vsearch_complete$type="dum_virus"

# run for zvirus
zpars=lapply(1:nrow(hgrid),function(x) hfit(x, set = filtered, response="dum_zvirus", folds = 5))

## get results
zresults=data.frame(sapply(zpars,function(x) x$trainAUC),
                    sapply(zpars,function(x) x$testAUC),
                    sapply(zpars,function(x) x$spec),
                    sapply(zpars,function(x) x$sen),
                    sapply(zpars,function(x) x$wrow),
                    sapply(zpars,function(x) x$best))
names(zresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
zsearch_complete = merge(zresults,hgrid,by="row")

## save
zsearch_complete$type="dum_zvirus"

search=rbind.data.frame(vsearch_complete, zsearch_complete)
search$type=factor(search$type,levels=c("dum_virus","dum_zvirus"))
write.csv(search, "search_Complete.csv", row.names = FALSE)

## Do this for full datasets (set1 with 10 folds)
vpars_full=lapply(1:nrow(hgrid),function(x) hfit(x, set = full, response="dum_virus", folds = 10))

## get results
vresults_full=data.frame(sapply(vpars_full,function(x) x$trainAUC),
                    sapply(vpars_full,function(x) x$testAUC),
                    sapply(vpars_full,function(x) x$spec),
                    sapply(vpars_full,function(x) x$sen),
                    sapply(vpars_full,function(x) x$wrow),
                    sapply(vpars_full,function(x) x$best))
names(vresults_full)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
vsearch_full=merge(vresults,hgrid,by="row")

## save
vsearch_full$type="dum_virus"

## Full zvirus
zpars_full=lapply(1:nrow(hgrid),function(x) hfit(x, set = full, response="dum_zvirus", folds = 10))

## get results
zresults_full=data.frame(sapply(zpars_full,function(x) x$trainAUC),
                         sapply(zpars_full,function(x) x$testAUC),
                         sapply(zpars_full,function(x) x$spec),
                         sapply(zpars_full,function(x) x$sen),
                         sapply(zpars_full,function(x) x$wrow),
                         sapply(zpars_full,function(x) x$best))
names(zresults_full)=c("trainAUC","testAUC",
                       "spec","sen","row","best")

## combine and save
zsearch_full=merge(zresults_full,hgrid,by="row")

## save
zsearch_full$type="dum_zvirus"

# write csv
search=rbind.data.frame(vsearch_full, zsearch_full)
search$type=factor(search$type,levels=c("dum_virus","dum_zvirus"))
write.csv(search, "search_full.csv", row.names = FALSE)

#### Define brt function - at this point, I do not want to mess with the cites portion
# modify this to make it appropriate for the bat data
# Original takes:
# seed - the seeds to be used for the unique splits, this also dictates the number of splits
# response - indicating what the response variable is

# what I am adding:
# data - control the dataset I plan to use because I have full and complete versions
# parameter_df - subsets for that response and then filters by best AUC, use it to pull values to ref for brt
# full - cv folds differ between full and filtered - maybe I could just do a ifelse statement that takes yes/no

brts <- function(data_df, seed, response, parameter_df, full){
  
  # new data
  ndata <- data_df
  
  # correct the response and remove raw response variables
  ndata <- ndata %>% 
    mutate(response = !!sym(response)) %>%
    select(-c("virus", "zvirus", "dum_zvirus", "dum_virus"))
  
  # parameter dataset filter by response and select the best parameters by AUC, resulting in a single row of data
  # reference it for the gbm parameters
  params <- parameter_df %>%
    filter(type == response) %>%
    filter(testAUC == max(testAUC))
  
  ## use rsample to split
  set.seed(seed)
  split <- initial_split(ndata,prop=0.7,strata="response")
  
  ## test and train for data
  train <- training(split)
  test <- testing(split)
  
  # pull response test and train
  yTrain <- train$response
  yTest <- test$response
  
  # If else statement for cv.folds by dataset 
  folds <- ifelse(full=="yes", 10, 5)
  
  ## parameters from parameter_df - search grid data
  nt <- params$n.trees
  shr <- params$shrinkage
  int.d <- params$interaction.depth
  
  ## BRT
  set.seed(1)
  gbmOut <- gbm(response ~ . ,data=train,
                n.trees=nt,
                distribution="bernoulli",
                shrinkage=shr,
                interaction.depth=int.d,
                n.minobsinnode=4,
                cv.folds=folds,class.stratify.cv=TRUE,
                bag.fraction=0.5,train.fraction=1,
                n.cores=1,
                verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter <- gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds <- predict(gbmOut,test,n.trees=best.iter,type="response")
  
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
  ## I removed the portion for cites  but this portion seems important for when variance of ytest is zero?
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
  
  # ## predict with cites
  # preds <- predict(gbmOut,data,n.trees=best.iter,type="response")
  # pred_data <- data[c("tip",'treename',"fam","gen","hPCR","competence")]
  # pred_data$pred <- preds
  # pred_data$type <- response
  # 
  # ## predict with mean cites
  # pdata <- data_df
  # pdata$cites <- mean(pdata$cites)
  # pred_data$cpred <- predict(gbmOut,pdata,n.trees=best.iter,type="response")
  # 
  # ## sort
  # pred_data <- pred_data[order(pred_data$pred,decreasing=T),]
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ",auc_test,sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              roc=perf,
              rinf=bars,
              #predict=pred_data,
              traindata=train,
              testdata=test,
              seed=seed))
  
}

#### Running function on datasets
# hyper parameter grid data
search_filtered <- read.csv("/Volumes/BETKE 2021/bathaus/search_Complete.csv")
search_full <- read.csv("/Volumes/BETKE 2021/bathaus/search_full.csv")

#### apply across specified number of splits smax
# Filtered dataset
smax=100
fvirus_brts <- lapply(1:smax,function(x) brts(data_df = filtered, seed = x,response = "dum_virus", parameter_df = search_filtered, full = "no"))
fzvirus_brts <- lapply(1:smax,function(x) brts(data_df = filtered, seed=x,response= "dum_zvirus", parameter_df = search_filtered, full = "no"))

# write to files
#setwd("~/Desktop/hantaro/data/clean files")
saveRDS(fvirus_brts,"fvirus brts.rds")
saveRDS(fzvirus_brts,"fzvirus brts.rds")

# full dataset
virus_brts <- lapply(1:smax,function(x) brts(data_df = full, seed = x,response = "dum_virus", parameter_df = search_full, full = "yes"))
zvirus_brts <- lapply(1:smax,function(x) brts(data_df = full, seed=x,response="dum_zvirus", parameter_df = search_full, full = "yes"))

# write to files
#setwd("~/Desktop/hantaro/data/clean files")
saveRDS(virus_brts,"virus brts.rds")
saveRDS(zvirus_brts,"zvirus brts.rds")





































