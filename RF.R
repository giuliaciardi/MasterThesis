# Machine Learning with H2o
# Load the H2O library and start up the H2O cluster locally on your machine
library(h2o)
library(caret)
h2o.init(nthreads = -1, #Number of threads -1 means use all cores on your machine
         max_mem_size = "8G")  #max mem size is the maximum memory to allocate to H2O

h2o.flow()

#import data
data <- read.csv("C:/Users/giuli/Google Drive/TESI/MyData_ACVD.csv", stringsAsFactors = T)
as.data.frame(data)
data$Enrichment..0.control..1.case.<- as.factor(data$Enrichment..0.control..1.case.)
data <- data[,-c(1,2)]

# remove high correlated
trainDescr <- data[,-1]
descrCorr <- cor(trainDescr)
highCorr <- findCorrelation(descrCorr, 0.85)
trainDescr <- trainDescr[, -highCorr]
ncol(trainDescr) #287 non high correlated predittori rimasti
colnames(trainDescr)
data <- data[, -highCorr]

# h2o
data <- as.h2o(data)
dim(data)
h2o.levels(data$Enrichment..0.control..1.case.)

# Partition the data into training and test sets
splits <- h2o.splitFrame(data = data, 
                         ratios = 0.8 ,  #partition data into 80%, 20%
                         seed = 1)

#setting a seed will guarantee reproducibility
train <- splits[[1]]
test <- splits[[2]]

# Take a look at the size of each partition
nrow(train)  # 325
nrow(test)  # 80

# Identify response and predictor variables
y <- "Enrichment..0.control..1.case."
x <- setdiff(names(data), y)
nfolds <- 10


# H2O's Random Forest (RF) implements a distributed version of the standard 
# Random Forest algorithm and variable importance measures.
# First we will train a basic Random Forest model with default parameters. 
# A seed is required for reproducibility

rf_fit1 <- h2o.randomForest(x = x,
                            y = y,
                            training_frame = train,
                            nfolds = nfolds,
                            model_id = "rf_fit1",
                            seed = 1) 

rf_fit1@model$model_summary

#Cross Validation metrics
rf_fit1@model$cross_validation_metrics
rf_fit1@model$cross_validation_metrics_summary

rf_fit1@model$variable_importances[1:20,]

h2o.auc(rf_fit1, xval = T)
h2o.auc(rf_fit1, train = T)

#variable importance table
vi <- as.data.frame(h2o.varimp(rf_fit1))
vi20 <- vi[1:20,]; vi20

## Create predictions using our latest RF model against the test set.
finalRf_predictions<-h2o.predict(object = rf_fit1,
                                 newdata = test)
finalRf_predictions

finalRF_test <- h2o.performance(model = rf_fit1,
                                newdata = test)
finalRF_test@metrics$MSE

h2o.confusionMatrix(finalRF_test)
h2o.auc(finalRF_test)
h2o.sensitivity(finalRF_test, thresholds = 0.44)
h2o.specificity(finalRF_test, thresholds = 0.44)


## We could further experiment with deeper trees or a higher percentage of
##  columns used (mtries).

############## CARTESIAN GRID SEARCH #####################
# RF Hyperparamters
# ntrees_opt = c(50, 100, 200, 500, 1000)
# mtries_opt <- c(8, 10, 15,20)
# max_depth_opt <- c(5, 10, 15, 20, 25)
# sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
# col_sample_rate_per_tree_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
# 
# hyper_params <- list(ntrees = ntrees_opt,
#                      mtries = mtries_opt,
#                      max_depth = max_depth_opt)
#                      sample_rate = sample_rate_opt,
#                      col_sample_rate_per_tree = col_sample_rate_per_tree_opt)
# 
# rf_gridC <- h2o.grid("randomForest", 
#                      grid_id = "mygridC",
#                      x = x, 
#                      y = y,
#                     training_frame = train,
#                     seed = 123,
#                     nfolds = nfolds,
#                     keep_cross_validation_predictions = F,                    
#                     hyper_params = hyper_params)

#importo grid di Luigi
# load("~/Dropbox (Kiwi)/Kiwi Team Folder/SHARED/PROGETTI/Microbioma/TESI/my_grid.RData")
# 
# rf_sorted_grid <- h2o.getGrid(grid_id = "mygridC",sort_by = "accuracy", decreasing = T)
# print(rf_sorted_grid)

# best_model <- h2o.getModel(rf_sorted_grid@model_ids[[2]])
# summary(best_model)
# 
# best_rf_perf <- h2o.performance(model = best_model,
#                                   newdata = test)
# 
# h2o.confusionMatrix(best_rf_perf)
# h2o.auc(best_rf_perf)
# 
# Look at the hyperparameters for the best model
# print(best_model@model[["model_summary"]])


############### RANDOM ################
# Search a random subset of these hyper-parmameters. Max runtime 
# and max models are enforced, and the search will stop after we 
# don't improve much over the best 5 random models.
ntrees_opt = c(50, 100, 200, 500, 1000)
mtries_opt <- c(8, 10, 15,20)
max_depth_opt <- c(5, 10, 15, 20, 25)

hyper_params <- list(ntrees = ntrees_opt,
                     mtries = mtries_opt,
                     max_depth = max_depth_opt)

search_criteria = list(strategy = "RandomDiscrete",
                       max_runtime_secs = 600, #early stopping for RandomD
                       max_models = 50,  #early stopping for RandomD
                       stopping_metric = "AUTO", #early stopping for RandomD
                       stopping_tolerance = 0.001, 
                       stopping_rounds = 5, #after best 5 random models
                       seed = 123)

rf_gridR <- h2o.grid("randomForest",
                      x = x, y = y,
                     grid_id = "mygrid",
                     training_frame = train,
                     seed = 123,
                     nfolds = nfolds,
                     keep_cross_validation_predictions = TRUE,                    
                     hyper_params = hyper_params,
                     search_criteria = search_criteria)

rf_modelsR <- lapply(rf_gridR@model_ids, function(model_id) h2o.getModel(model_id))

rf_sorted_gridR <- h2o.getGrid(grid_id = "mygrid", sort_by = "accuracy", decreasing = T)
print(rf_sorted_gridR)

best_modelR <- h2o.getModel(rf_sorted_gridR@model_ids[[1]])
summary(best_modelR)

# Now let's evaluate the best model found performance on a test set
# so we get an honest estimate of top model performance
best_rf_perfR <- h2o.performance(model = best_modelR,
                                newdata = test)

h2o.confusionMatrix(best_rf_perfR)
h2o.auc(best_rf_perfR)
h2o.accuracy(best_rf_perfR)
h2o.sensitivity(best_rf_perfR,thresholds = 0.399634438119829)
h2o.specificity(best_rf_perfR,thresholds = 0.399634438119829)

# Look at the hyperparameters for the best model
print(best_modelR@model[["model_summary"]])


# ciclo for per ottenere 5 volte le prime 20 variabili per importanza
rf_20predictors_5 = list()
for (i in 1:5) {
  rf_fitbest_for_i<- h2o.randomForest(x = x,
                                      y = y,
                                      training_frame = data,
                                      nfolds = nfolds,
                                      ntrees = 200, #best tune
                                      mtries = 8, #best tune
                                      max_depth = 15) #best tune
  vi_i <- as.data.frame(h2o.varimp(rf_fitbest_for_i))
  vi20_i <- vi_i[1:20,] 
  vi20_i$i <- i  #  keep track of which iteration
  rf_20predictors_5[[i]] <- vi20_i # add it to list
}

big_data_rf = do.call(cbind, rf_20predictors_5)

write.csv(big_data_rf, file = "RF_varimp20_5")


### All done, shutdown H2O    
h2o.shutdown(prompt=FALSE)
