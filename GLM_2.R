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
#data$Enrichment..0.control..1.case.<- as.factor(data$Enrichment..0.control..1.case.)
data <- data[,-c(1,2)]

########### sparsity ###########
# spars <- apply(X = data[,-1], 2, FUN = function(c) (sum(c > 0)/405*100))
# hist(spars,
#      xlim = c(0,100),
#      nclass = 100,
#      xlab = "Percent of samples", ylim = c(0,30), ylab = "Number of MLGs",
#      col = "red",
#      main = "Numero di MLG presenti in una data percentuale dei campioni")

# rimuovo high correlated
trainDescr <- data[,-1]
descrCorr <- cor(trainDescr)
highCorr <- findCorrelation(descrCorr, 0.85)
trainDescr <- trainDescr[, -highCorr]
ncol(trainDescr) #287 non high correlated predittori rimasti
colnames(trainDescr)
data <- data[, -highCorr]

#h2o
data <- as.h2o(data)
dim(data)
#h2o.levels(data$Enrichment..0.control..1.case.)

# Partition the data into training, validation and test sets
splits <- h2o.splitFrame(data = data, 
                         ratios = 0.8 ,  #partition data into 80%, 20%
                         seed = 1)
#setting a seed will guarantee reproducibility
train <- splits[[1]]
test <- splits[[2]]

# Take a look at the size of each partition
# Notice that h2o.splitFrame uses approximate splitting not exact splitting (for efficiency)
# so these are not exactly 80%, 20% of the total rows
nrow(train)  # 325
nrow(test)  # 80

# Identify response and predictor variables
y <- "Enrichment..0.control..1.case."
x <- setdiff(names(data), y)
family <- "binomial"
nfolds <- 10


# 1. Let's start with a basic binomial Generalized Linear Model (manually)
# By default, h2o.glm uses a regularized elastic net model
# Next we will do some automatic tuning by passing in a validation frame and setting 
# `lambda_search = True`.  Since we are training a GLM with regularization, we should 
# try to find the right amount of regularization (to avoid overfitting).  The model 
# parameter, `lambda`, controls the amount of regularization in a GLM model and we can 
# find the optimal value for `lambda` automatically by setting `lambda_search = TRUE` 
# and passing in a cross validation (which is used to evaluate model performance using a 
# particular value of lambda).

glm_fit <- h2o.glm(x = x, 
                   y = y, 
                   training_frame = train,
                   model_id = "glm_fit",
                   nfolds = nfolds,
                   family = family,
                   seed = 1,
                   standardize = T,
                   alpha = 0.5) #amount of regularization

#default alpha 0.5, amount of mixing between the two
#Lambda Regularization strength
summary(glm_fit)
glm_fit@parameters$alpha
glm_fit@parameters$lambda
glm_fit@parameters

#Cross Validation metrics
glm_fit@model$cross_validation_metrics
glm_fit@model$model_summary
glm_fit@model$lambda_best

# Let's compare the performance of the GLM on the test set
glm_perf <- h2o.performance(model = glm_fit,
                            newdata = test)
glm_perf
h2o.confusionMatrix(glm_perf)
h2o.sensitivity(glm_perf,thresholds = 0.456379868259154)
h2o.specificity(glm_perf,thresholds = 0.456379868259154)


#individuo predittori importanti
predictors <- glm_fit@model$coefficients_table #alcuni 0
active_pred <- subset(predictors, coefficients!=0)
inactive_pred <- subset(predictors, coefficients==0)


################# RANDOM GRID SEARCH ###########################

# This is set to run fairly quickly, increase max_runtime_secs 
# or max_models to cover more of the hyperparameter space.
# Also, you can expand the hyperparameter space of each of the 
# algorithms by modifying the hyper param code below.

search_criteria2 = list(strategy = "RandomDiscrete", #early stopping for RandomD
                        stopping_metric = "MSE", #early stopping for RandomD
                        stopping_tolerance = 0.0001,
                        stopping_rounds = 50,
                        seed = 123)

alpha_opt <- seq(0, 1, by = 0.1)
lambda_opt <-  seq(0, 1, by = 0.01)
hyper_params2 <- list(alpha = alpha_opt,
                      lambda = lambda_opt)

# Train and validate a grid of GLMs
glm_grid2 <- h2o.grid("glm", x = x, y = y,
                      grid_id = "glm_grid2",
                      training_frame = train,
                      nfolds = nfolds,
                      family = family,
                      seed = 1,
                      hyper_params = hyper_params2,
                      search_criteria = search_criteria2)

glm_gridperf2 <- h2o.getGrid(grid_id = "glm_grid2", 
                             sort_by = "mse")
print(glm_gridperf2)

best_glm2 <- h2o.getModel(glm_gridperf2@model_ids[[1]])
print(best_glm2@model[["model_summary"]])

#try best model
glm_fitbest <- h2o.glm(x = x, 
                       y = y, 
                       training_frame = train,
                       model_id = "glm_fitbest",
                       nfolds = nfolds,
                       family = family,
                       seed = 1,
                       standardize = T,
                       lambda = 0.16, #best value from grid
                       alpha = 0.05) #best value from grid


summary(glm_fitbest)
glm_fitbest@parameters$alpha
glm_fitbest@parameters$lambda

#Cross Validation metrics
glm_fitbest@model$cross_validation_metrics
glm_fitbest@model$model_summary
glm_fitbest@model$lambda_best

# Let's compare the performance of the GLM on the test set
glmbest_perf <- h2o.performance(model = glm_fitbest,
                                newdata = test)

glmbest_perf
h2o.confusionMatrix(glmbest_perf)
h2o.sensitivity(glmbest_perf,thresholds = 0.325579995761446)
h2o.specificity(glmbest_perf,thresholds = 0.325579995761446)

#individuo predittori importanti sull'intero dataset
# usando il miglior modello

glm_fitbest_all<- h2o.glm(x = x, 
                          y = y, 
                          training_frame = data,
                          model_id = "glm_fitbest_all",
                          nfolds = nfolds,
                          family = family,
                          standardize = T,
                          lambda = 0.16, #best value from grid
                          alpha = 0.05) #best value from grid

predictors <- glm_fitbest_all@model$standardized_coefficient_magnitudes
active_pred <- subset(predictors, coefficients!=0)
inactive_pred <- subset(predictors, coefficients==0)
dim(active_pred)
dim(inactive_pred)

active_pred[1:20,]


# ripeto 5 volte il miglior algoritmo per ottenere una matrice con le prime 20 MLG

glm_20predictors_5 = list()
for (i in 1:5) {
  glm_fitbest_for_i<- h2o.glm(x = x, 
                              y = y, 
                              training_frame = data,
                              #model_id = "glm_fitbest_all",
                              nfolds = nfolds,
                              family = family,
                              standardize = T,
                              lambda = 0.16, #best value from grid
                              alpha = 0.05) #best value from grid
  
  predictors_i <- glm_fitbest_for_i@model$standardized_coefficient_magnitudes
  active_pred_i <- as.data.frame(subset(predictors_i, coefficients!=0)[1:20,])
  active_pred_i$i <- i  # maybe you want to keep track of which iteration produced it?
  glm_20predictors_5[[i]] <- active_pred_i # add it to your list
}
big_data = do.call(cbind, glm_20predictors_5)



h2o.shutdown(prompt=FALSE)
