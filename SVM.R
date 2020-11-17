#SVM
library(caret)
library(h2o)
library(plyr)
library(dplyr)
library(e1071)
library(pROC)
h2o.init(nthreads = -1)

#upload dataset
data <- read.csv("MyData_ACVD.csv", stringsAsFactors = T)
as.data.frame(data)
data$Enrichment..0.control..1.case. <- as.factor(data$Enrichment..0.control..1.case.)
data <- data[,-c(1,2)]

#Qualità del dato

########### COR ###################
#High correlated vars
trainDescr <- data[,-1]
descrCorr <- cor(trainDescr)
highCorr <- findCorrelation(descrCorr, 0.85) #396 high correlated
trainDescr <- trainDescr[, -highCorr]
ncol(trainDescr) #287 non high correlated predittori rimasti
colnames(trainDescr)
data <- data[, -highCorr]

###### NZV ############
#The function nearZeroVar can be used to identify near zero-variance predictors in a dataset
# the percentage of unique values is less than 20% and
#the ratio of the most frequent to the second most frequent value is greater than 20
# nzv <- nearZeroVar(trainDescr, uniqueCut = 20, freqCut = 20, names = F, saveMetrics = F)
# nzv_pred <- rownames(nzv[which(nzv$nzv==T),]); nzv_pred
# 
# data <- data[,-nzv] #284 predittori rimasti

# # Partition the data into training and test set
# # setting a seed will guarantee reproducibility
splits <- h2o.splitFrame(data = as.h2o(data), 
                         ratios = 0.8 , #partition data into 80%, 20%
                         seed = 1)
train <- as.data.frame(splits[[1]])
test <- as.data.frame(splits[[2]])

# Take a look at the size of each partition
dim(train)
dim(test)
table(test$Enrichment..0.control..1.case.)

# # ######### LINEAR #######################
# library(e1071)
# svmfit = svm(Enrichment..0.control..1.case. ~ ., data = train, kernel = "linear", cost = 10, scale = T)
# # #You tell SVM that the kernel is linear, the tune-in parameter cost is 10, and scale equals false. 
# 
# #gives its summary
# print(svmfit)
# # #you can see that the number of support vectors is 6:
# # #they are the points that are close to the boundary or on the wrong side of the boundary.
# 
# # ########## NON lINEAR ####################
# fit = svm(Enrichment..0.control..1.case. ~ ., data = train, scale = FALSE, kernel = "radial", cost = 5)
# print(fit)
# # 
# # 
# # ################### GRID SEARCH #######################
# # We can further improve our SVM model and tune it so that the error is even lower. 
# # We will now go deeper into the SVM function and the tune function. 
# # We can specify the values for the cost parameter and epsilon which is 0.1 by default. 
# # A simple way is to try for each value of epsilon between 0 and 1 (I will take steps of 0.01) 
# # and similarly try for cost function from 4 to 2^9 (I will take exponential steps of 2 here). 
# # I am taking 101 values of epsilon and 8 values of cost function. 
# # I will thus be testing 808 models and see which ones performs best. 
# # The code may take a short while to run all the models and find the best version.
# svm_tune <- tune(svm, y ~ x, data = train,
#                  ranges = list(epsilon = seq(0,1,0.01), cost = 2^(2:9)))
# print(svm_tune)
# #Printing gives the output:
# #Parameter tuning of ???svm???:
# # - sampling method: 10-fold cross validation 
# #- best parameters:
# # epsilon cost
# #0.09 256
# #- best performance: 2.569451
# #This best performance denotes the MSE. The corresponding RMSE is 1.602951 which is square root of MSE
# 
# 
# #The best model
# best_mod <- svm_tune$best.model
# best_mod_pred <- predict(best_mod, train) 
# #We can now see the improvement in our model by calculating its RMSE error using the following code
# error_best_mod <- train$y - best_mod_pred 
# # this value can be different on your computer
# # because the tune method randomly shuffles the data
# best_mod_RMSE <- sqrt(mean(error_best_mod^2)) # 1.290738
# 
# #We can also plot our tuning model to see the performance of all the models together
# plot(svm_tune)
# 
# # This plot shows the performance of various models using color coding. 
# # Darker regions imply better accuracy. 
# # The use of this plot is to determine the possible range where we can narrow down our search to and try further tuning if required. 



####### TRAINING E TESTING LINEAR ##########
#trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
trctrl <- trainControl(method = "cv", number = 10)
set.seed(3233)

svm_Linear <- train(Enrichment..0.control..1.case. ~ . , data = train, method = "svmLinear",
                    trControl = trctrl,
                    preProcess = c("center", "scale"))

svm_Linear

# test_pred <- predict(svm_Linear, newdata = test)
# confusionMatrix(test_pred, test$Enrichment..0.control..1.case.)

########## GRID SEARCH ##############
#grid for C cost
grid <- expand.grid(C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5, 10, 100, 500))
set.seed(3233)
svm_Linear_Grid <- train(Enrichment..0.control..1.case. ~., data = train, method = "svmLinear",
                           trControl=trctrl,
                           preProcess = c("center", "scale"),
                           tuneGrid = grid)

svm_Linear_Grid

svm_Linear_Grid$


plot(svm_Linear_Grid)
#The above plot is showing that our classifier is giving best accuracy starting from C = 0.01. 

test_pred_grid <- predict(svm_Linear_Grid, newdata = test)
confusionMatrix(test_pred_grid, test$Enrichment..0.control..1.case., positive = "1")



####### TRAINING E TESTING NON LINEAR: RBF kernel ##########
set.seed(3233)
svm_Radial <- train(Enrichment..0.control..1.case. ~., data = train, method = "svmRadial",
                      trControl=trctrl,
                      preProcess = c("center", "scale"),
                      tuneLength = 10)
svm_Radial
plot(svm_Radial)

test_pred_Radial <- predict(svm_Radial, newdata = test)
confusionMatrix(test_pred_Radial, test$Enrichment..0.control..1.case., positive = "1")

########## GRID SEARCH ##############
#grid for C cost e sigma
grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
                                     0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                           C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                                 1, 1.5, 2,5,10))
set.seed(3233)
svm_Radial_Grid <- train(Enrichment..0.control..1.case. ~., data = train, method = "svmRadial",
                           trControl=trctrl,
                           preProcess = c("center", "scale"),
                           tuneGrid = grid_radial)

svm_Radial_Grid
plot(svm_Radial_Grid)
plot(svm_Radial_Grid, plotType = "level")

svm_Radial_Grid$bestTune
subset(svm_Radial_Grid$results, sigma == 0.01 & C == 1.5)

#Let's check our trained models accuracy on the test set.
test_pred_Radial_Grid <- predict(svm_Radial_Grid, newdata = test)
confusionMatrix(test_pred_Radial_Grid, test$Enrichment..0.control..1.case. , positive = "1")


####### TRAINING E TESTING POLINOMIAL KERNEL ##########
set.seed(3233)
svm_Poly <- train(Enrichment..0.control..1.case. ~., data = train, method = "svmPoly",
                    trControl=trctrl,
                    preProcess = c("center", "scale"))
svm_Poly
plot(svm_Poly)
plot(svm_Poly, metric = "Kappa")
svm_Poly$bestTune
#A level plot of the accuracy values 
plot(svm_Poly, plotType = "level")
resampleHist(svm_Poly)

test_pred_Poly <- predict(svm_Poly, newdata = test)
confusionMatrix(test_pred_Poly, test$Enrichment..0.control..1.case.)

##########GRID SEARCH##############
#grid for C cost, scale e polynomial degree
grid_poly <- expand.grid(scale = c(0,0.01, 0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                        C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,1, 1.5, 2, 5, 10),
                        degree = c(1,2,3))
set.seed(3233)
svm_Poly_Grid <- train(Enrichment..0.control..1.case. ~., data = train, method = "svmPoly",
                         trControl=trctrl,
                         preProcess = c("center", "scale"),
                         tuneGrid = grid_poly)

svm_Poly_Grid
plot(svm_Poly_Grid)
svm_Poly_Grid$bestTune
subset(svm_Poly_Grid$results, degree == 2 & scale == 0.01 & C== 0.25)

test_pred_Poly_Grid <- predict(svm_Poly_Grid, newdata = test)
confusionMatrix(test_pred_Poly_Grid, test$Enrichment..0.control..1.case. , positive = "1")

#The SVM model hasn't' a built-in variable importance score
#the varImp can still be used to get scores
#For SVM classification models, the default behavior is to compute the area under the ROC curve.

#ROC curve variable importance
roc_imp_lin <- varImp(svm_Linear_Grid, scale = T)
roc_imp_pol <- varImp(svm_Poly_Grid, scale = T)
roc_imp_rad <- varImp(svm_Radial_Grid, scale = T)

plot(roc_imp_lin, top = 20, main = "Importanza delle variabili SVM con kernel lineare")
plot(roc_imp_pol, top = 20, main = "Importanza delle variabili SVM con kernel polinomiale di grado 2")
plot(roc_imp_rad, top = 20, main = "Importanza delle variabili SVM con kernel RBF")

# #### compare 3 models ######
# models <- list(svmlin = svm_Linear_Grid, svmpol = svm_Poly_Grid, svmrad = svm_Radial_Grid)
# testPred <- predict(models, newdata = test)
# 
# predValues <- extractPrediction(models, testX = test[,-1], testY = test[,1])
# testValues <- subset(predValues, dataType == "Test")
# trainValues <- subset(predValues, dataType == "Training")
# head(testValues)
# head(trainValues)

### All done, shutdown H2O    
h2o.shutdown(prompt=FALSE)
