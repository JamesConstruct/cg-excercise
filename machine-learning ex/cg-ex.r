library(reshape2)
library(caret)
library(rpart.plot)
library(class)
library(pander)
library(neuralnet)
library(randomForest)


set.seed(1024) # reproducibilty
setwd("/")


normalize = function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

predictor.load = function(path) {

    cat(paste("Loading from path: ", path, "\n\n"))
    df = read.table(path, header=T)
    # df <- as.data.frame(lapply(df, normalize))

    rownames(df) = df$id
    df = df[,!colnames(df) %in% c("chr_start_end","id")]
    unknown = df[is.na(df$status),]   # unlabeled data
    df = df[!is.na(df$status),]       # labeled data
    df$status = factor(df$status)
    #print("Labeled data")
    #print(table(df$status))

    # Split training and testing
    train = createDataPartition(y = df$status, p=0.8, list=F)
    training = df[train,]
    #print("Training dimensions")
    #print(dim(training))
    #print("Training status")
    #print(table(training$status))
    testing = df[-train,]
    #print("Testing dimensions")
    #print(dim(testing))
    #print("Testing status")
    #print(table(testing$status))

    out = list(df, train, training, testing, unknown)
    names(out) = list("df", "train", "training", "testing", "unknown")

    return(out)

}

predictor.build.prototype = function(method, data, number=10, repeats=3, tuneLength=10) {

    trctrl = trainControl(method="repeatedcv", number=number, repeats=repeats)
    fit = train(status ~ ., data=data$training, method=method, trControl=trctrl, tuneLength=tuneLength)
    return(fit)

}

predictor.build.tree = function(data) {

    trctrl = trainControl(method="repeatedcv", number=10, repeats=3)
    fit = train(status ~ ., data=data$training, method="rpart", parms = list(split="information"), trControl=trctrl, tuneLength=10, control = list(maxdepth = 10))
    prp(fit$finalModel, tweak=1.2)

    return(fit)

}
predictor.build.knn = function(data) {

    return(predictor.build.prototype("knn", data))

}

predictor.build.nb = function(data) {

    return(predictor.build.prototype("naive_bayes", data))

}

predictor.build.rf = function(data) {

    return(predictor.build.prototype("rf", data))

}

predictor.build.svm = function(data) {

    return(predictor.build.prototype("svmLinear", data))

}

predictor.build.lda = function(data) {

    return(predictor.build.prototype("lda2", data))

}

predictor.build.nn = function(data) {

    fit = neuralnet(status ~ ., data=data$training, hidden=c(2, 2),err.fct = "ce",algorithm = "rprop+", lifesign = 'full', act.fct = "logistic", linear.output = FALSE, stepmax = 100000)
    return(fit);

}

predictor.build.rf.tune = function(data) {
    
    trctrl = trainControl(method="repeatedcv", number=10, repeats=3)
    model = randomForest(formula = status ~ ., data = data$training, trControl=trctrl)

    plot(model)
    varImpPlot(model)

    model_tuned <- tuneRF(
        x=data$training[,-15], # exclude status
        y=data$training$status, 
        ntreeTry=260,
        mtryStart=2, 
        stepFactor=1.5,
        improve=0.01,
        trace=FALSE # don't show real-time progress
    )
    # print(model)

}

predictor.build.rfImproved = function(data) {
    
    trctrl = trainControl(method="repeatedcv", number=10, repeats=3)
    model = randomForest(formula = status ~ ., data = data$training, trControl=trctrl, ntree = 500, mtry=3, tuneLength=100)

    return(model)

}


predictor.evaluate = function(fit, data) {

    pred = predict(fit, newdata=data$testing)
    cm <- table(data.frame(prediction=pred, truth=data$testing$status))
    pander(cm, caption="Confusion Matrix")

    pred = predict(fit, newdata=data$testing)
    print(caret::confusionMatrix(pred, data$testing$status, positive="1"))

}

predict.nn = function(fit, newdata) {   # convert linear NN output to binary

    p = compute(nn, newdata)
    prob = p$net.result
    pred <- ifelse(prob>0.5, 1, 0)
    return(pred)

}

predictor.predict = function(fit, data) {

    # Apply to unknown data
    pred = predict(fit, newdata=data$unknown)
    print("Outcome prediction of unlabeled data")
    print(table(pred))
    data$unknown$status = pred

    # Write the predicted data to a file
    data$df = rbind(data$df, data$unknown)
    data$df$id = rownames(data$df)
    write.table(data$df, "predictions.tsv", quote=F, row.names=F, col.names=T, sep="\t")

}

predictor.testAll = function() {

    dataPath = "deletion.tsv.gz"
    data = predictor.load(dataPath)

    tree = predictor.build.tree(data) # 0.9512
    knn = predictor.build.knn(data)  # 0.9338
    nb = predictor.build.nb(data) # 0.9443
    rf = predictor.build.rf(data) # 0.9582
    svm = predictor.build.svm(data) # 0.9547
    lda = predictor.build.lda(data) # 0.9512
    rfImp = predictor.build.rfImproved(data)
    # nn = predictor.build.nn(data)   # requires to change status from factor to numeric
    predictor.evaluate(tree, data)
    predictor.evaluate(knn, data)
    predictor.evaluate(nb, data)
    predictor.evaluate(rf, data)
    predictor.evaluate(svm, data)
    predictor.evaluate(lda, data)
    predictor.evaluate(rfImp, data)

}

predictor.processAll = function() {

    dataPath = "deletion.tsv.gz"
    data = predictor.load(dataPath)

    rfImp = predictor.build.rfImproved(data)
    predictor.evaluate(rfImp, data)

    predictor.predict(rfImp, data)

}



#      SIMPLE RF EVALUATION
#           Reference
# Prediction   0   1
#          0  45   2
#          1   5 235
                                          
#                Accuracy : 0.9756          
#                  95% CI : (0.9504, 0.9901)
#     No Information Rate : 0.8258          
#     P-Value [Acc > NIR] : 8.599e-16       
                                          
#                   Kappa : 0.9132          
                                          
#  Mcnemar's Test P-Value : 0.4497          
                                          
#             Sensitivity : 0.9916          
#             Specificity : 0.9000          
#          Pos Pred Value : 0.9792          
#          Neg Pred Value : 0.9574          
#              Prevalence : 0.8258          
#          Detection Rate : 0.8188          
#    Detection Prevalence : 0.8362          
#       Balanced Accuracy : 0.9458          
                                          
#        'Positive' Class : 1 


#    IMPROVED RF EVALUATION

#           Reference
# Prediction   0   1
#          0  47   2
#          1   3 235
                                          
#                Accuracy : 0.9826          
#                  95% CI : (0.9598, 0.9943)
#     No Information Rate : 0.8258          
#     P-Value [Acc > NIR] : <2e-16          
                                          
#                   Kappa : 0.939           
                                          
#  Mcnemar's Test P-Value : 1               
                                          
#             Sensitivity : 0.9916          
#             Specificity : 0.9400          
#          Pos Pred Value : 0.9874          
#          Neg Pred Value : 0.9592          
#              Prevalence : 0.8258          
#          Detection Rate : 0.8188          
#    Detection Prevalence : 0.8293          
#       Balanced Accuracy : 0.9658          
                                          
#        'Positive' Class : 1 