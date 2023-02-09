#!~/opt/bin/ Rscript

# Keunwan Park

options<-commandArgs(trailingOnly=T)

if(length(options) < 2) stop("Invalid argument number\n\nRscript [].r [DrugPair train set] [DrugPair test set]\n")
if(!is.na(options[1])) train_data_fn=options[1]
if(!is.na(options[2])) test_data_fn=options[2]

ntree=200 
if(!is.na(options[3])) ntree=as.numeric(options[3])

train_prefix = tools::file_path_sans_ext(basename(train_data_fn))
test_prefix = tools::file_path_sans_ext(basename(test_data_fn))

dir_path = dirname(train_data_fn)
print(dir_path)

train_data <- read.table(train_data_fn)
test_data <- read.table(test_data_fn)

train_data_name <- train_data[,c(1,2)]
test_data_name <- test_data[,c(1,2)]

train_data <- train_data[,-c(1,2)]
test_data <- test_data[,-c(1,2)]

lidx = ncol(train_data)

test_pos_idx <- which(test_data[,lidx]==1)
test_neg_idx <- which(test_data[,lidx]==0)

train_data[,lidx] = as.factor(train_data[,lidx])
test_data[,lidx] = as.factor(test_data[,lidx])

colnames(train_data)[lidx]="class"
colnames(test_data)[lidx]="class"


weights = as.numeric(train_data$class)
pos_idx <- which(weights==2)
neg_idx <- which(weights==1)
# pos x 6 = neg
weights[pos_idx] = 1
neg_rate = round(length(pos_idx) / length(neg_idx),2)
#neg_rate = 0.15
weights[neg_idx] = neg_rate/2
print(paste("Negative rate :",neg_rate))

#################### training and validation #######################
require(ranger)
require(PRROC)

cat("Data preparation finished...now making a classification model by RANGER\n")

#data.rb <- ranger(class~., data=train_data, case.weights = weights, num.trees=ntree, probability=T,save.memory=T)
data.rb <- ranger(class~., data=train_data, case.weights = weights, num.trees=ntree, probability=T,save.memory=T, min.node.size=3)

cat("Ranger finished... Now saving the randomforest model\n")
model_name = paste(train_data_fn,".ranger_rf_model",sep="")
#save(data.rb, file = model_name) 

print(data.rb)

cat("Testing Train-model to Test data...\n")

test_pred <- predict(data.rb,test_data[,-lidx]) 

roc <- roc.curve( test_pred$predictions[test_pos_idx,2], test_pred$predictions[test_neg_idx,2])
pr <- pr.curve( test_pred$predictions[test_pos_idx,2], test_pred$predictions[test_neg_idx,2])

cat( paste("Test set PR/ROC AUC value :",round(pr$auc.integral,5),round(roc$auc,5),"\n") )


write.table(cbind(train_data_name,round(data.rb$predictions[,2],5),train_data[,lidx]), file=paste(dir_path,"/oob.",train_prefix,".",test_prefix,".score",sep=''),quote=F,col.names=F,row.names=F) 
write.table(cbind(test_data_name,round(test_pred$predictions[,2],5),test_data[,lidx]), file=paste(dir_path,"/result.",train_prefix,".",test_prefix,".score",sep=''),quote=F,col.names=F,row.names=F) 

#load(file = model_name)
#load(file = roc_name)






