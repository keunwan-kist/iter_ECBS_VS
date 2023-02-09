#!~/opt/bin/ Rscript

# Keunwan Park

options<-commandArgs(trailingOnly=T)

if(length(options) < 2) stop("Invalid argument number\n\nRscript [].r [DrugPair train set] [model file name output] [ntree]\n")
if(!is.na(options[1])) train_data_fn=options[1]
if(!is.na(options[2])) output_model_fn = options[2]

ntree=200
if(!is.na(options[3])) ntree = as.numeric(options[3])

train_prefix = tools::file_path_sans_ext(basename(train_data_fn))

train_data <- read.table(train_data_fn)
train_data_name <- train_data[,c(1,2)]
train_data <- train_data[,-c(1,2)]
lidx = ncol(train_data)
train_data[,lidx] = as.factor(train_data[,lidx])
colnames(train_data)[lidx]="class"

weights = as.numeric(train_data$class)
pos_idx <- which(weights==2)
neg_idx <- which(weights==1)
# pos x 6 = neg
weights[pos_idx] = 1
neg_rate = round(length(pos_idx) / length(neg_idx),2)
weights[neg_idx] = neg_rate
print(paste("Negative rate :",neg_rate))

#################### training and validation #######################
require(ranger)

cat("Data preparation finished...now making a classification model by RANGER\n")
data.rb <- ranger(class~., data=train_data, case.weights = weights, num.trees=ntree, probability=T,save.memory=T)

cat("Ranger finished... Now saving the randomforest model ---> [data.rb] \n")
save(data.rb, file = output_model_fn) 

#load(file = model_name)
#load(file = roc_name)






