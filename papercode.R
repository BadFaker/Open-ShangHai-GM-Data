require(caret)
require(AMORE)
require(genalg)
## 数据准备
dataset <- read.csv("C:/Users/Lich/Desktop/dataset.csv",header=T)
Sample <- dataset[c("data_field","dep_type","open_type","refresh_rate","file_type")]
Species <- dataset["target"]
Score <- dataset[c("z_attention", "z_utilization")]
AU <- dataset["score"]
# 设置70%为训练集30%为观测集
require(caret)
Sample_Index <- createDataPartition(dataset$target, p = 0.7 , list = FALSE)
m <- nrow(Sample_Index)
n <- nrow(Species)
Train_Meas <- Sample[Sample_Index,]
Train_Species <- Species[Sample_Index,]
Train_Score <- Score[Sample_Index,]
Train_AU <- AU[Sample_Index,]
Test_Meas <- Sample[-Sample_Index, ]
Test_Species <- Species[-Sample_Index,]
Test_Score <- Score[-Sample_Index,]
## 神经网络部分
#激活函数
tansig <- function(x){
  return((1-exp(-2*x))/(1+exp(-2*x)))
}
#神经网络输入层，隐藏层，输出层参数设置
nclass <- c(5, 11, 2)
#神经网络正向子函数
net <- function(data, parameter){
  data <- scale(data)
  c_d <- ncol(data)
  #提取权值
  w_hide <- matrix(parameter[1:(c_d * nclass[2])], nrow = c_d, ncol = nclass[2])
  w_out <- parameter[(c_d * nclass[2] + 1):(c_d * nclass[2] + nclass[2] * nclass[3] )]
  #计算各个节点值
  in_value <- as.matrix(data)
  hide_value <- tansig(in_value %*% w_hide)
  out_value <- tansig(hide_value %*% w_out)
  return(out_value)
}
# 计算误差子函数
err <- function(parameter){
  out_value <-  net(data_train, parameter);
  error <- 0.5 * sum((out_value-data_train[, 15])^2);
  return(error);
}	
# 模型效果检验函数
P <- function(pred_score){
  p = 0;  K = 0;  r_a = 0;  r_u = 0;
  for(i in 1:(n-m)){
    K = K + ((Test_Score[i,1] - pred_score[i,1])^2 + (Test_Score[i,2] - pred_score[i, 2])^2 )^0.5
    r_a = r_a + abs((Test_Score[i,1] -pred_score[i,1])/(Test_Score[i,1]))
    r_u = r_u + abs((Test_Score[i,2] -pred_score[i,2])/(Test_Score[i,2]))
    if((pred_score[i,1]*pred_score[i,2] - mean(Train_AU))*(Test_Species[i] - 1.5) <= 0 ){
      p = p + 1
    }
  }
  return(c(p/(n-m), K/(n-m), r_a/(n-m), r_u/(n-m)))
}
## 构建标准BP神经网络
for (i in 1:m) {
  Train_Meas[i,] <- as.numeric(as.vector(Train_Meas)[i,])
  Train_Species[i] <- as.numeric(as.vector(Train_Species)[i])
  Train_Score[i,] <- as.numeric(as.vector(Train_Score)[i,])
}
for (i in 1:(n-m)) {
  Test_Meas[i,] <- as.numeric(as.vector(Test_Meas)[i,])
  Test_Species[i] <- as.numeric(as.vector(Test_Species)[i])
  Test_Score[i,] <- as.numeric(as.vector(Test_Score)[i,])
}
# 进行训练 
net <- newff(n.neurons=c(5,11,2), learning.rate.global=1e-3, momentum.global=1e-3,
             error.criterium="LMS", Stao=NA, hidden.layer="tansig", 
             output.layer="purelin", method="BATCHgd")
result <- train(net, Train_Meas, Train_Score, error.criterium="LMS", 
                report= F, show.step=100, n.shows=5)
# 预测实验
pred_score <- sim(result$net, Test_Meas)

## 遗传算法部分
#生成随机染色体
rndDNA = function(n)
{
  return(runif(100,0,1)) 
}
# 监控函数
monitor <- function(obj){
  xlim = c(0, 1);
  ylim = c(0, 1);
  net <- newff(n.neurons=c(5,11,2), learning.rate.global=obj$population[1,1], momentum.global=x[1,2],
               error.criterium="LMS", Stao=NA, hidden.layer="tansig", 
               output.layer="purelin", method="BATCHgd")
  result <- train(net, Train_Meas, Train_Score, error.criterium="LMS", report= F, show.step=100, n.shows=5)
  fy <- sim(result$net, Test_Meas);
  plot(Test_Score,  pch=20);
  points(x,l,type='l',col='blue',lwd=2);
  points(fy,pch = "o",col='red');
}
# 适应度函数
f<-function(x){
  net <- newff(n.neurons=c(5,11,2), learning.rate.global=x[1], momentum.global=x[2],
               error.criterium="LMS", Stao=NA, hidden.layer="tansig", 
               output.layer="purelin", method="BATCHgd")
  result <- train(net, Train_Meas, Train_Score, error.criterium="LMS", report= F, show.step=100, n.shows=5)
  fy <- sim(result$net, Test_Meas)
  K = 0
  for(i in 1:nrow(fy)){
    K = K + ((fy[i,1]-Test_Score[i,1])^2 + (fy[i,2]- Test_Score[i,2])^2)^0.5
  }
  return (K) 
}
# 种群交叉函数
cross <- function(i,nGroup,crossGroup,prop)
{
  a = nGroup[i,];  b = crossGroup[i,];  n = length(a);  m = max(1,trunc(n*prop))
  tmpa = a;  st = sample(1:n,1);  ind = st:(st + m) ;    
  if (st + m>n) {
    bind = which(ind > n)
    ind[bind] = ind[bind] %% n + 1
  }
  cross = intersect(b, a[ind])
  tmpa[ind] = cross
  return(tmpa)
}
# 种群变异函数
ariation <- function(dna,prop){
  n = length(dna);  pos = which(runif(n)<prop)
  if (length(pos)==0){
    return(dna)
  }
  pos = sample(pos,1);  newind = sample((1:n)[-pos],1)
  if (pos>newind) {
    tmp = dna[newind:n]
    tmp = tmp[-(pos-newind+1)]
    dna[newind] = dna[pos]
    dna[(newind+1):n] = tmp
  }
  else
  {
    tmp = dna[1:newind]
    tmp = tmp[-pos]
    dna[newind] = dna[pos]
    dna[1:(newind-1)] = tmp
  }
  return(dna)
}
# 参数优化函数
ga_bp <- rbga(matrix(0,1,ncrol(parameter)), matrix(1,1,ncrol(parameter)), popSize =100, iters =1000, evalFunc = f, mutationChance =0.001, ,elitism = 0.3,verbose= T, monitorFunc = monitor)
ga_bp$population[1,]

## 作图部分
# 数据标注
x <- 0.01*c(1:100); l = mean(Train_AU) /x
plot(dataset$z_attention[which(dataset$target == 1)],dataset$z_utilization[which(dataset$target == 1)], pch = 19,
     xlim = c(0,1), ylim = c(0,1), 
     xlab = "z_attention", ylab = "z_utilizaion")
points(dataset$z_attention[which(dataset$target == 2)],dataset$z_utilization[which(dataset$target == 2)], pch = 19,col='red')
points(x,l,type = 'l',col = 'blue',lwd = 2, title(sub =  "l: Z(A)Z(U)-0.05618374 = 0"))
# 相关性检验
require(corrplot)
factor_cor <- cor(dataset[c("z_attention", "z_utilization", "data_field", "dep_type", "open_type", "refresh_rate", "file_type")],method = "spearman") 
corrplot(factor_cor, method = "number")
# 标准BP神经网络
plot(Test_Score, pch=20)
points(pred_score, col='red', pch = "o")
points(x,l,type='l',col='blue',lwd=2)
# 优化迭代误差
ga_bp <- rbga(matrix(0,1,ncrol(parameter)), matrix(1,1,ncrol(parameter)), popSize =100, iters =1000, evalFunc = f, mutationChance =0.001, elitism = 0.3, verbose= T, monitorFunc = monitor)
plot(ga_bp)
