#FINAL Project
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7093812/
setwd("D:\\JHU courses\\Gene Visiualization\\Final Project")
data <- read.table(file = "GSE137187_Count_normalized_matrix.txt" ,header =T, row.names = 1)
data <- as.data.frame(data) 

#Test for outlier samples and provide visual proof.  Remove these outliers. 

#knockout and control of BKLf8KO mouse islet for islet transcriptomic analysis 
#s961 insulin antagonist treatment

#take away the first column that has the gene symbol
dat <- data[,2:29]

#putting gene , knockout, and control into their own groups
symbol <- data[,1]
symbol <- as.data.frame(symbol)
symbol$ens <- rownames(dat)


#index of samples that are treated with S961
col.s961 <- grep("S961", colnames(dat))

#s961 treated
s961 <- subset(dat[col.s961])
#reg 
reg <- grep("S961",names(dat))
s961[1]

 
#cv mean
dat.pcat <- prcomp(t(dat)) 
plot(dat.pcat$x, main = "PCA plot of KO/Ctl ,S961/No-S961")
text(dat.pcat$x, label = dimnames(dat)[[2]],pos = 1, cex=0.8)

#Avg Corr of samples ; all samples
dat.cor <- cor(dat)
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of Control/KO samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

#CTL.2.S964

dat$CTL.2.S961 <- NULL


#filter out genes that have low expression values 

#filtering out genes 
??geneFilter
dat.rmeans <- apply(dat,1, mean)
dat.rmeans <- as.data.frame(dat.rmeans)

dat.sd <- apply(dat,1, sd)
dat.sd <-as.data.frame(dat.sd)
summary(dat.sd)
dat.sd[6183,]
dat.sd[2080,]

#cv = std/mean
cv = dat.sd/dat.rmeans
colnames(cv) <- "cv"
cv$cv
plot(1:nrow(cv),cv$cv)
summary(cv)
cv = abs(cv)
new.cv <- subset(cv, cv > 0.82 )


#removing the genes based on average cv
n <- match(rownames(new.cv), symbol$ens )
dat <- dat[n,]



hist(dat.rmeans$dat.rmeans)

hist(cv$cv, main= "CV of all Genes")
hist(new.cv$cv, main= "Genes meeting cv > 0.82")

#run t test and apply to all genes


t.test.all.genes <-function(x, s1,s2){
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  res <- t.test(x1,x2, paired  = FALSE)
  value <- as.numeric(res$p.value)
  #return(res)
  return (value)
  
}


res <- apply(dat,1, t.test.all.genes, s1= c(1:13) , s2 = c(14:27))
res.p <-  apply(dat,1, t.test.all.genes, s1= c(1:13) , s2 = c(14:27))


#subsetting pvals less than 0.05

res.pval <- subset(res.p, res.p <0.05)
res.pval <- as.data.frame(res.pval)

 
#genes of significance
#matching rownames to columnnames
res.name <- rownames(res.pval)
index <- match(res.name, symbol$ens)
p.sym <- symbol[index,]

#sign pval <0.05 data
filt.index <- match(res.name, dimnames(dat)[[1]])
dat <- dat[filt.index,]

#plot genes into a histogram
hist(res.pval$res.pval , main = "Genes of Significance p<0.05")

hist(res.p, main= "Unfiltered Pval for all genes with Student T")

#t-test with s961
colnames(s961)[1:8]
colnames(s961)[9:16]
s961.ttest <- apply(s961 , 1, t.test.all.genes, s1=c(1:8), s2=c(9:16))
s961.ttest.p <- subset(s961.ttest, s961.ttest <0.05)

#ANOVA test
??aov
aov.all.genes <- function(x,s1,s2) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  #x3 <- as.numeric(x[s3])
  fac <- c(rep('A',length(x1)), rep('B',length(x2)), rep('C',length(x3)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3))
  names(a.dat) <- c('factor','express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}


aov <- apply(dat.clust,1,s1= , s2=)


#subset  data by the genes that you determined  using clustering/ dendro
#clustering
# HCA heat map of top 100 genes
# adjust colors for heat map

dat.clust <- dat	

#red = upreg green= downreg
hm.rg <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000","#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")

heatmap(as.matrix(filt.dat),col=hm.rg)

in.dat<- match(rownames(dat.clust), symbol$ens)
symbol[in.dat,1]

rownames(dat.clust) <- in.dat
rownames(dat.clust)
dat.hca  <- hclust(dist(dat.clust, method = "manhattan"), method = "median")
plot(dat.hca)
rect.hclust(dat.hca, k=2, border = 'red')
clusters <- cutree(dat.hca,2)
clusters <- as.data.frame(clusters)

install.packages('factoextra')
library(factoextra)
dat <- na.omit(dat)
dat <- scale(dat)
distance <- get_dist(dat) 
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

s <- grep("KO", colnames(dat))
s <- dat[,s]
t <- grep("S961", colnames(s) )
s <-s[,t]
k <- kmeans(s , centers = 2 , nstart = 25)
fviz_cluster(k, data = s)
k2 <- kmeans(dat, centers = 2, nstart = 25)
str(k2)
fviz_cluster(k2 , data= dat)



#
#s961 treated

#reg 

s <- kmeans(s961.c, centers = 2, iter.max = 10)
plot(s , col = s$cluster, cex = 1)
points(s$centers)

s961.c <- grep("S961",names(dat)) 
s961.c <- as.matrix(dat[,s961.c])
#kmeans
install.packages("NbClust")
library(NbClust)
nb <- NbClust(filt.dat, diss= NULL, distance = "euclidean", method = "kmeans")
dd <- cbind(filt.dat, clusters$clusters)
cl <- kmeans(dat.clust,centers = 2, iter.max = 20)



# discriminant genes from classification
library(MASS)
names(dat)
clas <- names(dat)
clas[grep("CTL.\\d.S961",clas)] <- rep("c.t",length(clas[grep("CTL.\\d.S961",clas)]))
clas[grep("CTL.\\d",clas)] <- rep("c.r",length(clas[grep("CTL.\\d",clas)]))
clas[grep("KO.\\d.S961",clas)] <- rep("ko.t",length(clas[grep("KO.\\d.S961",clas)]))
clas[grep("KO.\\d",clas)] <- rep("ko.r",length(clas[grep("KO.\\d",clas)])) 


datx <- as.data.frame(t(dat))
datx <- data.frame(clas,datx)
dat.lda <- lda(clas~.,datx)
dat.pred <- predict(dat.lda,filt.dat)
plot(dat.pred$x,bg=as.numeric(factor(clas)),pch=21,col=1,ylab="Discriminant function",axes=F,xlab="Score",main="Discriminant function for Ko/Control Treated/Non treated-S9651 dataset")
axis(1)
axis(2)
text(dat.pred$x , label = dat.pred$class, pos = 1, cex = 0.9) 
table(dat.pred$class, clas) 






