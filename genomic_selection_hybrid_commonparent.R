library(rrBLUP)
hY <- read.table("hybrid_BLUP_rmoutlier/inbred_2017_hybrid_1kwt_rmoutlier_BLUP_withgeno.csv",sep=',',head=T)
hG <- read.table("ICRISAT_hybrids_GTPL_withpheno.raw",sep='\t',head=T)
mf <- read.table("B_and_R.txt",sep='\t',head=T)

rownames(hY) <- hY[,1]
rownames(hG) <- hG[,1]

a <- colnames(hG)
l <- length(x)

g <- c()
for (n in 6:15){
  g <- c(g,2^n)
}
mylen <- length(a)-1
for (i in c(g,mylen)){
  m <- i
  write.table(m,file=paste("hybrid-hybrid_nonparent-parent_",m,"1kwt.txt",sep=""))
  for (j in 1:20){
    set.seed(200+j)
    b <- sample(a[-1],m)
    set.seed(100+j)
    B <- sample(na.omit(mf$Bline),1)
    set.seed(10+j)
    R <- sample(na.omit(mf$Rline),1)
    sB <- rownames(hY)[grep(B,rownames(hY))]
    sR <- rownames(hY)[grep(R,rownames(hY))]
    test <- c(sB,sR)
    test <- unique(test)
    predgeno <- hG[test,]
    predgeno <- predgeno[,c('Taxa',b)]
    predpheno <- hY[test,]
    train1 <- rownames(hY)[is.na(pmatch(rownames(hY),test))]
    q <- 2+j
    set.seed(q)
    trainY <- hY[sample(train1),]
    t <- cut(seq(1,nrow(trainY)),breaks=5,labels=FALSE)
    for (k in 1:5){
      testIndexes <- which(t==k,arr.ind=TRUE)
      traingeno <- hG[-testIndexes,]
      traingeno <- traingeno[,c('Taxa',b)]
      trainpheno <- hY[as.character(rownames(traingeno)),]
      ans <- kinship.BLUP(y=trainpheno[,2],G.train = traingeno[,-1],G.pred = predgeno[,-1],K.method = "RR")
      mycor <- cor(ans$g.pred,predpheno[,2])
      write.table(mycor,file=paste("hybrid-hybrid_nonparent-parent_",m,"1kwt.txt",sep=""),append=TRUE)
      print (mycor)
    }
  }
}