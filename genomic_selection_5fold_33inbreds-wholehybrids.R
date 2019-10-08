library(rrBLUP)
hY <- read.table("hybrid_BLUP_rmoutlier/inbred_2017_hybrid_1kwt_rmoutlier_BLUP_withgeno.csv",sep=',',head=T)
hG <- read.table("ICRISAT_hybrids_GTPL_withpheno.raw",sep='\t',head=T)
iG <- read.table("ICRISAT.Pearl.Millet.tGBS.SNPs_Only-20160714.Recalled.MCR20.normal.GTPL.beagle.withpheno.parent.raw",sep='\t',head=T)
iY <- read.table("inbred_BLUP_rmoutlier/inbred_2017_1kwt_rmoutlier_BLUP_withgeno.csv",sep='\t',head=T)

a <- colnames(hG)
rownames(hG) <- hG[,1]
rownames(hY) <- hY[,1]
rownames(iG) <- iG[,1]
rownames(iY) <- iY[,1]

predpheno <- hY

g <- c()
for (n in 6:15){
  g <- c(g,2^n)
}
iY <- iY[rownames(iG),]
mylen <- length(a)-1
for (i in c(g,mylen)){
  m <- i
  write.table(m,file=paste("5fold-33inbreds-hybrid_",m,"1kwt.txt",sep=""))
  for (j in 1:20){
    set.seed(200+j)
    b <- sample(a[-1],m)
    sG <- hG[,c('Taxa',b)]
    predgeno <- sG[as.character(predpheno$Taxa),]
    isG <- iG[,c('Taxa',b)]
    q <- 2+j
    set.seed(q)
    isY <- iY[sample(rownames(iY)),]
    t <- cut(seq(1,nrow(isY)),breaks=5,labels=FALSE)
    for (k in 1:5){
      testIndexes <- which(t==k,arr.ind=TRUE)
      trainpheno <- isY[-testIndexes, ]
      traingeno <- isG[as.character(trainpheno$Taxa),]
      ans <- kinship.BLUP(y=trainpheno[,2],G.train = traingeno[,-1],G.pred = predgeno[,-1],K.method = "RR")
      mycor <- cor(ans$g.pred,predpheno[,2])
      write.table(mycor,file=paste("5fold-33inbreds-hybrid_",m,"1kwt.txt",sep=""),append=TRUE)
      print (mycor)
    }
  }
}
