library(rrBLUP)
prob <- read.table("Imputed_hybrid_genotype_tGBS.filtermafboth.raw",sep=' ',head=T)
real <- read.table("Impute_hybrid_tGBS_realSNPmethods.filtermafbyprobability.raw",sep=' ',head=T)
wY <- read.table("inbred_2017_hybrid_PH_rmoutlier_BLUP.csv",sep=',',head=T)
rownames(prob) <- prob[,1]
rownames(real) <- real[,1]
rownames(wY) <- wY[,1]
idx <- intersect(wY$Genotype,rownames(prob))
wY <- wY[idx,]
prob <- prob[idx,]
real <- real[idx,]
#write.table(wY,file='inbred_2017_hybrid_PH_rmoutlier_BLUP_withgeno.csv')
a <- colnames(prob)
for (i in c(1:100)){
    set.seed(200+i)
    b <- sample(a[-1],2000)
    outfile1 <- paste('probability_',i,'_pred.txt',sep='')
    outfile2 <- paste('probability_',i,'_beta.txt',sep='')
    outfile3 <- paste('realSNP_',i,'_pred.txt',sep='')
    outfile4 <- paste('realSNP_',i,'_beta.txt',sep='')
    traingeno <- prob[,c('taxa',b)]
    predgeno <- prob[,c('taxa',b)]
    rownames(predgeno) <- prob[,1]
    realgeno <- real[,c('taxa',b)]
    predrealgeno <- real[,c('taxa',b)]
    rownames(predrealgeno) <- real[,1]
    ans <- kinship.BLUP(y=wY[,2],G.train = traingeno[,-1],G.pred = predgeno[,-1],K.method = "RR")
    ans1 <- kinship.BLUP(y=wY[,2],G.train = realgeno[,-1],G.pred = predrealgeno[,-1],K.method = "RR")
    print (ans$g.pred[1])
    kk <- as.data.frame(ans$g.pred)
    kk1 <- cbind(real$taxa,kk)
    ll <- as.data.frame(ans1$g.pred)
    ll1 <- cbind(real$taxa,ll)
    write.table(kk1,file=outfile1)
    write.table(ans$beta,file=outfile2)
    write.table(ll1,file=outfile3)
    write.table(ans1$beta,file=outfile4)
}
