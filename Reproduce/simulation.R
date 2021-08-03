
library(matlabr)
library(R.matlab) 
library(MendelianRandomization) 
library(psych) 
library(psychTools)
library(factoextra)

######################### FUNCTION #########################

## FUNCTION: Main function to simulate data (based on July 16, 2021 meeting).
## Input: 
##    n: sample size
##    p: number of SNPs
##    b: non-zero coefficients (effect size) for B vector 
##    num.sig: number of SNPs with effect on mediator
##    q: number of observed mediators (i.e. M and M0 matrix)
##    rho: the correlation parameters to generate toeplitz covariance matrix for SNP's (i.e. X)
##    a: causal effect of Mf on outcome (i.e. Y)
##    sigma1: the standard deviation of normal random noises when generating M and M0
##    sigma2: the standard deviation of normal random noises when generating B_star

## Return: A list with all the simulated data and related infoation. 
simu.mvmr <- function(n=500, p=20, b=1, num.sig=10, q=40, a=1, 
                      sigma1=1, sigma2=1, rho=0) {
  X <- matrix(NA, n, p)
  snp_prob <- rep(1/3, 3) ## might change later on
  for (i in 1:p) {
    one_sum <- as.vector(rmultinom(n=1, size=n, prob=snp_prob))
    one_snp <- c(rep(0,one_sum[1]), rep(1,one_sum[2]), rep(2,one_sum[3]))
    one_snp <- one_snp[sample(1:n, n)]
    X[,i] <- one_snp ## n*p matrix
  }
  
  B <- c(rep(b,num.sig), rep(0,p-num.sig)) ## vector of length p
  Mf <- X%*%B + rnorm(n, sd=sigma1) ## n*1 matrix
  
  B_star <- matrix(rep(B,q/2), p, q/2) + rnorm(p*q/2, sd=sigma2) ## p*(q/2) matrix 
  B0_star <- matrix(0, p, q/2) ## p*(q/2) matrix
  
  M <- X %*% B_star + rnorm(n*q/2, sd=sigma1) ## n*(q/2) matrix 
  M0 <- X %*% B0_star + rnorm(n*q/2, sd=sigma1) ## n*(q/2) matrix 
  
  Y <- Mf * a + rnorm(n) ## vector of length n
  
  output <- list(X=X, B=B, Mf=Mf, B_star=B_star, B0_star=B0_star, 
                 M=M, M0=M0, Y=Y)
  return(output)
}



## FUNCTION: Run Factor Analysis, extrct g-score, and run linear regression with 
##           g-score and single SNP
## Input: 
##    dataM: mediators data
##    dataX: SNP data
## Return: A matrix with linear regression information 
fa.and.lm <- function(dataM, dataX) {
  #### Calculate g-factor of mediators
  searchN=fa.parallel(dataM, main="Parallel Analysis", plot=F)
  factorN=searchN[["nfact"]]
  fares<-omega(dataM, fm="pc", nfactors = factorN, plot=F) #rotate = "Promax", "oblimin",
  # summary(fares)
  gscore=fares$scores[,"g"]
  
  #### linear regression of gscore and single SNP
  gscore.lm=matrix(NA, nrow = ncol(dataX), ncol = 4) 
  colnames(gscore.lm)=c("Beta_exp", "Beta_SE_exp","Tvalue_exp", "Pvalue_exp")
  for (i in 1:ncol(dataX)) {
    fm=summary(lm(gscore~dataX[,i] ))
    gscore.lm[i,1]=fm$coefficients[2,1]
    gscore.lm[i,2]=fm$coefficients[2,2]
    gscore.lm[i,3]=fm$coefficients[2,3]
    gscore.lm[i,4]=fm$coefficients[2,4]
  }
  return(gscore.lm)
}



## FUNCTION: Run Principal Component Analysis, extrct first PC, and run linear regression  
##           with first PC and single SNP
## Input: 
##    dataM: mediators data
##    dataX: SNP data
## Return: A matrix with linear regression information 
pca.and.lm <- function(dataM, dataX) {
  #### Estimate first pc of mediators
  pcres <- prcomp(dataM, center=TRUE, scale=TRUE)
  pc1 <- pcres$x[,1]
  
  #### linear regression of first pc and single SNP
  pc1.lm=matrix(NA, nrow = ncol(dataX), ncol = 4 ) 
  colnames(pc1.lm)=c("Beta_exp", "Beta_SE_exp","Tvalue_exp", "Pvalue_exp")
  for (i in 1:ncol(dataX)) {
    fm=summary(lm(pc1~dataX[,i]))
    pc1.lm[i,1]=fm$coefficients[2,1]
    pc1.lm[i,2]=fm$coefficients[2,2]
    pc1.lm[i,3]=fm$coefficients[2,3]
    pc1.lm[i,4]=fm$coefficients[2,4]
  }
  
  return(pc1.lm)
}


######################### Parameter setting  #########################
## generate seed for each replication
set.seed(2021)
seed_vec <- sample(1:1e4, 1000, replace=F)

## parameters that change
start_rep <- 1
end_rep <- 500
a <- 1 ## causal effect on outcome, 0.5 or 1
print(paste("a =", a))

## parameters that are fixed
n <- 500
p <- 20
b <- 2
num_sig <- 10
q <- 40
sigma1 <- 1
sigma2 <- 1
alpha <- 1e-5 ## p-value cut-off for mediator~SNP p-value


start_time <- Sys.time()
for (idx in start_rep:end_rep) {
  print(paste(idx, "-th replication starts", sep=""))
  set.seed(seed_vec[idx])
  
  ######################### 1. Simulate data #########################
  sim_data <- simu.mvmr(n=n, p=p, num.sig=num_sig, q=q, a=a, b=b, 
                        sigma1=sigma1, sigma2=sigma2)
  
  dataX <- sim_data$X
  dataY <- as.vector(sim_data$Y)
  dataM <- cbind.data.frame(sim_data$M, sim_data$M0)
  
  ######################### 2. MR-single IV #########################
  
  ## 2a. mediator ~ SNP regression
  colnames(dataM)[1:(q/2)]=paste("M", 1:(ncol(dataM)/2), sep="")
  colnames(dataM)[(q/2+1):q]=paste("M0_", 1:(ncol(dataM)/2), sep="")
  
  ulmres.med <- vector("list", q)
  for (j in 1:q) {
    med=dataM[,j]
    temp.res=matrix(NA, p, 4) 
    colnames(temp.res)=c("Beta_exp", "Beta_SE_exp","Tvalue_exp", "Pvalue_exp")
    for (i in 1:p) {
      fm=summary(lm(med~dataX[,i]))
      temp.res[i,1]=fm$coefficients[2,1]
      temp.res[i,2]=fm$coefficients[2,2]
      temp.res[i,3]=fm$coefficients[2,3]
      temp.res[i,4]=fm$coefficients[2,4]
    }
    
    ulmres.med[[j]]=temp.res
    names(ulmres.med)[[j]]=colnames(dataM)[j]
  }
  
  snp_med_pval <- sapply(ulmres.med, function(x){x[,4]})
  snp_med_neglogP <- apply(snp_med_pval, 2, function(x){-log10(x)})
  
  ## 2b. outcome ~ SNP regression
  ulmres.out=matrix(NA, p, 4) 
  colnames(ulmres.out)=c("Beta_out", "Beta_SE_out","Tvalue_out", "Pvalue_out")
  for (i in 1:ncol(dataX)) {
    fm=summary(lm(dataY ~ dataX[,i] ))
    ulmres.out[i,1]=fm$coefficients[2,1]
    ulmres.out[i,2]=fm$coefficients[2,2]
    ulmres.out[i,3]=fm$coefficients[2,3]
    ulmres.out[i,4]=fm$coefficients[2,4]
  }
  
  
  ## 2c. univariate MR for each SNP-mediator pair
  uniMR.1snp=c() 
  uniMR.1snp.pvalue <- matrix(NA, p, q)
  uniMR.1snp.fullresults <- vector("list", p) ## MR fitting results for all SNP vs. all fa
  for (h in 1:p) {  ## MR using single SNP
    uniMR=data.frame(matrix(NA, q, 5))
    colnames(uniMR)=c("Estimate","StdError","CILower","CIUpper","Pvalue")
    rownames(uniMR)=names(ulmres.med)
    uniMR.fullresults <- vector("list", q) ## MR fitting results for one SNP vs. all fa
    for (k in 1:q) {
      MRobj=mr_ivw(mr_input(bx = ulmres.med[[k]][h,1], bxse=ulmres.med[[k]][h,2],
                            by = ulmres.out[h,1], byse = ulmres.out[h,2]))
      uniMR[k,1]=MRobj$Estimate
      uniMR[k,2]=MRobj$StdError
      uniMR[k,3]=MRobj$CILower
      uniMR[k,4]=MRobj$CIUpper
      uniMR[k,5]=MRobj$Pvalue
      
      uniMR.fullresults[[k]]=MRobj
      names(uniMR.fullresults)[k]=names(ulmres.med)[k]
    }
    uniMR.1snp[[h]]=uniMR
    uniMR.1snp.pvalue[h,] <- uniMR$Pvalue
    uniMR.1snp.fullresults[[h]]=uniMR.fullresults
  }
  
  uniMR.1snp.neglogP <- apply(uniMR.1snp.pvalue, 2, function(x){-log10(x)})
  colnames(uniMR.1snp.neglogP) <- colnames(dataM)
  
  ## B-H correction 
  uniMR_bh <- matrix(p.adjust(uniMR.1snp.pvalue, method="BH"), p, q)
  uniMR_bh_cut <- (uniMR_bh < alpha)
  uniMR_sel_col <- which(colSums(uniMR_bh_cut)>0)
  uniMR_sel_row <- which(rowSums(uniMR_bh_cut)>0)
  
  comb2_neglogP <- uniMR.1snp.neglogP ## combine two p-values
  comb2_neglogP[snp_med_pval>alpha] <- 0
  comb2_neglogP <- comb2_neglogP + matrix(runif(p*q,0,1e-10), p, q)

  
  ######################### 3. LAS submatrix #########################
  ## save the data for passing to Matlab 
  write.table(comb2_neglogP, row.names=F, col.names=F, sep=",",
              file=paste("/home-net/home-3/kehj@umd.edu/my_code/MR_multi/for_matlab/las_input_",
                         idx, ".csv", sep=""))
  
  ## use LAS matlab function in the mtba package 
  code <- c("cd '/your/own/directory/';", 
            paste("W=readtable('las_input_", idx, ".csv');", sep=""),
            "W1=table2array(W);",
            "addpath('/your/own/directory/mtba/');",
            "res=LAS(W1, 1);", 
            "rowSel=res.RowxNum;",
            "colSel=res.NumxCol;",
            paste("save('sub_matrix_",idx,".mat','rowSel','colSel','res')", sep=""))
  
  run_matlab <- run_matlab_code(code, verbose=FALSE)
  
  ## get the results back from Matlab
  fa.cluster <- readMat(paste("/your/own/directory/sub_matrix_",
                              idx, ".mat", sep=""))
  
  subdataM <- dataM[,fa.cluster$colSel==1]
  subdataX <- dataX[,fa.cluster$rowSel==1]
  
  las_sel_col <- which(fa.cluster$colSel==1)
  las_sel_row <- which(fa.cluster$rowSel==1)
  print(las_sel_col)
  print(las_sel_row)
  
  ######################### 4. PCA and factor analysis #########################
  ## 4.1 all mediators: factor analysis + lm
  gscore_lm_all <- fa.and.lm(dataM, dataX)
  
  ## 4.2 selected mediators (subset): factor analysis + lm
  gscore_lm_sub <- fa.and.lm(subdataM, dataX)
  
  ######################### 6. MR: multiple IVs #########################
  ## 5.1 all mediators: gscore and SNP
  mmr_gs_all <- mr_allmethods(mr_input(bx = gscore_lm_all[,1], bxse=gscore_lm_all[,2],
                                       by = ulmres.out[,1], byse = ulmres.out[,2] ), 
                              method = "all")
  
  ## 5.2  selected mediators (subset): gscore and SNP
  mmr_gs_sub <- mr_allmethods(mr_input(bx = gscore_lm_sub[,1], bxse=gscore_lm_sub[,2],
                                       by = ulmres.out[,1], byse = ulmres.out[,2]), 
                              method = "all")
  
  print(mmr_gs_all@Values[4,]) ## beta estiamtes by IVW
  print(mmr_gs_sub@Values[4,])
  
  ######################### Save the results #########################
  save(mmr_gs_all, mmr_gs_sub, uniMR.1snp.fullresults, 
       uniMR.1snp.pvalue, uniMR_bh, fa.cluster, 
       uniMR_sel_col, uniMR_sel_row, las_sel_col, las_sel_row,
       file=paste("/your/own/directory/July31_simumvmr_idx",idx,
                  "_n",n,"_p",p,"_q",q,"_b",b,"_a",a,"_sigma",10*sigma1,"_alpha",-log10(alpha),
                  ".RData", sep=""))
}
end_time <- Sys.time()
print(end_time - start_time)
