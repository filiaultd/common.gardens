##heritbility functions

library(coxme)
library(HLMdiag)
library(lme4)

## ---- H2lmekin

H2lmekin=function(model, g="id"){
    varID=model$vcoef[g]
    varexp=model$vcoef
    varexp=varexp[names(varexp)!=g]
    H2=unlist(varID)/sum(unlist(varID), unlist(varexp),model$sigma^2)
    return(H2)
}
## ---- end-of-H2lmekin


## ---- varid

var_ID=function(m){
    require(HLMdiag)
    vc=varcomp.mer(m)
    H2=vc[2]/sum(vc)
    return(as.numeric(H2))
}
## ---- end-of-varid


get_CI=function(model, var_ID, nsim){
    boot=bootMer(model, var_ID, nsim=nsim, type = "parametric", use.u = F, ncpus=6, parallel="multicore")
    bias=mean(boot$t, na.rm=T)-mean(boot$t0, na.rm=T)
    CI=tryCatch(boot.ci(boot, type=c("perc")), error=function(){return(NA)})
    CI1=as.numeric(CI$percent[1,4])-bias
    CI2=as.numeric(CI$percent[1,5])-bias
    return(c(CI1, CI2))
}

panel.cor <- function(x, y, digits = 2, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor.test(x, y, method="spearman", use = "pairwise.complete.obs")
    txt1 <- format(c(r$estimate, 0.123456789), digits = digits)[1]
    txt2 <- format(c(r$p.value, 0.123456789), digits = digits)[1]
    text(0.5, 0.5, paste(txt1, "\n", txt2, sep=""), cex = 1)
}


## ---- H2

H2=function(Y, K, id){
    data=data.frame(Y=Y, id=id)
    model=lmekin(Y~ (1|id) ,data, varlist=K, method="REML")
    return(H2lmekin(model, g="id"))
}

## ---- end-of-H2

H2nokin=function(Y,id,cov=NULL){
    if(is.null(cov)){
        data=data.frame(Y=Y, id=id)
        m=lmer(Y~ (1|id) ,data)}else{
            data=data.frame(Y=Y, id=id, cov=cov)
            m=lmer(Y~ cov + (1|id) ,data)}
    return(var_ID(m))}


H2bootNokin=function(Y, id, cov=null, nsim=1000, ncpus=4){
    if(is.null(cov)){
        data=data.frame(Y=Y, id=id)
        m=lmer(Y~ (1|id) ,data)
        b=bootMer(m, nsim=nsim, ncpus=ncpus, FUN=var_ID, parallel="snow")}else{
            data=data.frame(Y=Y, id=id, cov=cov)
            m=lmer(Y~ cov+(1|id) ,data)
            b=bootMer(m, nsim=nsim, ncpus=ncpus, FUN=var_ID, parallel="snow")}
    return(b$t)
}
