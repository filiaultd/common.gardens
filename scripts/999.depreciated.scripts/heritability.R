
## ---- herit
library(lme4)
library(HLMdiag)
library(boot)
library(gplots)

source("./scripts/heritability_functions.R")

##Variance explained by ID within

d11=readRDS(file="./data/d11_for_analysis.rds")
d12=readRDS(file="./data/d12_for_analysis.rds")

acclist=read.table("./data/acc_list.txt", h=T, sep="\t")
acclist=acclist[acclist$tubes<=200,]
###histograms of trait distributions.
##A0013##
##d11
pdf("./figures/dist_d11.pdf", paper="special", height=8, width=8)
par(mfrow=c(3, 4), mar=c(5, 5, 1, 1))
for (t in c("area", "perimeter", "max_diameter", "sdR", "circle_area", "stockiness", "color", "herbivory", "fecundity", "ows", "sss", "fitness")){
    hist(d11[,t], breaks=50, main=t, xlab="")
}
dev.off()

##d12
pdf("./figures/dist_d12.pdf", paper="special", height=8, width=8)
par(mfrow=c(3, 4), mar=c(5, 5, 1, 1))
for (t in c("area", "perimeter", "max_diameter", "sdR", "circle_area", "stockiness", "color", "FT", "fecundity", "ows", "sss", "fitness")){
    hist(d12[,t], breaks=50, main=t, xlab="")
}
dev.off()

##Heritability estimates by exp and by year.

###traits=c("area", "perimeter", "max_diameter", "sdR", "circle_area", "stockiness", "color", "herbivory", "FT", "fecundity", "ows", "sss", "fitness")

traits=c("area", "stockiness", "color", "FT", "fecundity", "ows", "sss", "fitness")

##set the family distribution for each trait.

fam=data.frame(traits, family=c("log", "log", "log", "log", "sqrt", "bn", "bn", "sqrt"))
exps=sort(c("ULL", "RAT", "ADA", "RAM"))

##loop over years, and experiments and fill a dataframe with the H2 and CI

res=data.frame(matrix(ncol=6))
colnames(res)=c("year", "experiment" ,"trait", "H2", "CI_low", "CI_high")
m=data.frame(id=unique(c(paste(d11$id), paste(d12$id))))
l=1
blups=data.frame(ID=acclist[,1])
for(year in c(2011, 2012)){
    for(e in exps){
        if(year==2011){
            sub=droplevels(d11[d11$exp==e,])
        }
        if(year==2012){
            sub=droplevels(d12[d12$exp==e,])
        }
        s=sub[,c("exp", "block","id", traits[traits%in%colnames(sub)]),]
        ##remove columns with just NAs
        s=s[,apply(s, 2, function(x){any(is.na(x)==F)})]
        ##remove herbivory for now because I can't make it fit.
        s=s[,colnames(s)!="herbivory"]
        ##fill in res
        res[l:(l+(ncol(s)-4)),1]=paste(year)
        res[l:(l+(ncol(s)-4)),2]=e
        res[l:(l+(ncol(s)-4)),3]=colnames(s)[4:ncol(s)]
        ##models and H2 estimates
        for(t in colnames(s)[4:ncol(s)]){
            f=fam[match(t, fam$traits), 2]
            if(f=="log"){
                model=formula(paste("log(", t, "+1)~block + (1|id)", sep=""))
                fit=lmer(model, data=s)}
            if(f=="g"){
                model=formula(paste(t, "~block + (1|id)", sep=""))
                fit=lmer(model, data=s)}
            if(f=="sqrt"){
                model=formula(paste("sqrt(", t, ")~block + (1|id)", sep=""))
                fit=lmer(model, data=s)}
            if(f=="bn"){
                model=formula(paste(t,"~block + (1|id)", sep=""))
                fit=glmer(model, data=s, family="binomial")}
            if(f=="nb"){
                model=formula(paste(t,"~block + (1|id)", sep=""))
                fit=glmer.nb(model, data=s)}
            H2=var_ID(fit);ci=get_CI(fit, var_ID, nsim=100)
            res[l,4:6]=c(H2, ci[1], ci[2]); l=l+1
            b=ranef(fit)$id; b=data.frame(id=row.names(b), t=b[,1])
            blups[,paste(year, "_",e, "_", t, sep="")]=b[match(blups[,1], b[,1]),2]  ##A0018##
        }
    }
}

write.table(res, "./res/H2_2011_2012.txt", sep="|",col.names=T, row.names=F, quote=F)
write.table(blups, "./res/blups.txt", sep="\t",col.names=T, row.names=F, quote=F)

res$H2=round(res$H2, 3)
res$CI_low=round(res$CI_low, 3)
res$CI_high=round(res$CI_high, 3)

##making a figure.

res=read.table("./res/H2_2011_2012.txt", sep="|",h=T)
res=res[order(res[,3]),]
res=res[order(res[,2]),]
res=res[order(res[,1]),]

cols=c("dodgerblue1", "dodgerblue4", "gold3", "firebrick4")[as.numeric(as.factor(res$experiment))]
pdf("./figures/H2.pdf", paper="special", width=12, height=5)
par(mar=c(8, 5, 1, 1))
plot(0, 0, xaxt="n", type="n", xlim=c(0.5, nrow(res)+0.5), ylim=c(0, 1), ylab="H2", xlab="")
for (i in 1:nrow(res)){
    rect(i-0.5,0,i+0.5, res[i, 4], col=cols[i])
    segments(i,res[i, 5],i, res[i, 6], col=1)
}
axis(1, labels=res[,3], at=1:nrow(res), las=2, cex.lab=0.6)
legend("topright", legend=levels(as.factor(res$experiment)), col=c("dodgerblue1", "dodgerblue4", "gold3", "firebrick4"), pch=15)
segments(which(res[,1]==2012)[1]-0.5, 0,which(res[,1]==2012)[1]-0.5, 1)
text((which(res[,1]==2012)[1]-0.5)/2, 0.9, "2011", cex=3)
text(which(res[,1]==2012)[1]-0.5+(nrow(res)-which(res[,1]==2012)[1]-0.5)/2 , 0.9, "2012", cex=3)
dev.off()

## ---- end-of-herit

## ---- heritfit

##H2 estimate for just fitness

exps=c("ULL", "RAT", "ADA", "RAM")
res=data.frame(matrix(ncol=6))
colnames(res)=c("year", "experiment" ,"trait", "H2", "CI_low", "CI_high")
m=data.frame(id=unique(c(paste(d11$id), paste(d12$id))))
l=1
for(year in c(2011, 2012)){
    for(e in exps){
        if(year==2011){
            sub=droplevels(d11[d11$exp==e,])
        }
        if(year==2012){
            sub=droplevels(d12[d12$exp==e,])
        }
        s=sub[,c("exp", "block","id", "fitness")]##traits[traits%in%colnames(sub)]),]
        ##fill in res
        res[l:(l+(ncol(s)-4)),1]=paste(year)
        res[l:(l+(ncol(s)-4)),2]=e
        res[l:(l+(ncol(s)-4)),3]=colnames(s)[4:ncol(s)]
        ##models and H2 estimates
        model=formula(paste("log(fitness+1)~block + (1|id)", sep=""))
        fit=lmer(model, data=s)
        H2=var_ID(fit);ci=get_CI(fit, var_ID, nsim=1000)
        res[l,4:6]=c(H2, ci[1], ci[2]); l=l+1
    }
}

write.table(res, "./res/H2_fitness.txt", sep="|",col.names=T, row.names=F, quote=F)
res$H2=round(res$H2*100, 2)
res$CI_low=round(res$CI_low*100, 2)
res$CI_high=round(res$CI_high*100, 2)
write.table(res, "./res/H2_fitness_percent_rounded.txt", sep="|",col.names=T, row.names=F, quote=F)

## ---- end-of-heritfit


## ---- h2fitNS
##calculate residual means per after regressing the block effect.

source("./scripts/custom_functions.R")

d11=readRDS(file="./data/d11_for_analysis.rds")
d12=readRDS(file="./data/d12_for_analysis.rds")

##restric to fitness
traits=c("fitness")

d11$year=2011
d12$year=2012

d11=d11[, c("year", "exp", "block", "id", "fitness")]
d12=d12[, c("year", "exp", "block", "id", "fitness")]

d=rbind(d11, d12)
d$region[d$exp%in%c("ULL", "RAT")]="S"
d$region[d$exp%in%c("RAM", "ADA")]="N"
d$year=as.factor(d$year)
d=d[, c("year", "region", "exp","block", "id", "fitness")]
d=na.omit(d)
##recode block

###within region, remove exp effect and block within ##A0022##

m=data.frame(id=unique(c(paste(d11$id), paste(d12$id))))
l=1
for (R in c("N", "S")){
    dat=d[d$region==R,]
    model=formula(paste("fitness~0+block%in%exp%in%year", sep=""))
    fit=glm(model, data=dat)
    r=fit$residuals+mean(fit$coefficients) ##here I the mean effect of the experimental units to keep a scale to compare the two regions
    b=by(r, dat$id[is.na(dat[,"fitness"])==F], betterMean, N=3)
    b=data.frame(id=names(b), as.numeric(b))
    m[,paste(R, "_fitness", sep="")]=b[match(m[,1], b[,1]),2]
    l=l+1
}

write.table(m, "./res/means_NvsS.txt", sep="\t", col.names=T, row.names=F, quote=F)

## ---- end-of-h2fitNS
