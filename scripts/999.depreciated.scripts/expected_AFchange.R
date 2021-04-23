library(snpStats)
library(plyr)

##read teh phenotypic data from 2011 and 2012

d11=readRDS(file="./data/d11_for_analysis.rds")
d12=readRDS(file="./data/d12_for_analysis.rds")

d11=droplevels(d11[d11$id!=".",])
d12=droplevels(d12[d12$id!=".",])

##compute accession level values for ows 2011, fec 2011 and ows 2012, in each region
South=c("ULL", "RAT")
North=c("RAM", "ADA")

regions=list(South=South, North=North)


acc=read.table("./data/acc_list.txt", h=T, sep="\t")
res=data.frame(id=paste(acc[1:100,1]))

##ows 2011
for(r in names(regions)){
    subd11=d11[d11$exp%in%regions[[r]], c("exp", "block", "id","ows", "fecundity")]
    ##compute binomial probability of survival per accession
    m1=glm(ows~exp*id, data=subd11[is.na(subd11$ows)==F,], family = binomial)
    ##get fitted values and average per accession
    f=fitted.values(m1)
    d=data.frame(id=subd11$id[is.na(subd11$ows)==F], f)
    avd=ddply(d, "id", function(x){mean(x[,2])})
    ##make relative
    avd[,2]=avd[,2]/max(avd[,2])
    row.names(avd)=paste(avd[,1])
    res[,paste("surv_2011_", r, sep="")]=avd[paste(res$id),2]
}

##fecundity 2011

for(r in names(regions)){
    subd11=d11[d11$exp%in%regions[[r]], c("exp", "block", "id","ows", "fecundity")]
    ##compute binomial probability of survival per accession
    m2=glm(fecundity~exp*id, data=subd11[is.na(subd11$ows)==F,], family = gaussian)
    ##get fitted values and average per accession
    f=fitted.values(m2)
    d=data.frame(id=subd11$id[is.na(subd11$fecundity)==F], f)
    avd=ddply(d, "id", function(x){mean(x[,2])})
     ##make relative
    avd[,2]=avd[,2]/max(avd[,2])
    row.names(avd)=paste(avd[,1])
    res[,paste("fec_2011_", r, sep="")]=avd[paste(res$id),2]
}

##survival in 2012

for(r in names(regions)){
    subd12=d12[d12$exp%in%regions[[r]], c("exp", "block", "id","ows", "fecundity")]
    ##compute binomial probability of survival per accession
    m3=glm(ows~exp*id, data=subd12[is.na(subd12$ows)==F,], family = binomial)
    ##get fitted values and average per accession
    f=fitted.values(m3)
    d=data.frame(id=subd12$id[is.na(subd12$ows)==F], f)
    avd=ddply(d, "id", function(x){mean(x[,2])})
    ##make relative
    avd[,2]=avd[,2]/max(avd[,2])
    row.names(avd)=paste(avd[,1])
    res[,paste("surv_2012_", r, sep="")]=avd[paste(res$id),2]
}

saveRDS(res, "./res/CG_pred_ows_fec.rds")

##read in the SNP matrix

bed="./GWA/snps/sweden_200_MAF10.bed"
fam="./GWA/snps/sweden_200_MAF10.fam"
bim="./GWA/snps/sweden_200_MAF10.bim"
snps=read.plink(bed, bim, fam)

##choose a random SNP with reasonnable frequency to run tests

w=152563
x=as.numeric(matrix(snps$genotypes[,w]))-1
names(x)=snps$fam[,1]
afi=sum(x, na.rm=T)/(length(x)*2)

x[x==1]=NA
## N dispersed=120 per accession
Ndisp=120

pred=data.frame(id=res[,1], South=NA, North=NA)

for(r in c("South", "North")){
    ##compute new SNP frequencies 
    S=res[, c(1, grep(r, colnames(res)))]
                                        #S$snp=x[paste(S[,1])]
    S=na.omit(S)
    ##filter of first winter. 
    a=S[,paste("surv_2011_", r, sep="")] ##relative abundance after the first winter
    ##seed production * survival
    b=a*S[,paste("fec_2011_", r, sep="")]
    b=b/max(b)
    ##survival second winter
    d=b*S[,paste("surv_2012_", r, sep="")]
    d=d/max(d)
    names(d)=paste(S[,1])
    pred[,r]=d[paste(pred$id)]
}

##now compute final SNP frequencies per region.

snp=x[paste(pred[,1])]
r="North"

##pop=paste(unlist(
pred=na.omit(pred)

popS=paste(unlist(dlply(pred[,c("id", "South")], "id", function(x){return(rep(x[,1],round(as.numeric(x[,2]), 2)*10))})))
popN=paste(unlist(dlply(pred[,c("id", "North")], "id", function(x){return(rep(x[,1],round(as.numeric(x[,2]), 2)*10))})))

##compute frequencies from the SNP matrix

ifreq=apply(snps$genotypes, 2, function(x){(sum(x==03)+sum(x==02))/(length(x)*2)})
saveRDS(ifreq, "./res/initial_freq_200acc.rds")

snpS=snps$genotypes[popS,]
predSouth=apply(snpS, 2, function(x, L){(sum(x=="03")+sum(x=="02")/2)/(L)}, L=2*length(popS))
saveRDS(predSouth, "./res/pred_AF_South.rds")

snpN=snps$genotypes[popN,]
predNorth=apply(snpN, 2, function(x, L){(sum(x=="03")+sum(x=="02")/2)/(L)}, L=2*length(popN))
saveRDS(predNorth, "./res/pred_AF_North.rds")


##compare with real AF data

predSouth=readRDS("./res/pred_AF_South.rds")
predNorth=readRDS("./res/pred_AF_North.rds")

load("./007.calculate.afd.from.start/pol.AFD.out.Rdata")
af=pol.AFD.out
rm(pol.AFD.out)

af$snp=paste(af$Chromosome, "_",af$Position, sep="")
af$predS=predSouth[paste(af$snp)]
af$predN=predNorth[paste(af$snp)]

cor.test(af$N.mean, af$predN-ifreq[paste(af$snp)])
cor.test(af$S.mean, af$predN-ifreq[paste(af$snp)])

png("./figures/comp.png")
plot(af$S.mean, af$predS-ifreq[paste(af$snp)], cex=0.5, pch=16)
dev.off()
