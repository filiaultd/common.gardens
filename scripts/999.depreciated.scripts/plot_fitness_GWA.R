library(readr)

fit=read.table("./res/means_NvsS.txt", h=T, sep="\t")

acclist=read.table("./data/acc_list.txt", h=T, sep="\t")

##read in the ancestral state of SNPs
## g=gzfile("./data/snpeff.csv.gz.gz")
## anc=read_delim(g, col_names=F, delim="\t", na=c("NA", "<NA>"))
## anc=anc[,1:5]
## anc=as.data.frame(anc)
## colnames(anc)=c("chr", "pos", "ref", "alt", "anc")
## ##anc=anc[is.na(anc$anc)==F,]
## anc$chr=gsub("Chr", "", anc$chr)
## anc$rs=paste(anc$chr, "_", anc$pos, sep="")

## saveRDS(anc, "./data/anc.rds")

#anc=readRDS("./data/anc.rds")

##read in the genotypes for the 200 swedes

bedpath="./GWA/snps/sweden_200_MAF10"
require(snpStats)
require(reshape)
#require(fpc)
fam <- paste(bedpath, ".fam", sep="")
bim <- paste(bedpath, ".bim", sep="")
bed <- paste(bedpath, ".bed", sep="")
snps=read.plink(bed, bim, fam)

map=snps$map

colnames(map)[2]="rs"

## m=merge(map, anc, all.x=T, all.y=F, by="rs", sort=F)
## m=m[,-3] ##remove the centimorgan col

## ##how many SNPs have ancestry info

## table(m$anc)
## sum(is.na(m$anc)) ##we're missing anc info for 327053
## sum(is.na(m$anc)==F)/nrow(m) ##we have ancestry info for 57%

## ##how many of our SNPs were not in the snpeff file

## sum(map$rs%in%anc$rs==F)/nrow(map) ##6.8% of our SNPs are not in the SNPeff file so redoing it wouldn"t improve ancestry assignation very much.

## ## Is the reference allele always the same in the two datasets

## ##among those do ref alleles match
## sum(m$allele.2!=m$ref, na.rm=T)/nrow(m[is.na(m$pos)==F,]) ##na rm removes the SNPs that were not in the SNPeff file from this calculation

## ##my GWA effects are calculated for the rare vs common allele, not from ref/alt.
## ##columns refereing the alleles for the GWAs are allele.1 and allele.2. Columns refering to the SNPeff file are ref alt. We can deal with ancestrality later

## m$flipeff=NA
## m$flipeff[m$allele.2==m$alt]=-1
## m$flipeff[m$allele.2==m$ref]=1

##remove rows of m with is.na(flipeff)==T
##for these SNPs, I currently not tell which one is the reference

##m=m[is.na(m$flipeff)==F,]

##read the the GWA results

file=file("./GWA/means_NvsS/output/mlmm_fit_NvsS.assoc.txt")
gwa=read_delim(file, delim="\t",col_names=T)
##remove SNPs that are not in m
##gwa=gwa[gwa$rs%in%m$rs,]
##make sure m and GWA are in the same order
##sum(gwa$rs==m$rs)==nrow(gwa)
##adjust effect to make it the effect of the alternate allele
##gwa$beta_1=gwa$beta_1*m$flipeff
##gwa$beta_2=gwa$beta_2*m$flipeff

gwa$score=-log10(gwa$p_wald)
##subset to top snps
th=quantile(gwa$score, 0.99)
gwa=gwa[gwa$score>=th,]


##explore the table of Col-0 effects in different quadrants of the previous graph
##for each SNP, define its effet in both region
gwa$eff=NA
gwa$eff[gwa$beta_1>0 & gwa$beta_2>0]="PP"
gwa$eff[gwa$beta_1>0 & gwa$beta_2<0]="PN"
gwa$eff[gwa$beta_1<0 & gwa$beta_2>0]="NP"
gwa$eff[gwa$beta_1<0 & gwa$beta_2<0]="NN"

table(gwa$eff)

##make a figures
fam=merge(snps$fam[,1:2],acclist, by.x="pedigree", by.y="lines", keep.x=T, keep.y=F, sort=F)
subg=snps$genotypes[,map$rs%in%gwa$rs]

pdf("./figures/coeff_NvsS.pdf", paper="special", width=5, height=5, pointsize=8)
layout(matrix(c(2, 1, 1, 3, 1, 1, 1, 1, 4, 1, 1, 5), ncol=4, byrow=T))
plot(gwa$beta_2, gwa$beta_1, pch=16, col="Dodgerblue", cex=0.5, ylab="Effect of alternative on Fitness in the North", xlab="Effect of alternative allele on Fitness in the South", cex.lab=1.4)
abline(0,0)
segments(0, -100,0, +100)
legend("topleft", paste("Number of SNPs=", nrow(gwa),"out of", nrow(snps$map)), bty="n")
for(q in c("PN", "PP", "NN", "NP")){
    x=subg[,gwa$eff==q]
    Nx=x[fam$region=="N Sweden",]
    Sx=x[fam$region=="S Sweden",]
    #f=apply(x, 2, function(x){(sum(x=="01")+sum(x=="02")/2)/length(x)})
    Nf=apply(Nx, 2, function(x){(sum(x=="01")+sum(x=="02")/2)/length(x)})
    Sf=apply(Sx, 2, function(x){(sum(x=="01")+sum(x=="02")/2)/length(x)})
    ##determine which allele is Northern, which is Southern for each SNP and polarize allele frequencing
    R=rep("N", ncol(x))
    R[Nf<Sf]="S"
    par(mgp=c(-1.2, -1.2, -1))
    pie(table(R), xlab=paste("N=",sum(gwa$eff==q), sep=""), cex.lab=2)
                                        #fPN=f;fPN[Nf<Sf]=1-f[Nf<Sf]
    #fPS=f;fPS[Sf<Nf]=1-f[Sf<Nf]
    #hist(fPN, ylim=c(-2000, 2000), col="Dodgerblue")
    #par(new = TRUE)
    #hist(fPS, ylim=c(2000, -2000), col="Firebrick")
    #par(new = FALSE)
#plot(dN, col="Dodgerblue", main=""); points(dS, type="l", col="Firebrick3")
}
dev.off()

##compare GWA coeffs with changes in AF in SR allelic freq

load("./data/ns.AFD.out.Rdata")
SR=ns.AFD.out
SR$rs=paste(SR$Chromosome, "_", SR$Position, sep="")
NSR=nrow(SR)

##remove changing allele freq where m<sd

SR=SR[(abs(SR$N.mean)>(10*SR$N.sd) & abs(SR$N.mean)>0.1) | (abs(SR$S.mean)>(10*SR$S.sd) & abs(SR$S.mean)>0.1),]

sub=SR[SR$rs%in%gwa$rs,]
SR$eff=NA
SR$eff[SR$N.mean>0 & SR$S.mean>0]="PP"
SR$eff[SR$N.mean<0 & SR$S.mean>0]="NP"
SR$eff[SR$N.mean>0 & SR$S.mean<0]="PN"
SR$eff[SR$N.mean<0 & SR$S.mean<0]="NN"

##make a figures
fam=merge(snps$fam[,1:2],acclist, by.x="pedigree", by.y="lines", keep.x=T, keep.y=F, sort=F)
subg=snps$genotypes[,map$rs%in%SR$rs]

pdf("./figures/AF_NvsS.pdf", paper="special", width=5, height=5, pointsize=8)
layout(matrix(c(2, 1, 1, 3, 1, 1, 1, 1, 4, 1, 1, 5), ncol=4, byrow=T))
plot(SR$S.mean,SR$N.mean, pch=16, col="Dodgerblue", cex=0.5, ylab="Change alternative allele in the North", xlab="Change alternative allele in the South", xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4))
points(sub$N.mean, sub$S.mean, pch=16, col="Gold2", cex=1)
abline(0,0, lty=2)
segments(0, -100,0, +100, lty=2)
legend("topleft", paste("Number of significant SNPs in SR data=", nrow(SR),"out of", NSR, "; Number of SNPs also in top gwa hits=", nrow(sub), "out of", nrow(gwa), "top gwa hits"), bty="n")
for(q in c("PN", "PP", "NN", "NP")){
    x=subg[,SR$eff==q]
    Nx=x[fam$region=="N Sweden",]
    Sx=x[fam$region=="S Sweden",]
    #f=apply(x, 2, function(x){(sum(x=="01")+sum(x=="02")/2)/length(x)})
    Nf=apply(Nx, 2, function(x){(sum(x=="01")+sum(x=="02")/2)/length(x)})
    Sf=apply(Sx, 2, function(x){(sum(x=="01")+sum(x=="02")/2)/length(x)})
    ##determine which allele is Northern, which is Southern for each SNP and polarize allele frequencing
    R=rep("N", ncol(x))
    R[Nf<Sf]="S"
    par(mgp=c(-1.2, -1.2, -1))
    pie(table(R), xlab=paste("N=",sum(SR$eff==q), sep=""), cex.lab=2)
}
dev.off()

##merge the results from gwa and SR (only keeping top GWA hits)

B=merge(gwa, SR[,c("rs", "N.mean", "N.sd", "S.mean", "S.sd")], by="rs")

##make figure of consistency between GWA an AF change


pdf("./figures/gwa_vs_AF.pdf", paper="special", height=5, width=5, pointsize=8)
layout(matrix(1:4, ncol=2, byrow=T))
par(ma r=c(4, 4, 1, 1))
plot(B$beta_1, B$N.mean,cex=0.5, pch=16, col="Dodgerblue", xlab="fecundity effect in the North",ylab="Allele frequency change in the North")
m=lm(B$N.mean~B$beta_1)
abline(m)
C=cor.test(B$beta_1, B$N.mean, method="pearson")
legend("topright", paste("r=",format(C$estimate, digits=3), "\np-value=", format(C$p.value, digits=3, scientific=T), "\nNumber of SNPs considered=",nrow(B)), bty="n")
m=lm(B$N.mean~B$beta_1)
abline(m)
plot(B$beta_1, B$S.mean,cex=0.5, pch=16, col="Dodgerblue", xlab="fecundity effect in the North",ylab="Allele frequency change in the South")
C=cor.test(B$beta_1, B$S.mean, method="pearson")
legend("topright", paste("r=",format(C$estimate, digits=3), "\np-value=", format(C$p.value, digits=3, scientific=T)), bty="n")
m=lm(B$S.mean~B$beta_1)
abline(m)
plot(B$beta_2, B$N.mean,cex=0.5, pch=16, col="Dodgerblue", xlab="fecundity effect in the South", ylab="Allele frequency change in the North")
C=cor.test(B$beta_2, B$N.mean, method="pearson")
legend("topright", paste("r=",format(C$estimate, digits=3), "\np-value=", format(C$p.value, digits=3, scientific=T)), bty="n")
m=lm(B$N.mean~B$beta_2)
abline(m)
plot(B$beta_2, B$S.mean,cex=0.5, pch=16, col="Dodgerblue", xlab="fecundity effect in the South",ylab="Allele frequency change in the South")
C=cor.test(B$beta_2, B$S.mean, method="pearson")
legend("topright", paste("r=",format(C$estimate, digits=3), "\np-value=", format(C$p.value, digits=3, scientific=T)), bty="n")
m=lm(B$S.mean~B$beta_2)
abline(m)
dev.off()

##quantify how much SNP overlap between GWA and SR

##re-readin the GWA results

file=file("./GWA/means_NvsS/output/mlmm_fit_NvsS.assoc.txt")
gwa=read_delim(file, delim="\t",col_names=T)
##make sure m and GWA are in the same order
gwa$score=-log10(gwa$p_wald)



##re-readin the SR data
load("./data/ns.AFD.out.Rdata")
SR=ns.AFD.out
SR$rs=paste(SR$Chromosome, "_", SR$Position, sep="")
NSR=nrow(SR)

##merge the two data sets and keep everything

B=merge(gwa, SR[,c("rs", "N.mean", "N.sd", "S.mean", "S.sd")], by="rs", all=F)
##subset to top snps

th=quantile(B$score, 0.99)

SR_signif=(abs(SR$N.mean)>5*SR$N.sd & abs(SR$N.mean)>0.1) | (abs(SR$S.mean)>5*SR$S.sd & abs(SR$S.mean)>0.1)
#    (abs(B$N.mean)>5*B$N.sd | abs(B$S.mean)>5*B$S.sd) & (abs(B$N.mean)>0.1 | abs(B$S.mean))

GWA_signif=B$score>=th

tab=table(SR_signif, GWA_signif)

C=chisq.test(tab)



SR=SR[abs(SR$N.mean)>2*SR$N.sd | abs(SR$S.mean)>2*SR$S.sd,]
sub=SR[SR$rs%in%gwa$rs,]


##subset to top snps
th=quantile(gwa$score, 0.95)
gwa=gwa[gwa$score>=th,]


###### Daniele - output raw gwa variable for analysis with AFDs.
gwa.out <- gwa
colnames(gwa.out)[8:9] <- c("beta_south", "beta_north")
gwa.out <- as.data.frame(gwa.out)
save(gwa.out, file="gwa.out.Rdata")
