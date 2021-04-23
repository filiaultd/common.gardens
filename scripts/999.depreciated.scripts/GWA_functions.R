##functions

##mahattan plots

manhattan=function(pval, yaxis="pval", layout=T){
    #pval$chr=as.numeric(substring(pval$snp, 1, 1))
    #pval$pos=as.numeric(gsub("[1-5][:-:_: :]+", "", pvals[,1]))
    for (c in 1:5){
        assign(paste("pv",c,sep=""),pval[pval$chr==c,c(1:ncol(pval))])
    }
    pv1$manhattan=pv1[,"pos"]
    pv2$manhattan=pv2[,"pos"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"pos"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"pos"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"pos"]+max(pv4$manhattan)
    pv1$col="#133A7F"
    pv2$col="#518FFF"
    pv3$col="#2774FF"
    pv4$col="#4F607F"
    pv5$col="#1F5CCC"
    pv=pval
    pv$manhattan=c(pv1$manhattan,pv2$manhattan,pv3$manhattan,pv4$manhattan,pv5$manhattan)
    pv$col=c(pv1$col,pv2$col,pv3$col,pv4$col,pv5$col)
    deb=min(pv$manhattan)
    fin=max(pv$manhattan)
    if(layout==T){
        m=matrix(c(2, 3, 0, 1), ncol=2, byrow=T)
        layout(m,c(2,30),c(10,2))}
    par(mar=c(0,0.7,0,0.5))
    plot(0,0,xlim=c(deb,fin),ylim=c(0,2),type="n",axes=FALSE,frame.plot=FALSE, xaxs="i", yaxs="i")
    for(c in 1:5){
        X=min(pv$manhattan[pv[,"chr"]==c])+(max(pv$manhattan[pv[,"chr"]==c])-min(pv$manhattan[pv[,"chr"]==c]))/2
        text(X,1,paste("Chromosome",c),cex=1, offset=0)
    }
    par(mar=c(0,0,0,0), ps=8)
    plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
    if(yaxis=="pval"){
        text(0.06,0.5,bquote(-log^10~(italic(p)-value)),srt=90,cex=1.5)
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        ymax=max(-log10(pv[,"pval"]), na.rm=T)[1]+1
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        plot(pv$manhattan,-log10(pv[,"pval"]),xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(1,ymax),axes=FALSE, cex =0.5 ,  col =pv$col ,pch=16, xaxs="i", yaxs="i")}
    if(yaxis=="adjr2"){
        text(0.06,0.5,bquote(adjusted ~ italic(r^2)),srt=90,cex=1.5)
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        ymax=max(pv[,"adjr2"], na.rm=T)+0.01
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        plot(pv$manhattan,pv[,"adjr2"],xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(0.01,ymax),axes=FALSE, cex =0.5 ,  col =pv$col ,pch=16, xaxs="i", yaxs="i")}
    axis(1,cex.axis=1,lwd=0.5,xaxp=c(deb,fin,10),labels=FALSE, tck=0)
    axis(2,cex.axis=1,lwd=0.5,tck=0.03)
    box(lwd=0.5)
}


##return just manhattan positions
manhattan_pos=function(pval){
    #pval$chr=as.numeric(substring(pval$snp, 1, 1))
                                        #pval$pos=as.numeric(gsub("[1-5][:-:_: :]+", "", pvals[,1]))
    for (c in 1:5){
        assign(paste("pv",c,sep=""),pval[pval$chr==c,c(1:ncol(pval))])
    }
    pv1$manhattan=pv1[,"pos"]
    pv2$manhattan=pv2[,"pos"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"pos"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"pos"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"pos"]+max(pv4$manhattan)
    pv1$col="#133A7F"
    pv2$col="#518FFF"
    pv3$col="#2774FF"
    pv4$col="#4F607F"
    pv5$col="#1F5CCC"
    pv=pval
    pv$manhattan=c(pv1$manhattan,pv2$manhattan,pv3$manhattan,pv4$manhattan,pv5$manhattan)
    return(pv$manhattan)
}


##manhattan with highlighted regions

manhattan2=function(pval, yaxis="pval", layout=T,  reg=NULL, winreg=20000){
    #pval$chr=as.numeric(substring(pval$snp, 1, 1))
    #pval$pos=as.numeric(gsub("[1-5][:-:_: :]+", "", pvals[,1]))
    for (c in 1:5){
        assign(paste("pv",c,sep=""),pval[pval$chr==c,c(1:ncol(pval))])
    }
    if(is.null(reg)==F){
         for (c in 1:5){
             assign(paste("reg",c,sep=""),reg[reg$chr==c,c(1:ncol(reg))])
    }}        
    pv1$manhattan=pv1[,"pos"]
    pv2$manhattan=pv2[,"pos"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"pos"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"pos"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"pos"]+max(pv4$manhattan)
    if(is.null(reg)==F){
        reg1$m1=reg1$pos-winreg
        reg1$m2=reg1$pos+winreg
        reg2$m1=reg2$pos-winreg+max(pv1$manhattan)
        reg3$m1=reg3$pos-winreg+max(pv2$manhattan)
        reg4$m1=reg4$pos-winreg+max(pv3$manhattan)
        reg5$m1=reg5$pos-winreg+max(pv4$manhattan)
        reg2$m2=reg2$pos+winreg+max(pv1$manhattan)
        reg3$m2=reg3$pos+winreg+max(pv2$manhattan)
        reg4$m2=reg4$pos+winreg+max(pv3$manhattan)
        reg5$m2=reg5$pos+winreg+max(pv4$manhattan)
    }
    pv1$col="#133A7F"
    pv2$col="#518FFF"
    pv3$col="#2774FF"
    pv4$col="#4F607F"
    pv5$col="#1F5CCC"
    pv=pval
    pv$manhattan=c(pv1$manhattan,pv2$manhattan,pv3$manhattan,pv4$manhattan,pv5$manhattan)
    reg=rbind(reg1, reg2, reg3, reg4, reg5)
    pv$col=c(pv1$col,pv2$col,pv3$col,pv4$col,pv5$col)
    deb=min(pv$manhattan)
    fin=max(pv$manhattan)
    if(layout==T){
        m=matrix(c(2, 3, 0, 1), ncol=2, byrow=T)
        layout(m,c(2,30),c(10,2))}
    par(mar=c(0,0.7,0,0.5))
    plot(0,0,xlim=c(deb,fin),ylim=c(0,2),type="n",axes=FALSE,frame.plot=FALSE, xaxs="i", yaxs="i")
    for(c in 1:5){
        X=min(pv$manhattan[pv[,"chr"]==c])+(max(pv$manhattan[pv[,"chr"]==c])-min(pv$manhattan[pv[,"chr"]==c]))/2
        text(X,1,paste("Chromosome",c),cex=1, offset=0)
    }
    par(mar=c(0,0,0,0), ps=8)
    plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
    if(yaxis=="pval"){
        text(0.06,0.5,bquote(-log^10~(italic(p)-value)),srt=90,cex=1.5)
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        ymax=max(-log10(pv[,"pval"]), na.rm=T)[1]+1
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        plot(0,0,xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(1,ymax),axes=FALSE, xaxs="i", yaxs="i", type="n")
        if(is.null(reg)==F){
            rect(reg$m1, rep(0, nrow(reg)), reg$m2, rep(50, nrow(reg)), col="gold",lwd=0)} 
        points(pv$manhattan,-log10(pv[,"pval"]), cex =0.5 ,  col =pv$col ,pch=16)}
    if(yaxis=="adjr2"){
        text(0.06,0.5,bquote(adjusted ~ italic(r^2)),srt=90,cex=1.5)
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        ymax=max(pv[,"adjr2"], na.rm=T)+0.01
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        plot(0,0,xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(1,ymax),axes=FALSE, xaxs="i", yaxs="i", type="n")
        if(is.null(reg)==F){
            rect(reg$m1, rep(0, nrow(reg)), reg$m2, rep(50, nrow(reg)), col="gold",lwd=0)} 
        points(pv$manhattan,pv[,"adjr2"], cex =0.5 ,  col =pv$col ,pch=16)}
    axis(1,cex.axis=1,lwd=0.5,xaxp=c(deb,fin,10),labels=FALSE, tck=0)
    axis(2,cex.axis=1,lwd=0.5,tck=0.03)
    box(lwd=0.5)
}


##manhattan bayesR
manhattan_bayesR=function(asso, K=NULL){
    for (c in 1:5){
        assign(paste("pv",c,sep=""),asso[asso$chr==c,c(1:ncol(asso))])
    }
    pv1$manhattan=pv1[,"pos"]
    pv2$manhattan=pv2[,"pos"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"pos"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"pos"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"pos"]+max(pv4$manhattan)
    pv1$col="#133A7F"
    pv2$col="#518FFF"
    pv3$col="#2774FF"
    pv4$col="#4F607F"
    pv5$col="#1F5CCC"
    pv=asso
    pv$manhattan=c(pv1$manhattan,pv2$manhattan,pv3$manhattan,pv4$manhattan,pv5$manhattan)
    pv$col=c(pv1$col,pv2$col,pv3$col,pv4$col,pv5$col)
    deb=min(pv$manhattan)
    fin=max(pv$manhattan)
    m=matrix(c(2, 3, 0, 1), ncol=2, byrow=T)
    layout(m,c(2,30),c(10,2))
    par(mar=c(0,0.7,0,0.5))
    plot(0,0,xlim=c(deb,fin),ylim=c(0,2),type="n",axes=FALSE,frame.plot=FALSE, xaxs="i", yaxs="i")
    for (c in 1:5){
        X=min(pv$manhattan[pv[,"chr"]==c])+((max(pv$manhattan[pv[,"chr"]==c])-min(pv$manhattan[pv[,"chr"]==c])))/2
        text(X,1,paste("Chromosome",c),cex=1, offset=0)
    }
    par(mar=c(0,0,0,0), ps=8)
    plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
    text(0.06,0.5,bquote(beta),srt=90,cex=1.5)
    par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
    ymax=max(abs(pv[,"beta"]), na.rm=T)
    par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
    if(is.null(K)==F){pv=pv[pv$k%in%K,]}
    plot(pv$manhattan,abs(pv[,"beta"]),xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(0,ymax),axes=FALSE, cex =0.5 ,  col =pv$col ,pch=16, xaxs="i", yaxs="i")
    axis(1,cex.axis=1,lwd=0.5,xaxp=c(deb,fin,10),labels=FALSE, tck=0)
    axis(2,cex.axis=1,lwd=0.5,tck=0.03)
    box(lwd=0.5)
}



##id manhattan
manhattanID=function(pval, yaxis="pval"){
    #pval$chr=as.numeric(substring(pval$snp, 1, 1))
    #pval$pos=as.numeric(gsub("[1-5][:-:_: :]+", "", pvals[,1]))
    require(Acinonyx)
    for (c in 1:5){
        assign(paste("pv",c,sep=""),pval[pval$chr==c,c(1:ncol(pval))])
    }
    pv1$manhattan=pv1[,"pos"]
    pv2$manhattan=pv2[,"pos"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"pos"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"pos"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"pos"]+max(pv4$manhattan)
    #pv1$col="#133A7F"
    #pv2$col="#518FFF"
    #pv3$col="#2774FF"
    #pv4$col="#4F607F"
    #pv5$col="#1F5CCC"
    pv=pval
    pv$manhattan=c(pv1$manhattan,pv2$manhattan,pv3$manhattan,pv4$manhattan,pv5$manhattan)
    #pv$col=c(pv1$col,pv2$col,pv3$col,pv4$col,pv5$col)
    deb=min(pv$manhattan)
    fin=max(pv$manhattan)
    if(yaxis=="adjr2"){
        iplot(pv$manhattan,pv[,"adjr2"])
        readline(prompt="Press [enter] to continue")
        return(pv[iset.selected(),])
    }
    if(yaxis=="pval"){
        iplot(pv$manhattan, pv[,"pval"])
        readline(prompt="Press [enter] to continue")
        return(pv[iset.selected(),])
    }
}


GWA_explorer2=function(pval_file, bedpath, figpath, threshold=5, win=10000, promoter=1){
    require(readr)
    require(GenomicRanges)
    require(TxDb.Athaliana.BioMart.plantsmart22)
    ##require(VariantAnnotation)
    require(snpStats)
    require(dbscan)
    require(reshape2)
    ##get genes from arabidopsis
    txdb=TxDb.Athaliana.BioMart.plantsmart22
    genes= genes(txdb)
    assoc=read_tsv(pval_file , col_names=T, progress=F)
    assoc=as.data.frame(assoc)
    colnames(assoc)[grep("p_", colnames(assoc))]="pval"
    assoc$plog=-log10(assoc$pval)
    #assoc=assoc[, c("chr", "rs", "ps", "af", "p_lrt","plog")]
    fam <- paste(bedpath, ".fam", sep="")
    bim <- paste(bedpath, ".bim", sep="")
    bed <- paste(bedpath, ".bed", sep="")
    snps=read.plink(bed, bim, fam)
      ##prepare a Genome wide track to be able to locate the peaks.
    for (c in 1:5){
        assign(paste("pv",c,sep=""),assoc[assoc$chr==c,c(1:ncol(assoc))])
    }
    pv1$manhattan=pv1[,"ps"]
    pv2$manhattan=pv2[,"ps"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"ps"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"ps"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"ps"]+max(pv4$manhattan)
    pv1$col="#133A7F"
    pv2$col="#518FFF"
    pv3$col="#2774FF"
    pv4$col="#4F607F"
    pv5$col="#1F5CCC"
    pv=rbind(pv1, pv2, pv3, pv4, pv5)
    deb=min(pv$manhattan)
    fin=max(pv$manhattan)
    assoc=pv
    top=assoc[assoc[,grep("plog", colnames(assoc))]>=threshold,grep("col", colnames(assoc), inv=T)]
    ##now get the candidate genes near top SNPs
    grgwa=with(top, GRanges(chr, IRanges(start=ps-win, end=ps+win)), pval=pval, score=y)
    ##look for overlap
    CG=subsetByOverlaps(genes, grgwa)
    ##for each chromosome, cluster top snps into peaks and plot
    for(chr in unique(top$chr)){
        topc=top[top$chr==chr,]
        if(nrow(topc)==1){grps=1}else{
            clust=dbscan(as.matrix(topc$ps), eps=2*win, minPts=1)
            grps=clust$cluster
        }
        for(g in unique(grps)){
            ug=topc[grps==g,]
            r=range(ug$ps)+c(-win, +win)
            r2=range(ug$ps)
            ##LD computation
            geno=snps$genotype[,snps$map$chromosome==chr & snps$map$position>=r[1] & snps$map$position<=r[2]]
            ##plotting
            jpeg(paste(figpath, "_", chr, "_",round(r[1]/1e6, 2), "_", round(r[2]/1e6, 2), "Mb.jpeg",sep=""),res=800, units="in", width=6, height=8)
            m=matrix(c(1:6, 6, 1:6, 7), ncol=2)
            layout(m, widths=c(5, 0.5), heights=c(2, 2,2,2, 1, 2, 3))
            par(mar=c(0, 4, 3, 2))
            plot(assoc$manhattan, assoc$plog, col=assoc$col, pch=16, cex=0.8, xlim=c(deb, fin), xaxt="n", xlab="",ylab="", type="p")
            points(top$manhattan, top$plog, col="gold3", pch=16, cex=1.2, xlim=c(deb, fin))
            legend("topleft", "A", bty="n")
            ##indicate the region zoomed on.
            rect(assoc$manhattan[assoc$chr==chr & assoc$ps==r2[1]], 0, assoc$manhattan[assoc$chr==chr & assoc$ps==r2[2]], max(abs(top$plog)), border="red", col=NULL, lwd=1)
            par(mar=c(0, 4, 3, 2))
            plot_genes(as.data.frame(CG), r=r)
            legend("topleft", "B", bty="n")
            axis(3)
            assocR=assoc[assoc$chr==chr & assoc$ps>=r[1] & assoc$ps<=r[2],]
            assocR$col[assocR$ps%in%top$ps]="gold3"
            assocR$lwd=0.2
            assocR$lwd[assocR$rs%in%top$ps]=3
            assocR$eff=-log10(assocR[, grep("pval", colnames(assocR))])
            par(mar=c(0,4, 0, 2))
            plot(assocR$ps, assocR$eff, xlim=r, ylim=range(assocR$eff), pch=16, xaxt="n", ylab=bquote(-log^10~(italic(p)-value)), xlab="", col=assocR$col)
            legend("topleft", "C", bty="n")
                                        #freq=col.summary(x)$MAF
            plot(assocR$ps, assocR$af, xlim=r, pch=16, xaxt="n", ylab="MAF", xlab="", col=assocR$col, yaxt="n", type="n")
            legend("bottomleft", "D", bty="n")
            segments(assocR$ps, 0 ,assocR$ps,  assocR$af, col=assocR$col, lwd=assocR$lwd)
            segments(assocR$ps[assocR$rs%in%top$rs], 0 ,assocR$ps[assocR$rs%in%top$rs], assocR$af[assocR$rs%in%top$rs], col="gold3", lwd=2)
            axis(2, at=c(0.1, 0.3), labels=c(0.1, 0.3))
            myPal <- colorRampPalette(c("deepskyblue4","firebrick3", "gold3"))
            link=ld(geno, geno, stats="R.squared", symmetric=T)
            ##link[lower.tri(link)]=NA
            diag(link)=NA
            link[is.na(link)]=0
            link[upper.tri(link)]=NA
            #tm=matrix(1, nrow=nrow(link), ncol=ncol(link))
            #tm[lower.tri(tm)]=0
            #tm=t(tm)
            s <- melt(link)
            s=s[is.na(s[,1])==F,]
            s[, 1]=as.numeric(gsub("[1-5]_", "", s[,1])) #as.numeric(substring(s[,1], 3, ))
            s[, 2]=as.numeric(gsub("[1-5]_", "", s[,2])) #as.numeric(substring(s[,1], 3, ))
            s=s[order(s[,2], decreasing=F),]
            s=s[order(s[,1], decreasing=F),]
            s=s[s[,1]!=s[,2],]
            cols=myPal(100)
            ux=sort(unique(c(s[,1], s[,2])))
            par(mar=c(0,4, 0, 2))
            plot(0, 0, type="n", xlim=r, ylim=c(0, 1), xaxt="n", yaxt="n", xlab="", ylab="")
            segments(seq(min(ux), max(ux), length=length(ux)), 0, ux, 1, lwd=0.1)
            par(mar=c(2,4, 0, 2))
            plot(0, 0, type="n", xlim=r , ylim=c(2+(length(ux)), 0), xaxt="n", yaxt="n", xlab="", ylab="")
            x=seq(min(ux), max(ux), length=length(ux))
            unit=diff(x)[1]
            for(u in 1:length(ux)){
                vals=as.numeric(apply(link, 2, function(x, u){return(as.numeric(na.omit(x))[u])}, u))
                x2=x[u:length(x)]-(u*unit/2)
                y=rep(u, length(x2))
                rect(x2-unit/2, y-0.5,x2+unit/2,y+0.5, col=cols[findInterval(vals, seq(0, 1, length=100))], border=NA, lwd=0)
            }
            legend("bottomleft", "E", bty="n")
            ##add a legend to LD plot
            par(mar=c(2,0, 0, 1))
            plot(0, 0, xlim=c(0,1), ylim=c(0,1), type="n", bty="n", axes=F, xlab="", ylab="")
            colors=myPal(200)
            for(i in 1:200){
                rect(0,(i-1)/200,0.3,i/200, col=colors[i], border=FALSE, lwd=0)
            }
            for(i in seq(0, 1, 0.2)){
                text(0.5, i, paste(i), cex=1)
            }
            dev.off()
        }
    }
}


## GWA_explorer3=function(assoc, bedpath, figpath, threshold=5, win=10000, promoter=1,fst=NULL, h12=NULL, fst_th=0.01){
##     require(readr)
##     require(GenomicRanges)
##     require(TxDb.Athaliana.BioMart.plantsmart22)
##     ##require(VariantAnnotation)
##     require(snpStats)
##     require(dbscan)
##     require(reshape2)
##     ##get genes from arabidopsis
##     txdb=TxDb.Athaliana.BioMart.plantsmart22
##     genes= genes(txdb)
##     ##assoc=read_tsv(pval_file , col_names=T, progress=F)
##     assoc=as.data.frame(assoc)
##     assoc$plog=-log10(assoc$pval)
##     #assoc=assoc[, c("chr", "rs", "ps", "pval","plog")]
##     fam <- paste(bedpath, ".fam", sep="")
##     bim <- paste(bedpath, ".bim", sep="")
##     bed <- paste(bedpath, ".bed", sep="")
##     snps=read.plink(bed, bim, fam)
##     #colnames(assoc)[grep("p_", colnames(assoc))]="pval"
##      ##prepare a Genome wide track to be able to locate the peaks.
##     for (c in 1:5){
##         assign(paste("pv",c,sep=""),assoc[assoc$chr==c,c(1:ncol(assoc))])
##     }
##     pv1$manhattan=pv1[,"ps"]
##     pv2$manhattan=pv2[,"ps"]+max(pv1$manhattan)
##     pv3$manhattan=pv3[,"ps"]+max(pv2$manhattan)
##     pv4$manhattan=pv4[,"ps"]+max(pv3$manhattan)
##     pv5$manhattan=pv5[,"ps"]+max(pv4$manhattan)
##     pv1$col="#133A7F"
##     pv2$col="#518FFF"
##     pv3$col="#2774FF"
##     pv4$col="#4F607F"
##     pv5$col="#1F5CCC"
##     pv=rbind(pv1, pv2, pv3, pv4, pv5)
##     deb=min(pv$manhattan)
##     fin=max(pv$manhattan)
##     assoc=pv
##     top=assoc[assoc[,grep("plog", colnames(assoc))]>=threshold,grep("col", colnames(assoc), inv=T)]
##     ##now get the candidate genes near top SNPs
##     grgwa=with(top, GRanges(chr, IRanges(start=ps-win, end=ps+win)), pval=pval, score=y)
##     ##look for overlap
##     CG=subsetByOverlaps(genes, grgwa)
##     ##for each chromosome, cluster top snps into peaks and plot
##     for(chr in unique(top$chr)){
##         topc=top[top$chr==chr,]
##         if(nrow(topc)==1){grps=1}else{
##             clust=dbscan(as.matrix(topc$ps), eps=2*win, minPts=1)
##             grps=clust$cluster
##         }
##         for(g in unique(grps)){
##             ug=topc[grps==g,]
##             r=range(ug$ps)+c(-win, +win)
##             r2=range(ug$ps)
##             ##LD computation
##             x=snps$genotype[,snps$map$chromosome==chr & snps$map$position>=r[1] & snps$map$position<=r[2]]
##             link=ld(x, x, stats="R.squared")
##             link[lower.tri(link)]=NA
##             ##plotting
##             jpeg(paste(figpath, "_", chr, "_",round(r[1]/1e6, 2), "_", round(r[2]/1e6, 2), "Mb.jpeg",sep=""),res=800, units="in", width=6, height=8)
##             if(is.null(fst)==T & is.null(h12)==T){rd=5}
##             if(is.null(fst)==F & is.null(h12)==F){rd=7}
##             if(is.null(fst)==T & is.null(h12)==F){rd=6}
##             if(is.null(fst)==F & is.null(h12)==T){rd=6}
##             m=matrix(c(1:rd, rd, 1:rd, rd+1), ncol=2)
##             layout(m, widths=c(5, 0.5), heights=c(rep(2, rd-2), 1, 2, 3))
##             par(mar=c(0, 4, 3, 2))
##             plot(assoc$manhattan, assoc$plog, col=assoc$col, pch=16, cex=0.8, xlim=c(deb, fin), xaxt="n", xlab="",ylab="", type="p")
##             points(top$manhattan, top$plog, col="gold3", pch=16, cex=1.2, xlim=c(deb, fin))
##             legend("topleft", "A", bty="n")
##             ##indicate the region zoomed on.
##             rect(assoc$manhattan[assoc$chr==chr & assoc$ps==r2[1]], 0, assoc$manhattan[assoc$chr==chr & assoc$ps==r2[2]], max(abs(top$plog)), border="red", col=NULL, lwd=1)
##             par(mar=c(0, 4, 3, 2))
##             plot_genes(as.data.frame(CG), r=r)
##             legend("topleft", "B", bty="n")
##             axis(3)
##             assocR=assoc[assoc$chr==chr & assoc$ps>=r[1] & assoc$ps<=r[2],]
##             assocR$col[assocR$ps%in%top$ps]="gold3"
##             assocR$lwd=0.2
##             assocR$lwd[assocR$rs%in%top$ps]=3
##             assocR$eff=-log10(assocR[, grep("pval", colnames(assocR))])
##             par(mar=c(0,4, 0, 2))
##             plot(assocR$ps, assocR$eff, xlim=r, ylim=range(assocR$eff), pch=16, xaxt="n", ylab=bquote(-log^10~(italic(p)-value)), xlab="", col=assocR$col)
##             if(is.null(fst)==F){
##                 plot(assocR$ps, assocR$fst, xlim=r, pch=16, xaxt="n", ylab="MAF", xlab="", col=assocR$col, type="p")
##             }
##             if(is.null(h12)==F){
##                 plot(h12$center[h12[,1]==chr], h12$h12[h12[,1]==chr], xlim=r, pch=16, xaxt="n", ylab="H12")
##             }
##             legend("topleft", "C", bty="n")
##                                         #freq=col.summary(x)$MAF
##             myPal <- colorRampPalette(c("deepskyblue4","firebrick3", "gold3"))
##             s <- melt(link)[melt(upper.tri(link))$value,]
##             s[, 1]=as.numeric(substring(s[,1], 3, ))
##             s[, 2]=as.numeric(substring(s[,2], 3, ))
##             s=s[order(s[,2], decreasing=F),]
##             s=s[order(s[,1], decreasing=F),]
##             cols=myPal(100)
##             ux=sort(unique(c(s[,1], s[,2])))
##             par(mar=c(0,4, 0, 2))
##             plot(0, 0, type="n", xlim=r, ylim=c(0, 1), xaxt="n", yaxt="n", xlab="", ylab="")
##             segments(seq(min(ux), max(ux), length=length(ux)), 0, ux, 1, lwd=0.1)
##             par(mar=c(2,4, 0, 2))
##             plot(0, 0, type="n", xlim=r , ylim=c(2+(length(ux)), 0), xaxt="n", yaxt="n", xlab="", ylab="")
##             x=seq(min(ux), max(ux), length=length(ux))
##             unit=diff(x)[1]
##             for(u in 1:length(ux)){
##                 vals=na.omit(as.numeric(apply(link, 1, function(x, u){return(as.numeric(na.omit(x))[u])}, u)))
##             x2=x[u:length(x)]-u*unit/2
##             y=rep(u, length(x2))
##             rect(x2-unit/2, y-0.5,x2+unit/2,y+0.5, col=cols[findInterval(vals, seq(0, 1, length=100))], border=NA, lwd=0)
##         }
##         legend("bottomleft", "D", bty="n")
##         ##add a legend to LD plot
##         par(mar=c(2,0, 0, 1))
##         plot(0, 0, xlim=c(0,1), ylim=c(0,1), type="n", bty="n", axes=F, xlab="", ylab="")
##         colors=myPal(200)
##         for(i in 1:200){
##             rect(0,(i-1)/200,0.3,i/200, col=colors[i], border=FALSE, lwd=0)
##             }
##         for(i in seq(0, 1, 0.2)){
##             text(0.5, i, paste(i), cex=1)
##         }
##         dev.off()
##     }
## }


##an additional funtion for plotting candidate genes
plot_genes=function(ugp, win=diff(range(c(ugp$start,ugp$stop)))/50, r){
    res=matrix(ncol=nrow(ugp), nrow=nrow(ugp))
    y=rep(1:5, length.out=nrow(ugp))
    ## for(i in 1:(nrow(ugp))){
    ##     for(j in i:nrow(ugp)){
    ##         ##do the genes overlap
    ##         g1=as.numeric(ugp[i,c("start", "end")])
    ##         g1=c(min(g1)-win, max(g1)+win)
    ##         g2=as.numeric(ugp[j,c("start", "end")])
    ##         g2=c(min(g2)-win, max(g2)+win)
    ##         if(any(g1>=min(g2) & g1<=max(g2))| any(g2>=min(g1) & g2<=max(g1))){
    ##             res[i,j]=1->res[j,i]}else{res[i, j]=0; res[j, i]=0}
    ##         ##abs(mean(as.numeric(g1))-mean(g2))
    ##     }
    ## }
    ## res=as.dist(res)
    ## track=cutree(hclust(res), h=0)
    track=y
    plot(0, 0, xlim=r, ylim=c(0.5, max(track)+1), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
    sens=ugp[ugp$strand=="+",]
    antisens=ugp[ugp$strand=="-",]
    arrows(x0=sens[,"start"], y0=track[ugp$strand=="+"], x1=sens[,"end"], y1=track[ugp$strand=="+"], angle=45, length=0.08, code=2)
    arrows(x0=antisens[,"end"], y0=track[ugp$strand=="-"], x1=antisens[,"start"], y1=track[ugp$strand=="-"], angle=45, length=0.08, code=2)
    text(apply(ugp[, c("start", "end")], 1, mean), track, label=ugp$gene, cex=0.8, adj=c(0.5, -1))
}






gwa_explorer=function(pval_file, g_folder, go_file, maf, threshold, win, promoter=1, get_GO=T,GO_enrich=F, plot_nodes=F, go_map_file="~/Documents/data/postdoc/geno_GWA/TAIR10/map_go.txt", N_genes=20){
  if(GO_enrich==T){require(topGO)}
  out_file=paste(pval_file, "_cg.txt", sep="")
  go_out=paste(pval_file, "_go.txt", sep="")
##pull out top snps and the candidate genes
system(paste("get_top_snps.py -m ",maf," -p ",pval_file," -t ", threshold, " > ",pval_file,"_snps.txt",sep=""))
system(paste("get_candidates.py -t ",pval_file,"_snps.txt  -g ", g_folder, " -w  ", win," -p ", promoter, "  > ", out_file, sep=""))
  ##get go annotations
  if(get_GO==T){system(paste("get_go.py -v 1 -g ", go_file, " -t ", out_file, " > ",go_out, sep=""))}
  ##read the outputs
  top=read.table(out_file, sep="\t", h=F)
  ##do the manhattan plot again, with some gene annotations.
  pval=read.table(pval_file, skip=1, sep=",")
  x=file(pval_file)
  z=readLines(x, n=1)
  N=as.integer(unlist(strsplit(unlist(strsplit(z, ";"))[1], "="))[2])
  min_N=ceiling(N*maf)
  pval=pval[pval[,7]>=min_N & pval[,7]<=(N-min_N) ,]
  g=top[top[,9], c(9,1, 10, 11)]
  g=droplevels(unique(g))
  g$score=apply(g, 1, function(x, y){max(-log10(y[y[,9]==x[1],3]))[1]}, y=top)
  g=g[order(g$score, decreasing=T),]
  close(x)
  if(nrow(g)>N_genes){g=g[1:N_genes,]}
  genes=g[,1]
  pdf(paste(pval_file, "_manhattan.pdf", sep=""), paper="special", width=10, height=3)
  a=manhattan(pval=pval,top=top, genes=genes, close_device=T)
  write.table(a, paste(pval_file, "_annotations.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
  if(GO_enrich==T){
    go=read.table(go_out, sep="\t", h=F, na.string="NA")
    geneID2GO <- readMappings(go_map_file)
    geneNames <- names(geneID2GO)
    myInterestingGenes=unique(top[,9])
    geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
    names(geneList) <- geneNames
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    allRes <- GenTable(GOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 20)
    if(plot_nodes==T){
    pdf(paste(pval_file, "_go_nodes.pdf", sep=""), paper="special", width=10, height=10)
    showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 20, useInfo = 'all')
    dev.off()
    }
    write.table(allRes, paste(pval_file,"_20enriched_goterms.txt", sep=""), row.names=F, col.names=T, quote=F)
  }
  allgenes=unique(top[,9])
  write.table(allgenes, paste(pval_file,"_all_genes.txt", sep=""), row.names=F, col.names=F, quote=F)
  genes_scatter(top=top, pval_file=pval_file)
  return("Done!")
}


##mahattan plots with candidates

manhattanCG=function(pval, yaxis="pval", cg=NULL){
    require(wordcloud)
    #pval$chr=as.numeric(substring(pval$snp, 1, 1))
    #pval$pos=as.numeric(gsub("[1-5][:-:_: :]+", "", pvals[,1]))
    for (c in 1:5){
        assign(paste("pv",c,sep=""),pval[pval$chr==c,c(1:ncol(pval))])
    }
    pv1$manhattan=pv1[,"pos"]
    pv2$manhattan=pv2[,"pos"]+max(pv1$manhattan)
    pv3$manhattan=pv3[,"pos"]+max(pv2$manhattan)
    pv4$manhattan=pv4[,"pos"]+max(pv3$manhattan)
    pv5$manhattan=pv5[,"pos"]+max(pv4$manhattan)
    pv1$col="#133A7F"
    pv2$col="#518FFF"
    pv3$col="#2774FF"
    pv4$col="#4F607F"
    pv5$col="#1F5CCC"
    pv=pval
    pv$manhattan=c(pv1$manhattan,pv2$manhattan,pv3$manhattan,pv4$manhattan,pv5$manhattan)
    pv$col=c(pv1$col,pv2$col,pv3$col,pv4$col,pv5$col)
    deb=min(pv$manhattan)
    fin=max(pv$manhattan)
    m=matrix(c(2, 3, 0, 1), ncol=2, byrow=T)
    layout(m,c(2,30),c(10,2))
    par(mar=c(0,0.7,0,0.5))
    plot(0,0,xlim=c(deb,fin),ylim=c(0,2),type="n",axes=FALSE,frame.plot=FALSE, xaxs="i", yaxs="i")
    for (c in 1:5){
        X=min(pv$manhattan[pv[,"chr"]==c])+((max(pv$manhattan[pv[,"chr"]==c])-min(pv$manhattan[pv[,"chr"]==c])))/2
        text(X,1,paste("Chromosome",c),cex=1, offset=0)
    }
    par(mar=c(0,0,0,0), ps=8)
    plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
    if(yaxis=="pval"){
        text(0.06,0.5,bquote(-log^10~(italic(p)-value)),srt=90,cex=1.5)
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        ymax=max(-log10(pv[,"pval"]), na.rm=T)+1
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        plot(pv$manhattan,-log10(pv[,"pval"]),xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(1,ymax),axes=FALSE, cex =0.5 ,  col =pv$col ,pch=16, xaxs="i", yaxs="i")}
    if(yaxis=="adjr2"){
        text(0.06,0.5,bquote(adjusted ~ italic(r^2)),srt=90,cex=1.5)
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        ymax=max(pv[,"adjr2"], na.rm=T)+0.01
        par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
        plot(pv$manhattan,pv[,"adjr2"],xlab="",ylab="",frame.plot=FALSE,xlim=c(deb,fin),ylim=c(0.01,ymax),axes=FALSE, cex =0.5 ,  col =pv$col ,pch=16, xaxs="i", yaxs="i")}
    axis(1,cex.axis=1,lwd=0.5,xaxp=c(deb,fin,10),labels=FALSE, tck=0)
    axis(2,cex.axis=1,lwd=0.5,tck=0.03)
    box(lwd=0.5)
    cg$mh1=0
    cg$mh1=0
    if(is.null(cg)==FALSE){
        cg$mh1[cg[,2]==1]=cg[cg[,2]==1,3]
        cg$mh2[cg[,2]==1]=cg[cg[,2]==1,4]
        cg$mh1[cg[,2]==2]=cg[cg[,2]==2,3]+max(pv1$manhattan)
        cg$mh2[cg[,2]==2]=cg[cg[,2]==2,4]+max(pv1$manhattan)
        cg$mh1[cg[,2]==3]=cg[cg[,2]==3,3]+max(pv2$manhattan)
        cg$mh2[cg[,2]==3]=cg[cg[,2]==3,4]+max(pv2$manhattan)
        cg$mh1[cg[,2]==4]=cg[cg[,2]==4,3]+max(pv3$manhattan)
        cg$mh2[cg[,2]==4]=cg[cg[,2]==4,4]+max(pv3$manhattan)
        cg$mh1[cg[,2]==5]=cg[cg[,2]==5,3]+max(pv4$manhattan)
        cg$mh2[cg[,2]==5]=cg[cg[,2]==5,4]+max(pv4$manhattan)
        for(r in 1:nrow(cg)){
            rect(cg[r, "mh1"], 0, cg[r,"mh2"], ymax-0.1*ymax, col="firebrick2", lwd=0.5, border="firebrick2")
        }
        cg2=cg
        ##[cg[,"cg"]!=".",]
        if(nrow(cg)>1){
        textplot(x=cg2$mh1, y=rep((ymax-0.05*ymax), nrow(cg2)),paste(cg2[,"cg"]), new=F, cex=1, show.lines=TRUE)}else{text(x=cg2$mh1, y=rep((ymax-0.05*ymax), nrow(cg2)),paste(cg2[,"cg"]))}
    }
}



zoom_manhattan=function(pval, yaxis="pval", chr=1, start=0, stop=100000, main=""){
    pv=pval[pval$chr==chr & pval$pos>=start & pval$pos<=stop,]
    col=c("#133A7F", "#518FFF", "#2774FF","#4F607F", "#1F5CCC")[chr]
    deb=min(pv$pos)
    fin=max(pv$pos)
    xlab=paste("Chromosome",chr, round(start/10000, 2),  "-", round(start/10000, 2))
    ylab=bquote(-log^10~(italic(p)-value))
    par(mar=c(4,4,1,0.5),ps=8,mgp=c(1.5,0.6,0))
    ymax=max(-log10(pv[,"pval"]), na.rm=T)+1
    plot(pv$pos,-log10(pv[,"pval"]),xlab=xlab,ylab=ylab,frame.plot=FALSE,xlim=c(deb,fin),ylim=c(0,ymax),axes=FALSE, cex =0.5 ,  col =col ,pch=16, xaxs="i", yaxs="i", main=main)
    axis(1,cex.axis=1,lwd=0.5)
    axis(2,cex.axis=1,lwd=0.5,tck=-0.03)
    box(lwd=0.5)
}


shuffle=function(v, chr){
    L=length(chr)
    ##v is the vector to shuffle, chr is the chromosome vector, and pos are the position on a chromosome
    rv=unlist(lapply(as.list(1:5), function(i, v, ch, rch, ror){
        x=v[ch==rch[i]]
        if(ror[i]==-1){x=rev(x)}
        return(x)
    }, v=1:L, ch=chr, rch=sample(1:5), ror=sample(c(1, -1), size=5, replace=T)))
    ##spin rv around by a random number
    sp=sample(1:L, size=1)
    rv=rv+sp
    rv[rv>L]=1:sp
    return(v[rv])
}




###running CAVIAR

CAVIAR_scan=function(pval_file, bedpath, respath, threshold=5, win_clust=20000, win_caviar=500){
    require(readr)
    require(snpStats)
    require(dbscan)
    require(reshape2)
    system(paste("mkdir ", respath, sep=""))
    assoc=read_tsv(pval_file , col_names=T, progress=F)
    assoc=as.data.frame(assoc)
    colnames(assoc)[grep("p_", colnames(assoc))]="pval"
    assoc$plog=-log10(assoc$pval)
    fam <- paste(bedpath, ".fam", sep="")
    bim <- paste(bedpath, ".bim", sep="")
    bed <- paste(bedpath, ".bed", sep="")
    snps=read.plink(bed, bim, fam)
    top=assoc[assoc[,grep("plog", colnames(assoc))]>=threshold,grep("col", colnames(assoc), inv=T)]
    for(chr in unique(top$chr)){
        topc=top[top$chr==chr,]
        if(nrow(topc)==1){grps=1}else{
            clust=dbscan(as.matrix(topc$ps), eps=2*win_clust, minPts=1)
            grps=clust$cluster
        }
        for(g in unique(grps)){
            ug=topc[grps==g,]
            r=range(ug$ps)+c(-win_clust, +win_clust)
            ##find the most associated SNP
            pt=ug[ug$plog==max(ug$plog),"ps"]
            ##MAKE A SMALLE WINDOW AROUND IT.
            rr=range(pt)+c(-win_caviar, +win_caviar)
            sr=assoc[assoc$chr==chr & assoc$ps>=rr[1] & assoc$ps<=rr[2],]
            sr$zscore=(sr$plog-mean(assoc$plog))/sd(assoc$plog)
            ##LD computation
            geno=snps$genotype[,snps$map$chromosome==chr & snps$map$position>=r[1] & snps$map$position<=r[2] & snps$map$position%in%sr$ps]
            link=ld(geno, geno, stats="Covar", symmetric=F)
            ##save input files for caviar
            zf=paste(respath,"zscore_", chr, "_",round(r[1]/1e6, 2), "_", round(r[2]/1e6, 2),".txt", sep="")
            write.table(sr[,c("rs", "zscore")], zf, sep="\t", col.names=F, row.names=F, quote=F)
            lf=paste(respath,"ld_", chr, "_",round(r[1]/1e6, 2), "_", round(r[2]/1e6, 2),".txt", sep="")
            write.table(link, lf, sep="\t", col.names=F, row.names=F, quote=F)
            ##define output
            outf=paste(respath,"caviar_out_", chr, "_",round(r[1]/1e6, 2), "_", round(r[2]/1e6, 2), sep="")
            outf2=paste(respath,"caviar_out_", chr, "_",round(r[1]/1e6, 2), "_", round(r[2]/1e6, 2), "_modelsearch", sep="")
            ##run CAVIARBF
            system(paste("caviarbf -o ", outf, " -r ", lf, " -z ", zf, " -a 0.1,0.2,0.4 -n 200 -c 5 -t 0 " , sep=""))
            system(paste("model_search -i ", outf, " -o ", outf2, " -m ",nrow(link), " -p 0 ",  sep=""))
        }
    }
}





##ph2 from gemma

ph2gemma=function(f){
    con=file(f)
    line=readLines(f) 
    ve=as.numeric(unlist(strsplit(line[grep("## ve estimate in the null model", line)], "="))[2])
    vg=as.numeric(unlist(strsplit(line[grep("## vg estimate in the null model", line)], "="))[2])
    return(vg/(ve+vg))
}    


##enrichment function

enrich=function(gwa, th=4, reg, NullSize=1000){
    top=gwa[gwa$score>=th,]
    grgwa=with(top, GRanges(chr, IRanges(start=pos, end=pos)), score=score)
    grsel=with(reg, GRanges(chr, IRanges(start=start, end=stop)), alpha=alpha, N_genes=N_genes)
    count=countOverlaps(grsel,grgwa, type="any")
    L=length(count[count>0])
    ##determine the size of each region
    reg$size=reg$stop-reg$start
    n=c()
    for(i in 1:NullSize){
        N=nrow(reg)
        R=gwa[sample(1:nrow(gwa), replace=F, size=N),c("chr", "pos")]
        R$start=R$pos-reg$size/2
        R$stop=R$pos+reg$size/2
        grR=with(R, GRanges(chr, IRanges(start=start, end=stop)))
        count=countOverlaps(grR,grgwa, type="any")
        n=append(n, length(count[count>0]))
    }
    e=ecdf(n)
    return(c(Overlaps=L, emp_pval=1-e(L)))
}


##function to perform GWA with Gemma, chr by chr, and for each chr, using a K matrix calculting excluding the chr investigated.

##test variables
## plink_prefix="./GWA/snps/sweden_final"
## ph2=readRDS("./res/ph2_hh.rds")
## ph2CI=read.table("./res/ph2_CI_albi.txt", h=T)
## ph2=data.frame(hh=names(ph2),ph2CI)
## traitsnames=c(paste(ph2[,1]), "fecundity")

## traits=c("F8", "F5", "B38", "fecundity")

## ldth=0.8
## snps=unlist(readRDS(paste("./GWA/prunned_snps_LDth_", ldth, ".rds", sep="")))
## pruned_SNPs=data.frame(chr=as.integer(substring(snps, 1, 1)),snps) 
## temp_folder="./smartGWAtest_temp"
## res_folder="./smartGWAtest_out"


## smartGWA=function(plink_prefix="./GWA/snps/sweden_smartGWA", traitsnames, pruned_SNPs, res_folder="./output_smartGWA", temp_folder="./tempGWA/"){
##     require(readr)
##     system(paste("mkdir ", res_folder, sep=""))
##     system(paste("mkdir ", temp_folder, sep=""))
##     chrs=1:5
##     for(chr in 2:5){
##         write.table(pruned_SNPs[pruned_SNPs$chr==chr,2], paste(temp_folder, "/SNPs_chr",chr, ".keep", sep=""), col.names=F,row.names=F, quote=F)
##         ##make a bed file with LD prunned SNPs for all other chr
##         system(paste("plink --bfile ", plink_prefix, " --make-bed --allow-no-sex --missing-code N,-9,0,NA --extract ", temp_folder, "/SNPs_chr", chr, ".keep -out ", temp_folder, "/sweden_chr", chr,  "_for_K", sep=""))
##         ##make a bed file with all SNPs on chr
##         system(paste("plink --bfile ", plink_prefix, " --make-bed --chr ",chr," --allow-no-sex --missing-code N,-9,0,NA -out ", temp_folder, "/sweden_chr", chr,  "_for_assoc", sep=""))
##         ##replace the fam file with the original one
##         system(paste("cp ", plink_prefix, ".fam ", temp_folder, "/sweden_chr", chr,  "_for_assoc.fam",sep=""))
##         ##run K computation with gemma
##         system(paste("gemma -gk 1 -bfile ",temp_folder, "/sweden_chr", chr, "_for_K -n 1 -o K_chr", chr, "excluded", sep=""))
##         for(trait in traits){
##             ##run association analysis with gemma
##             system(paste("gemma -lmm 2 -bfile ",temp_folder, "/sweden_chr", chr, "_for_assoc -n ", which(traitsnames==trait)," -k ./output/K_chr", chr, "excluded.cXX.txt -o lmm_", trait,"_chr", chr, sep=""))
##         }
##     }
##     for(trait in traits){
##         ##combine the association outputs for the 5 chr into one.
##         system(paste("cat ./output/lmm_", trait,"_chr", 1,".assoc.txt> ", res_folder, "/lmm_", trait, "allchr.assoc.txt", sep=""))
##         ##for the other chr, skip first line and apped to file
##         for(chr in 2:5){
##             system(paste("cat ./output/lmm_", trait,"_chr", chr,".assoc.txt | tail -n+2 >> ", res_folder, "/lmm_", trait, "allchr.assoc.txt", sep=""))
##         }
##         ##make the manhattan plot directly
##         file=paste(res_folder, "/lmm_", trait,"allchr.assoc.txt", sep="")
##         dat=read_tsv(file, col_names=T, progress=F)
##         p=as.data.frame(dat)
##         colnames(p)[3]="pos"
##         colnames(p)[9]="pval"
##         jpeg(paste(res_folder, "/lmm_",trait, "allchr.jpeg", sep=""), res=600, unit="in", height=4, width=8, pointsize=8)
##         manhattan(p)
##         dev.off()
##     }
## }


        
