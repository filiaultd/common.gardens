##here I assign accessions to a region based on the dominant admixture cluster they have

acc=read.table("./data/acc_list.txt", h=T, sep="\t")
acc=droplevels(acc[acc$tubes<=200,])
colnames(acc)[1]="id"
adm=read.table("./data/1001genomes-accessions and 1001genomes-admixture.csv", h=T, sep=",")

##choose K
K="K9"
##filter to keep only swedish lines. The K9 clustering seems to capture what we want and is the default on the website.
adm=droplevels(adm[adm$id%in%acc$id,c(1, 11,grep(K, colnames(adm)))])

##cluster 3 is the Northern sweden group
## cluster 6 is the Southern sweden group
adm=adm[,c(1,2, which(colnames(adm)%in%c(paste(K, ".ChrAll.3", sep=""), paste(K, ".ChrAll.6", sep=""))))]

##offcourse, only 186 of our lines are in the 1001 genomes dataset, so I'll be missing info.

m=merge(acc, adm, by="id")
table(m$region, m$group)

##I'll make a new group col which will re-assigne accessions that are not already assigned to Northern or Southern sweden groups based on which of the cluster 3 or 6 they have most.

m$newgroup=apply(m, 1, function(x){if(x[6]>x[7]){return("north_sweden")}else{return("south_sweden")}})

cs=acc[acc$region=="C Sweden",]

m=merge(cs, adm, by="id")


##This is very imperfect because some of Swedish line from the South are assigned to germany.



##determine ad-hoc which is the dominant Nothern and Southern
