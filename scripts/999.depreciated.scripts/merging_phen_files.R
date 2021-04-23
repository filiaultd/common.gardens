

## ---- phenmerge

##reading the randomization file.

acc=read.table("./data/acc_list.txt", sep="\t", h=T)


##reading orginal data files for each phenotype

##rosette data (for bothe years):
size=read.table("./data/rosette_size.txt", sep="\t", h=T)
size$experiment=toupper(size$experiment)
size$combi=paste(size$experiment, size$tray, size$row, size$col, sep="_")
color=read.table("./data/color_rosettes.txt", sep="\t", h=T)
color$exp=toupper(color$exp)
color$combi=paste(color$exp, color$tray, color$row, color$col, sep="_")

##I need to keep only one time point per year per experiment in this file.
size=size[(size$date=="20111119" & size$experiment=="ULL")==F,]
size=size[(size$date=="20121120" & size$experiment=="RAT")==F,]
color=color[(color$date=="20111119" & color$exp=="ULL")==F,]
color=color[(color$date=="20121120" & color$exp=="RAT")==F,]

##compare the area column from the color and size data
size$area_color=color$area[match(size$combi, color$combi)]
plot(size$area, size$area_color, col=as.factor(size$experiment))
##Here the area in pixel number (in color) and the area in cm^2 (size) doesn't have a correlation of one. it might be related to the difference in height at which the pictures were taken (Svante vs. Ben), the number of pixels on the cameras (Mia's Camera vs Nikon D60 later)...

##the stockiness is not in the size data.
size$stockiness=(4*pi*size$area/(size$perimeter^2))

#################################################
#######   2011 experiments   ####################
#################################################

##read the survival data

surv=read.table("./data/survival_2011.txt", h=T, sep="\t")
surv$id=acc$lines[match(surv$line, acc$tubes)]
surv$name=acc$name[match(surv$line, acc$tubes)]

##clean up surv (typos).

surv$spring[surv$spring==11]=1
surv$spring[surv$spring==10]=NA
surv$spring[surv$spring=="?"]=NA
surv$spring[surv$spring=="-"]=NA
surv=droplevels(surv)

##subset size and color to keep only the 2011 data

size2011=size[size$year==2011,]
color2011=color[color$year==2011,]

##read the herbivory data scored in RAT

herb=read.table("./data/rat.snail.2011.csv", sep=",", h=T)[,1:5]

##read the last fecundity dataset for 2011

fecundity=read.table("./data/fecundity_2011.txt", sep="\t", h=T)
##swopping rows and cols for fecundity to be consistent with the other datasets
#x=colnames(fecundity)[c(1:3, 5, 4, 6:ncol(fecundity))]
#colnames(fecundity)=x

##now merge all this, based on surv.

surv$combi=paste(surv$exp, surv$tray, surv$row, surv$column, sep="_")
fecundity$combi=paste(fecundity$exp, fecundity$tray, fecundity$row, fecundity$col, sep="_")
herb$combi=paste("RAT", herb$tray, herb$row, herb$column, sep="_")

##now built the data table.

data2011=surv

##having the planting date in the data would be nice.

dates=data.frame(expand.grid(c("ADA", "RAM", "ULL", "RAT"), c("A", "B", "C")))
dates=dates[order(dates[,1]),]
colnames(dates)=c("exp", "block")
dates$planting=as.Date(c("2011-08-08","2011-08-10", "2011-08-12", "2011-08-07", "2011-08-09", "2011-08-11", "2011-08-31","2011-09-02", "2011-09-4", "2011-09-01", "2011-09-03", "2011-09-05"))
dates$comb=paste(dates$exp, dates$block, sep="_")
comb=paste(data2011$exp,data2011$block, sep="_")
data2011$planting_date=dates$planting[match(comb, dates$comb)]

##match in the phenotypes from size2011
data2011$rosette_date=size2011[match(data2011$combi, size2011$combi),"date"]
data2011$area=size2011[match(data2011$combi, size2011$combi),"area"]
data2011$perimeter=size2011[match(data2011$combi, size2011$combi),"perimeter"]
data2011$max_diameter=size2011[match(data2011$combi, size2011$combi),"max_diameter"]
data2011$sdR=size2011[match(data2011$combi, size2011$combi),"sdR"]
data2011$circle_area=size2011[match(data2011$combi, size2011$combi),"circle_area"]
data2011$stockiness=size2011[match(data2011$combi, size2011$combi),"stockiness"]
##match in the color phenotype
data2011$color=color2011[match(data2011$combi, color2011$combi),"color"]
##match in the herbivory scores
data2011$herbivory=herb[match(data2011$combi, herb$combi),"slug"]
##match in the fecundity estimates
data2011$fecundity=fecundity$fecundity[match(data2011$combi, fecundity$combi)]
##make the rosette_date a working date column.
x=data2011$rosette_date
y=as.Date(paste(substring(x, 1, 4),substring(x, 5, 6), substring(x, 7,8), sep="-"))
data2011$rosette_date=y
##reorder the columns a little
data2011=data2011[,c("exp", "block", "tray", "row", "column", "line", "id", "name", "planting_date", "errors", "fall", "spring", "sampled", "rosette_date", "area", "perimeter", "max_diameter", "sdR", "circle_area", "stockiness","color", "herbivory", "fecundity")]
##save it!
saveRDS(data2011, file="./data/data2011.rds")
write.table(data2011, "./data/data2011.txt", col.names=T, row.names=F, quote=F, sep="\t")


#################################################
#######   2012 experiments   ####################
#################################################

##read the survival data
surv=read.table("./data/survival_2012.txt", h=T, sep="\t")
surv$id=acc$lines[match(surv$line, acc$tubes)]
surv$name=acc$name[match(surv$line, acc$tubes)]

##put it in better shape
surv$combi=paste(surv$exp, surv$tray, surv$row, surv$column, sep="_")
surv$block[surv$tray<=27]="A"
surv$block[surv$tray>=28 & surv$tray<=54]="B"
surv$block[surv$tray>=55]="C"

##subset size to keep only data for the 2012 experiments
size2012=size[size$year==2012,]
color2012=color[color$year==2012,]

##read some flowering time data, from the North (FTN) and the South (FTS)
FTS=read.table("./data/FT_2012_South_dates.txt", sep="\t", h=T)
FTN=read.table("./data/FT_2012_North_dates.txt", sep="\t", h=T)

##read the 2012 fecundity estimates
fecundity=read.table("./data/fecundity_2012.txt", sep="\t", h=T)
##swopping rows and cols for fecundity.
x=colnames(fecundity)
#[c(1:3, 5, 4, 6:ncol(fecundity))]
#colnames(fecundity)=x
fecundity$combi=paste(fecundity$exp, fecundity$tray, fecundity$row, fecundity$col, sep="_")

##combine all data based on surv.

data2012=surv

##clean up the surv column

data2012$survival_03162013[data2012$survival_03162013==11]=1
data2012$survival_03162013[data2012$survival_03162013==","]=NA
data2012=droplevels(data2012)

##add the fall flowering and survival column from the phenotyping Rod did on the images.

for(e in c("ADA", "RAM", "ULL", "RAT")){
if(e=="ADA"){fbw=cbind(e, read.table(paste("./data/fbw_2012_", e, ".txt", sep=""), sep="\t", h=T))}else{fbw=rbind(fbw, cbind(e, read.table(paste("../all_phenotypes/fbw_2012/fbw_2012_", e, ".txt", sep=""), sep="\t", h=T)))}
}
fbw$combi=paste(fbw$e, fbw$tray, fbw$row, fbw$column, sep="_")


##clean it up

fbw[fbw$ft_fall=="9","ft_fall"]=0
fbw[fbw$ft_fall=="no photo","ft_fall"]=NA

##use this flowering before winter as a fall column. It's the same as what was done in 2011 but photos were phenotype the same evening for that.
data2012$fall=fbw[match(data2012$combi, fbw$combi),"ft_fall"]

##add planting date

dates=data.frame(expand.grid(c("ADA", "RAM", "ULL", "RAT"), c("A", "B", "C")))
dates=dates[order(dates[,1]),]
colnames(dates)=c("exp", "block")
dates$planting=as.Date(c("2012-08-08","2012-08-10", "2012-08-12", "2012-08-07", "2012-08-09", "2012-08-11", "2012-08-31","2012-09-02", "2012-09-4", "2012-09-01", "2012-09-03", "2012-09-05"))
dates$comb=paste(dates$exp, dates$block, sep="_")
comb=paste(data2012$exp,data2012$block, sep="_")
data2012$planting_date=dates$planting[match(comb, dates$comb)]

##add the rosette data
data2012$rosette_date=size2012[match(data2012$combi, size2012$combi),"date"]
data2012$area=size2012[match(data2012$combi, size2012$combi),"area"]
data2012$perimeter=size2012[match(data2012$combi, size2012$combi),"perimeter"]
data2012$max_diameter=size2012[match(data2012$combi, size2012$combi),"max_diameter"]
data2012$sdR=size2012[match(data2012$combi, size2012$combi),"sdR"]
data2012$circle_area=size2012[match(data2012$combi, size2012$combi),"circle_area"]
data2012$stockiness=size2012[match(data2012$combi, size2012$combi),"stockiness"]
##add the color data
data2012$color=color2012[match(data2012$combi, color2012$combi),"color"]
##add the fecundity data
data2012$fecundity=fecundity$fit[match(data2012$combi, fecundity$combi)]

##add flowering time in the spring data (not available for 2011)
colnames(FTS)[5]="col"
FT=rbind(FTN, FTS)

FT$combi=paste(FT$exp, FT$tray,FT$row, FT$col, sep="_")
data2012$flowering_date=FT[match(data2012$combi, FT$combi),"flowering_date"]
data2012$FT=FT[match(data2012$combi, FT$combi),"FT"]

##make the rosette_date a working date column.

x=data2012$rosette_date
y=as.Date(paste(substring(x, 1, 4),substring(x, 5, 6), substring(x, 7,8), sep="-"))
data2012$rosette_date=y

##reorder col

data2012=data2012[,c("exp", "block", "tray", "row", "column", "line", "id", "name", "planting_date", "errors", "fall", "survival_03162013", "epi", "PC_sampled", "rosette_date", "area", "perimeter", "max_diameter", "sdR", "circle_area", "stockiness", "color", "flowering_date", "FT", "fecundity")]

colnames(data2012)[match("survival_03162013", colnames(data2012))]="spring"

##save it:
saveRDS(data2012, file="./data/data2012.rds")
write.table(data2012, "./data/data2012.txt", col.names=T, row.names=F, quote=F, sep="\t")

## ---- end-of-phenmerge


