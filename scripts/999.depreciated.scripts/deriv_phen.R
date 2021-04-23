
##A0001##
##read the 2011 and the 2012 data

## ---- ows
d11=readRDS("./data/data2011.rds")
d12=readRDS("./data/data2012.rds")


##add the ows column.
d11$ows=NA
d11$ows[d11$fall%in%c(1, 2,4, 5) & d11$spring%in%c(1, 2, 4, 5)]=1
d11$ows[d11$fall%in%c(1, 2,4, 5) & d11$spring%in%c(0, 3)]=0 ##A0002##

d12$ows=NA
d12$ows[d12$fall%in%c(1, 2, 4, 5) & d12$spring%in%c(1, 2, 4, 5)]=1
d12$ows[d12$fall%in%c(1, 2, 4, 5) & d12$spring%in%c(0, 3) & d12$epi==F]=0 ##A0003##

## ---- end-of-ows

## ---- sss

##add the sss column
##A0004##
d11$sss=NA
d11$sss[d11$fall%in%c(1, 2, 4, 5) & is.na(d11$fecundity)==F & d11$sampled==F]=1
d11$sss[d11$fall%in%c(1, 2, 4, 5) & is.na(d11$fecundity)==T  & d11$sampled==F]=0
d12$sss=NA
d12$sss[d12$fall%in%c(1, 2,  4, 5) & is.na(d12$fecundity)==F & d12$PC_sampled==F & d12$epi==F]=1
d12$sss[d12$fall%in%c(1, 2,  4, 5) & is.na(d12$fecundity)==T  & d12$PC_sampled==F & d12$epi==F]=0

## ---- end-of-sss

## ---- verifowssss
###Do some checking

t1=table(d11[,c("ows", "sss")], exclude=NULL)
t2=table(d12[, c("ows", "sss")], exclude=NULL)

kable(data.frame(t1), caption="Number of data points for each value of ows and sss in the 2011 experiments")
kable(data.frame(t2), caption="Number of data points for each value of ows and sss in the 2012 experiments")

## ---- end-of-verifowssss


##table(d11[is.na(d11$fecundity)==F,c("ows", "sss")])
##table(d12[is.na(d12$fecundity)==F,c("ows", "sss")])

##checking the empty wells. ##A0006##

 #d11[d11$line=="empty" & d11$fall!=0,]
#d11[d11$line=="empty" & d11$fall!=0,]

## ---- cleanowssss
##removing theoritically "empty" pots from the datasets ##A0007##

d11=droplevels(d11[d11$line!="empty",])
d12=droplevels(d12[d12$line!="empty",])

t1=data.frame(table(ows=d11$ows, sss=d11$sss))
t2=data.frame(table(ows=d12$ows, sss=d12$sss))

kable(t1,caption="overwinter survival and survival to seed set in 2012") 
kable(t2,caption="overwinter survival and survival to seed set in 2013")

d11=d11[(d11$ows==0 & d11$sss==1)==F,]
d12=d12[(d12$ows==0 & d12$sss==1)==F,]

t1=data.frame(table(ows=d11$ows, sss=d11$sss))
t2=data.frame(table(ows=d12$ows, sss=d12$sss))

kable(t1,caption="overwinter survival and survival to seed set in 2012, after clean-up") 
kable(t2,caption="overwinter survival and survival to seed set in 2013, after clean-up")

## ---- end-of-cleanowssss

##cleaning up the rosette size data.

## ---- rosettes1

##A0008##
kable(table(d11[is.na(d11$area)==F,"exp"]), caption="number of rosettes measured per experiment in 2011. Keeping in mind that this dataset only includes one date but others may be available (see specific org file for analysis of the rosette data). In 2011 we didn't get rosette data because plants grew too big and were overlapping too much for photo analysis.")
kable(table(d12[is.na(d12$area)==F,"exp"]), caption="number of rosettes measured per experiment in 2012. Keeping in mind that this dataset only includes one date but others may be available (see specific org file for analysis of the rosette data). In 2011 we didn't get rosette data because plants grew too big and were overlapping too much for photo analysis.")

## ---- end-of-rosettes1

## ---- rosettes2
##A0009##
d11=d11[(is.na(d11$area)==F & d11$fall%in%c(1, 2, 4, 5)==F)==F,]
d12=d12[(is.na(d12$area)==F & d12$fall%in%c(1, 2, 4, 5)==F)==F,]
## ---- end-of-rosettes2

##rosette data.frame
##A0010##
## ---- rosettesage
d11$age_rosettes=as.Date(d11$rosette_date)-as.Date(d11$planting_date)
d12$age_rosettes=as.Date(d12$rosette_date)-as.Date(d12$planting_date)
## ---- end-of-rosettesage

## ---- addfittrait
##setting a "fitness" column
d11$fitness=d11$fecundity
d11$fitness[d11$sss==0 & is.na(d11$sss)==F]=0

d12$fitness=d12$fecundity
d12$fitness[d12$sss==0 & is.na(d12$sss)==F]=0

## ---- end-of-addfittrait

##removing the controls
## ---- controlsandco
##A0012##
c11=d11[d11$line%in%c("C1", "C2"),]
d11=d11[d11$line!="empty",]
d11=d11[d11$line%in%c("C1", "C2")==F,]
d11=d11[d11$error==".",]

c12=d12[d12$line%in%c("C1", "C2"),]
d12=d12[d12$line!="empty",]
d12=d12[d12$line%in%c("C1", "C2")==F,]
##quick fix of the error column in the 2012 data.
d12$error[d12$error==""]="."
d12$error[d12$error!="."]="error"
##and finish removing lines with errors
d12=droplevels(d12[d12$error==".",])

##removing Edi-0, Col-FRI-FLC, Col-FRI
##A0017##
d11=d11[d11$line%in%c("201", "202", "203")==F,]
d12=d12[d12$line%in%c("201", "202", "203")==F,]

## ---- end-of-controlsandco
#}}}

## ---- savedatasets
##save a final cleaned up version of the data

saveRDS(d11, file="./data/d11_for_analysis.rds")
saveRDS(d12, file="./data/d12_for_analysis.rds")


## ---- end-of-savedatasets
