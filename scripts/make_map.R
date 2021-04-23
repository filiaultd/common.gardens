

## ---- map
### mods DLF 04 July 19
### needed new API key for google

### mods DLF 24 Oct 19
### update to new geo locations

library(ggmap)
library(grid)

#clim=read.table("./data/worldclim_swedish_acc.txt", sep="\t", h=T, stringsAsFactor=T )
#clim_exp=read.table("./data/worldclim_swedish_exp.txt" , sep="\t", h=T)
#colnames(clim)[1]="id"

clim=read.table("./data/bioclim.v2.200.experimental.lines.txt", stringsAsFactors=FALSE)
clim_exp=read.table("./data/bioclim.v2.200.experimental.sites.txt", stringsAsFactors=FALSE)
colnames(clim)[1]="id"


#CG=clim_exp[clim_exp$experiments=="Common garden",c(2:4)]
#exp.sites <- clim_exp[clim_exp$experiments%in%c("Common garden","selection exp"),]
exp.sites <- clim_exp
exp.sites <- exp.sites[c(1,3,4,8:15),]
exp.sites <- exp.sites[order(exp.sites$lat, decreasing=TRUE),]
### so some of these are nearly duplicates.  I am going to consolodate them to make plotting and nomenclature a bit easier.  We can give a 
### clearer description of how far apart they are in the M and M...
exp.sites <- exp.sites[c(1,2,4,6,8,11),]
#exp.sites$simple.names <- c("N3","N2","N1","S1","S2","S3")
exp.sites$simple.names <- c("NB","NA","NM","SU","SR","ST")
#exp.sites$exp.type <- c("SnR","CG, SnR", "CG", "CG", "CG, SnR", "SnR")

### get map coordinates
map.bb=c(left=min(clim_exp$lon)-3, bottom=min(clim_exp$lat)-0.5,right=max(clim_exp$lon)+1.7, top=max(clim_exp$lat)+0.5)
map.s <- c(left=12.6, bottom=55.3, right=14.5, top=56.2)
map.n <- c(left=18, bottom=62.75, right=18.6, top=63)

## get map files from google
apikey <- Sys.getenv("GOOGLE_API_KEY") 
register_google(key = apikey)
mapImageData1 <- get_stamenmap(map.bb, color = "bw",source = "stamen", maptype = "toner-lite", zoom=6, api_key=apikey, crop=T, force=T)
mapImageS <- get_stamenmap(map.s, color = "bw", source="stamen", maptype = "toner-lite", zoom=8,api_key=apikey, crop=T, force=T)
mapImageN <- get_stamenmap(map.n, color = "bw",source = "stamen", maptype = "toner-lite", zoom=8, api_key=apikey, crop=T, force=T)
  

## make maps
map=ggmap(mapImageData1, darken=0.1)+geom_point(aes(x=lon, y=lat), data=clim, cex=1.2, col="dodgerblue4", alpha=0.5)+geom_point(aes(x=lon, y=lat), data=exp.sites, col="firebrick2", alpha=0.5, cex=2.5)+labs(title="A", x="longitude", y="latitude")+theme(plot.title = element_text(hjust=0))
map <- map + annotate("rect", xmin=map.s[1], xmax=map.s[3], ymin=map.s[2], ymax=map.s[4], alpha=0, color="black", size=0.25)
map <- map + annotate("rect", xmin=map.n[1], xmax=map.n[3], ymin=map.n[2], ymax=map.n[4], alpha=0, color="black", size=0.25)

mapN=ggmap(mapImageN, darken=0.1) + geom_point(aes(x=lon, y=lat), data=clim, cex=2, col="dodgerblue4", alpha=0.5)+ geom_point(aes(x=lon, y=lat), data=exp.sites, cex=6, col="firebrick2", alpha=1) + labs(title="B", x="longitude", y="latitude")+theme(plot.title = element_text(hjust=0))
mapN <- mapN + geom_text(aes(x=lon, y=lat, label=simple.names), data=exp.sites, fontface="bold", size=3.5)

mapS=ggmap(mapImageS, darken=0.1) + geom_point(aes(x=lon, y=lat), data=clim, cex=2, col="dodgerblue4", alpha=0.5)+ geom_point(aes(x=lon, y=lat), data=exp.sites, cex=6, col="firebrick2", alpha=1) + labs(title="C", x="longitude", y="latitude")+theme(plot.title = element_text(hjust=0))
mapS <- mapS + geom_text(aes(x=lon, y=lat, label=simple.names), data=exp.sites, fontface="bold", size=3.5 )
#mapS <- mapS + annotate("text", x = exp.sites$lon, y = exp.sites$lat, label = exp.sites$exp.type , color="orange", size=3, hjust=2)

#pdf("./figures/map_exp.DLF.2OCT19.pdf" , paper="special", width=7, height=5, pointsize=6)
pdf("./figures/map_exp.DLF.21JUN20.pdf" , paper="special", width=6, height=5, pointsize=6)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(map, vp = vplayout(1:2, 1))
print(mapN, vp = vplayout(1, 2))
print(mapS, vp = vplayout(2, 2))
dev.off()


## ---- end-of-map
