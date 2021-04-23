#### R functions for geo analysis


#############################################
####### prep.SAF - prepare starting genotype matrix
################################################

prep.SAF <- function(all.gts.file, min.MAF, max.NAF){
  ### MAF.min is minimum minor allele frequency
  ### NAF.max is maximum NA frequency
  ### 200 is hard coded as starting accession number
  n.accs <- 200
  
  all.gts <-read.csv(all.gts.file,stringsAsFactors=FALSE,colClasses="numeric",header=TRUE)
  ### remove all the resequenced local samples
  all.gts <- all.gts[,1:202]
  
  pos <- all.gts[,1:2] 
  ## this is the data set used for genotyping these samples (from Fernando)
  all.gts <- all.gts[,-c(1:2)]
  colnames(all.gts)<- gsub("X","",colnames(all.gts))
  ugh <- as.matrix(all.gts)
  all.gts <- ugh
  rm(ugh)
  ### coding of this is 0=ref, 1=alt, -1=NA, 2=het
  ### new coding to make 0.5=het and NA=NA
  all.gts[which(all.gts==-1)]<- NA
  all.gts[which(all.gts==2)] <- 0.5
  
  all.gts.alt <- apply(all.gts,1, function(x){sum(x,na.rm=TRUE)})
  all.gts.na <- apply(all.gts,1,function(x){length(x[is.na(x)])})
  #all.gts.sample <- n.accs-all.gts.na
  all.gts.alt.f <- all.gts.alt/n.accs
  all.gts.na.f <- all.gts.na/n.accs
  all.gts.ref.f <- 1-(all.gts.alt.f+all.gts.na.f)
  
  all.gts.f <- cbind(all.gts.ref.f, all.gts.alt.f, all.gts.na.f)
  colnames(all.gts.f) <- c("ref.f","alt.f","na.f")
  
  #### these are starting allele frequencies, including NA positions
  ### remove high NA percentage
  na.rem <- which(all.gts.f[,3] > max.NAF,)
  
  ### remove low MAF
  ## check alleles separately - NAs will make it so that these don't add up to 1
  #check ref allele
  ref.rem <- which(all.gts.f[,1] < min.MAF)
  #check alt allele
  alt.rem <- which(all.gts.f[,2] < min.MAF)
  ##remove these from data
  all.rem <- unique(c(na.rem, ref.rem,alt.rem))
  all.rem <- all.rem[order(all.rem)]
  
  all.gts.filtered <- all.gts[-all.rem,]
  all.pos.filtered <- pos[-all.rem,]
  out.dat <- cbind(all.pos.filtered, all.gts.filtered)
  return(out.dat)
}

###############################################
### following is from Matt Horton ######
### takes data frame with "latitude" and "longitude"
### to get geographic center
###############################################

determineGeographicCenter <- function(population){
  latitudes <- population[,"latitude"] * (pi/180);
  longitudes <- population[,"longitude"] * (pi/180);
  xs <- cos(latitudes) * cos(longitudes);
  ys <- cos(latitudes) * sin(longitudes);
  zs <- sin(latitudes);
  # for clarity (for other projects relying on this code), include weights, although they are of course 1 for this project.
  weights <- rep(1, length(xs));
  x <- sum(xs * weights)/sum(weights);
  y <- sum(ys * weights)/sum(weights);
  z <- sum(zs * weights)/sum(weights);
  longitude <- atan2(y, x);
  hypote_moose <- sqrt(x^2 + y^2);
  latitude <- atan2(z, hypote_moose);
  latitude <- latitude * (180/pi);
  longitude <- longitude * (180/pi);
  return(list("latitude"=latitude, "longitude"=longitude));
}

###############################################
### following is from Matt Horton ######
### lat and long need to be in radians for this one
### haversine Distance between two points
##################################################

haversineDistance <- function(lat1, long1, lat2, long2, radiusOfTheEarth=6371){
  to.radians <- function(measure){measure*(pi/180)}
  lat1 <- to.radians(lat1)
  lat2 <- to.radians(lat2)
  long1 <- to.radians(long1)
  long2 <- to.radians(long2)
  deltaLongitude <- abs(long2 - long1);
  y <- sqrt((cos(lat2) * sin(deltaLongitude))^2 + ((cos(lat1) * sin(lat2)) - (sin(lat1) * cos(lat2) * cos(deltaLongitude)))^2);
  x <- (sin(lat1) * sin(lat2)) + (cos(lat1) * cos(lat2) * cos(deltaLongitude));
  deltaSig <- atan2(y, x);
  d <- radiusOfTheEarth * deltaSig;
  return(d)}


##################################################
### plotting a value across the genome by chromosome
### takes a dataframe with Chromosome and Position columns
### as well as a column what you're trying to plot
###########################################

genome.plot <- function(updata, plot.var){
	chr.lengths <- c(30427671,19698289,23459830,18585056,26975502)
	chr.add <- c(0,cumsum(chr.lengths))[1:5]
	max.bp <- sum(chr.lengths)
	chr.colors <- c("blue","dodgerblue", "blue", "dodgerblue", "blue")
	chr.mids <- chr.add + (chr.lengths/2)
	up.s <- split(updata, updata$Chromosome)
	plot.col <- which(colnames(updata)==plot.var)
	plot(updata$Position,updata[,plot.col], xlim=c(0,max.bp), type="n", xlab="Chromosome", ylab=plot.var, xaxt="n")
	axis(1,at=chr.mids,labels=c(1:5))
	for(up.chr in 1:5){
	  #print(up.chr)
	  up.c <- up.s[[up.chr]]
	  up.add <- chr.add[up.chr]
	  up.c$Position.plot <- up.c$Position + up.add
	  points(up.c$Position.plot, up.c[,plot.col], col=chr.colors[up.chr])
	}
}

#######################################################
#### also a variant that adds points to an already-drawn genome plot
#######################################################

genome.points <- function(updata, plot.var, point.color){
	chr.lengths <- c(30427671,19698289,23459830,18585056,26975502)
	chr.add <- c(0,cumsum(chr.lengths))[1:5]
	max.bp <- sum(chr.lengths)
	up.s <- split(updata, updata$Chromosome)
	plot.col <- which(colnames(updata)==plot.var)
	for(up.chr in 1:5){
	  up.c <- up.s[[up.chr]]
	  up.add <- chr.add[up.chr]
	  up.c$Position.plot <- up.c$Position + up.add
	  points(up.c$Position.plot, up.c[,plot.col], col=point.color, pch=19)
	}
}
