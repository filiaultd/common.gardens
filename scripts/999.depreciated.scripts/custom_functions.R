

##a function to calculate a mean, ignoring NAs, but requiring at least 3 reps.
betterMean=function(x, N=3){x=na.omit(x);if(length(x)>=N){return(mean(x))}else{return(NA)}}
