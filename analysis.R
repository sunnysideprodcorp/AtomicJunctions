source("makeTraces.R")
library(ggplot2)

#FIND CROSSINGS OF A PARTICULAR CONDUCTANCE VALUE WITHOUT CHANGING DISPLACEMENT OF TRACE
findG <-function(alignG, trace)
{
  trace <- trace - alignG
  f <- which(trace>0, arr.ind=TRUE)
  if (length(f)>0) {  
    f[length(f)]
  }
  else{
    0
  }
}

findJumps2 <- function(){
x = c(1:100)
plot(x)
}
  
#FIND PLACES WITHIN SINGLE CONDUCTANCE PLATEAU WHERE THERE ARE SIGNIFICANT CHANGES/JUMPS BETWEEN CONFIGURATIONS
findJumps <- function(beginG=.3, endG=1.5, trace, lookahead = 20, sensitivityDivisor = 100, loq=0){

  if(is.na(beginG)){
    beginG = .3
  }
  if(is.na(endG)){
    endG = 1.5
  }
  
    #find beginning and end of conductance plateau
  startP <- findG(beginG, trace$Conductance)
  endP <- findG(endG, trace$Conductance)
  
  #select only portion of conductance signature corresponding to conductance plateau
  traceLocal <- subset(trace, Displacement > point2X(trace, startP - lookahead) & Displacement < point2X(trace, endP + lookahead), select= c(Displacement, Conductance))
  y <- traceLocal$Conductance
  
  #identify significant events by looking at sharp changes in difference of look-behind and look-ahead moving averages
  for(i in lookahead:(length(traceLocal$Conductance)-lookahead)){
    y[i] <- mean(traceLocal$Conductance[(i-lookahead):i]) - mean(traceLocal$Conductance[i:(i+lookahead)])
  }
  y <- y*y
  x <- findLevels(y, (max(y, na.omit<-TRUE)/sensitivityDivisor))
  
  #move crossing points to adjust back to full time series and not subset used for analysis in for-loop
  x <- x + startP - lookahead
  compressedLevels <- compressLevels(startP, endP, lookahead, x)
  
  compressedLevels <- c(startP, compressedLevels, endP)   
  compressedLevels = sapply(compressedLevels, function(x) point2X( trace, x))
  
  #if in loquacious mode, show fittings to see whether sensitivity is properly adjusted to find individual steps
  if(loq == 1)
  {
    #plot(trace$Conductance, type = "l", xlim=c(startP-2*lookahead, endP+2*lookahead))
    plot(trace, type = "l")

  for(i in 1:length(compressedLevels))
  {
    abline(v = compressedLevels[i], col="blue")
  }
  
  abline(h=beginG, col="red")
  abline(h=endG, col="red")
  
  
 }
 
  compressedLevels
}

#DETERMINE ALL POINTS WHERE A PUTATIVELY CONTINUOUS TIME SERIES PASSES A VALUE EITHER RISING OR DESCENDING
findLevels <- function(trace, level)
{
  trace2 <- trace
  trace2 <- trace2-level
  x <- which(diff(sign(trace2))!=0)  
}

#REMOVES REDUNDANT MARKS OF SIGNIFICANT EVENTS BECAUSE CURRENT METHODOLOGY IDENTIFIES TWO POINTS FOR EACH EVENT NOT ONE
compressLevels <- function(startP, endP, lookahead, levelsV)
{
  #remove values that are too close to the ends and so probably associated with beginning and end of full plateau 
  levelsV <- levelsV[which(!(abs(levelsV - startP)<(2*lookahead)))]  
  levelsV <- levelsV[which(!(abs(levelsV - endP)<(2*lookahead)))]
  
  #remove values that are "doubles", that is two indicators so close they indicate the same event
  levelsVlag <- c(0, levelsV)
  levelsV <- c(levelsV, 0)
  lagV <- levelsV - levelsVlag  
  levelsV <- levelsV[which(abs(lagV)<2*lookahead)]
  levelsV <- levelsV[which(levelsV !=0)]
}
