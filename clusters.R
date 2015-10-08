source("makeTraces.R")
library(ggplot2)

numPconst       <- 4000
max.steps.const <- 5
numTraces       <- 1000 #adjust this based on how many samples you would like
clusterSize     <- 3 #adjust this based on results of plotCluster()


# clustering methods
makeClusters <- function(numRep = 100, numP = numPconst, max.steps = max.steps.const, correction.offset=10, updateProgress = NULL, lowG = .1, highG = 1.3, sensitivityDivisor = 80 )
{
  cat("correction offset is ", correction.offset)
  
  #preallocate vectors to hold scalar measurements for each event
  #this is faster than extending vector many times
  mean.G.step <- rep(NA, max.steps*numRep)
  length.step <- rep(NA, max.steps*numRep)
  noise.G     <- rep(NA, max.steps*numRep)
  trace.array <- rep(NA, max.steps*numRep)
  moving.sum  <- 1
  r <- makeTrace()
  lengthT <- max(r$Displacement)
  deltaX  <- lengthT/numP
  df2     <- data.frame(Displacement=numeric(), Conductance = numeric())
  
  #iterating through each break junction pull
  for(i in 1:numRep){
    #simulating a break junction and identifying all 'events' within conductance plateau
    trace <- makeTrace()
    jumps <- findJumps(highG, lowG, trace, sensitivityDivisor)
    jumps <- sapply(jumps, function(x) x2point(trace, x))
    #iterating through each conductance event 
    for(j in 2:length(jumps)){
      #the jumps identify the start and end of each step within conductance plateau
      end   <- jumps[j] - correction.offset
      start <- jumps[j-1]  + correction.offset
      cat("start is ", start, " and end is ", end)
      #compute relevant scalar qualities of conductance traces
      length.step[moving.sum] <- (end - start)
      mean.G.step[moving.sum] <- mean(r$Conductance[start:end])
      noise.G[moving.sum]     <- sd(r$Conductance[start:end])
      trace.array[moving.sum] <- paste(r$Conductance[start:end], collapse = ",")
      
      #save vectorized time serie for displacement and conductance into a data frame where all time series are compiled without looking to individual trace identity
      #this is used for making 2d plots of time series from end of step backwards
      displacement            <- rev(c(1:length(r$Conductance[start:end])))*-1*deltaX
      df3                     <- data.frame(Displacement=displacement, Conductance = r$Conductance[start:end])
      
      #df2 holds all conductance-displacement pairs but does not retain info about individual step identity 
      df2                     <- rbind(df2, df3)      
      moving.sum              <- moving.sum + 1 
    }
    
    if(is.function(updateProgress)){
       if(i%%(round(numRep/4))==0) {
         updateProgress(detail = "generating data", value=((i/round(numRep/4))*.5))
       }
    }
    
  }	

  #cut off na's from preallocated trace.array
  trace.array <- trace.array[!is.na(length.step)]
  mean.G.step <- mean.G.step[!is.na(length.step)]
  noise.G     <- noise.G[!is.na(length.step)]
  length.step <- length.step[!is.na(length.step)]
  
  
  
  #df1 tracks all scalar qualities of traces and also the conductance step time series, which is saved as a string
  df1 <- data.frame(Noise = noise.G, AvgC = mean.G.step, Length=length.step, Trace = I(trace.array), check.names = FALSE, row.names = NULL) 
  
  list(df1, df2)
}


#plot within sum of squares of k means fits to number of groups to assess optimal group
#rule of thumb is to look for an "elbow" in the data where contribution to fit of each group goes down
#this function comes from http://www.statmethods.net/advstats/cluster.html
plotClusters <-function(df)
{
  #scale data since units and types of scalar quantities for clustering are entirely different
  #they are all equally weighted, but this is one thing that could be changed
  #for example weigh more reliable quantities more
  mydata <- scale(df)
  for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                       centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
}




bestCluster <-function(df, k, df2, updateProgress = NULL){

  
  if(is.function(updateProgress)){
      updateProgress(detail = "computing clusters", value=.6)
  }
  fit    <- kmeans(x=scale(df), centers=k) 
  updateProgress(detail = "done computing clusters", value=.8)
  mydata <- data.frame(df2, Cluster = fit$cluster)
}


#FIND PLACES WITHIN SINGLE CONDUCTANCE PLATEAU WHERE THERE ARE SIGNIFICANT CHANGES/JUMPS BETWEEN CONFIGURATIONS
findJumps2 <- function(beginG, endG, trace, lookahead = 20, sensitivityDivisor = 100, loq=0){
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
  
}
