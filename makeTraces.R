
require(audio)
require(tuneR)




#MAKE STEPS WITHIN A CONDUCTANCE STEP

makeSteps <- function(lowG, highG, num.G.groups = 3, max.steps = max.steps.const, mean.step.p = 200, sd.step.p = 50, temp=0, lowT=1){
  
  #decide how many within conductance plateau events there are for this plateau and segment conductance into boxes
  mid               <- rnorm(2, mean = (lowG + highG)/2, sd = (lowG + highG)/5)
  conductance.boxes <- seq(from=min(mid), to=max(mid), length.out = num.G.groups)
  numSteps          <- sample(1:max.steps, 1)
  
  #choose conductance and step length for each step
  #enforce G values that are physically realistic, no higher or lower than cutoff values
  stepsG <- rnorm(numSteps, mean = (lowG + highG)/2, sd = (lowG + highG)/5)
  stepsG <- replace(stepsG, stepsG<lowG, lowG)
  stepsG <- replace(stepsG, stepsG>highG, highG)
  stepsL <- floor(abs(rnorm(n = numSteps, mean = mean.step.p, sd=sd.step.p)))
  
  
  trace <- c(1:1000)
  meanG = mean(c(lowG, highG))
  trace = rnorm(n=1000, mean = meanG, sd = meanG/10 )
  
  #prepare entire molecular step
  if(lowT){
  trace  <- c(1:sum(stepsL))
  startp <- 1
  for(i in 1:length(stepsL)){
    endp               <- startp + stepsL[i] - 1
    trace[startp:endp] <- rnorm(n=stepsL[i], mean = stepsG[i], sd = stepsG[i]/12)
    startp             <- startp + stepsL[i]
  }
  }

  
  trace
}

# MAKES TRACES WITH A MOLECULAR STEP
makeTrace <- function(molG=10^-4, molSD=10^-4, excursion=5, GDist=.5, numP=numPconst, multiples=FALSE, temp=0, shiny=0, lowT=1){
  
  #make 3,2,1  G0 steps that have within-plateau events if temp = 0, otherwise they are relatively flat
  g3 <- makeSteps(lowG=2.1, highG=2.5, temp=temp, lowT=lowT)
  g2 <- makeSteps(lowG=1.4, highG=min(g3), temp=temp, lowT=lowT)
  g1 <- makeSteps(lowG=.5, highG=min(g2), temp=temp, lowT=lowT)
  
  
  #length measured in points
  g1stepL <- length(g1)
  g2stepL <- length(g2)
  g3stepL <- length(g3)
  
  #conversion between nm and indices
  delX <- numP/excursion
  dt   <- seq(from = 0 , to = excursion, excursion/numP)
  
  #store conductance time series in g (G for conductance)
  g      <- c(1:numP + 1)
  offset <- 0
  
  #assign 3,2,1 G0 steps generated earlier individually to single time series
  g[(1 + offset) : (g3stepL + offset)] <- g3 
  offset                               <- g3stepL + offset
  
  g[(1 + offset) : (g2stepL + offset)] <- g2 
  offset                               <- g2stepL + offset
  
  g[(1 + offset) : (g1stepL + offset)] <- g1 
  offset                               <- g1stepL + offset
  
  #metallic junction breaks to electronics noise level
  g[(offset):numP + 1] <- rnorm(n=(numP-offset+1),mean=10^-5, sd = 10^-5)
  
  #add 1/f noise aka pink noise
  p <- noise(kind = c("pink"))
  g <- g + p@left[1:length(g)]/10
  
  
  #subtract out negative noise so all values positive
  m <- min(g)
  g <- g - m +.000001
  g <- remove.repeat.NA(trace.smooth(g, type="moving-average", width=5))
  df <- data.frame(Displacement=dt,Conductance=g, check.names = FALSE, row.names = NULL)  
}

#REMOVES NA ONCE TRACE HAS BEEN REALIGNED
#from stackoverflow http://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value
remove.repeat.NA <- function(x) {   # repeats the last non NA value. Keeps leading NA
  ind <- which(!is.na(x))         # get positions of nonmissing values
  
  if(is.na(x[1]))                 # if it begins with a missing, add the 
    ind <- c(1,ind)               # first position to the indices
  
  t   <- diff(c(ind, length(x) + 1)) 
  vec <- rep(x[ind], times = t) # diffing the indices + length yields how often 
  # they need to be repeated
  
  tt     <- tail(t, n=1)
  ind    <- which(!is.na(vec))
  new    <- abs(diff(c(length(vec)+1, ind)))
  new[1] <- ind[1]
  vec2   <- rep(vec[ind], times = new)
  
}


#SMOOTHES TRACES
#from stack overflow http://stackoverflow.com/questions/9988292/how-to-smooth-a-curve
trace.smooth<-function(trace, type="Savitsky-Golay", width=10){
  
  if(type=="lowess"){
    smooth.trace<-with(trace, lowess(x=1:length(trace),
                                     y=trace,
                                     f=width/length(trace),
                                     delta=width/2))$y
  }
  
  if(type=="moving-average"){
    moving_average<-function(width=10){
      moving.average<-rep(1,width)/width
      return(moving.average)
    }
    
    moving.average<-moving_average(width)
    
    smooth.trace<-filter(trace, moving.average)
    
  }
  
  if(type=="Savitsky-Golay"){
    # Savitsky-Golay smoothing function 
    savistsky_golay<-function(width=10){
      x<-1:width-width/2
      y<-max(x^2)-x^2
      sg<-y/sum(y)
      return(sg)
    }
    
    sg<-savistsky_golay(width)
    
    smooth.trace<-filter(trace, sg)
  }
  
  return(smooth.trace)
}




#FINDS PLACE WHERE TRACE CROSSES ALIGNG
alignG <-function(alignG, trace, loq=0)
{
  trace2 <- trace$Conductance - alignG
  f <- which(trace2>0, arr.ind=TRUE)
  offset <- trace$Displacement[f[length(f)]] 
  trace$Displacement <-   trace$Displacement - offset
  if(loq==1){
    plot(trace)
    abline(h=alignG)
    abline(v=0)
  }  
  trace
}



#tests alignment
showAlignment <-function(Gval = .5){
  r <- makeTrace()
  index <- alignG(Gval, r$Conductance)
  plot(r)
  abline(h=Gval, col=4, lty=3)
  abline(v= point2X(r, index), col=4,lty=3)
}


#converts row number to displacement
point2X <-function(trace, point){
  displacement_range = max(trace$Displacement) - min(trace$Displacement)
  point/(length(trace$Displacement)-1)*displacement_range + min(trace$Displacement)
}

#converts displacement to row number
x2point <-function(trace, x){
  displacement_range = max(trace$Displacement) - min(trace$Displacement)
  floor( (x - min(trace$Displacement))/displacement_range * length(trace$Displacement))
}



#MAKE 2D HISTOGRAM
make2DHist <- function(N=50, alignG = .2, log=TRUE, method=1, createNew = TRUE, df=na, lowT=1, updateProgress = NULL)
{

  if(is.function(updateProgress)){
      updateProgress(detail = "generating data", value=.5)
  }
  if(createNew)
  {
    v  <- replicate(N, alignG(alignG, makeTrace(lowT=lowT)), simplify=FALSE)
    df <- do.call(rbind, v)
  }

  if(is.function(updateProgress)){
    updateProgress(detail = "compiling data", value=.5)
  }
  
  if (method == 1){
    #first method uses ggplot
    d <- ggplot(df,aes(Displacement, Conductance)) + xlim(-3,.5)
    d + stat_bin2d(bins=500)	
  }
  else{
    if(method==2){
      #second method using KernSmooth
      z <- bkde2D(df, .5)
      persp(z$fhat)
    }
    else{
      #third method used plain old hist2d
      xbins <- .1
      ybins <- .1
      
      y <- hist2d(df, nbins = c(100,100))
    }
  }
  
}


