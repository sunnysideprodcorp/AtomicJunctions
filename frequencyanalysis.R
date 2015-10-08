require(audio)
require(tuneR)



plot_freq_dep<-function(trace, beginG = 1.3, endG = .3)
{
  beginPoint = findG(beginG, trace$Conductance) + 100
  endPoint = findG(endG, trace$Conductance) - 100
  fft_r = fft(trace$Conductance[beginPoint:endPoint])
  ylim = quantile(Re(fft_r), c(.05, .99))
  cat("ylim is ", ylim)
  par(mfrow=c(3,1))
  cat("lengths are ", length(c(beginPoint:endPoint)), length(fft_r))
  plot(abs(Re(fft_r)), ylim=ylim, xlim = c(1,(length(fft_r)/2)))
  plot(c(beginPoint:endPoint), trace$Conductance[beginPoint:endPoint])
  plot(c(1:length(trace$Conductance)), trace$Conductance)
   cat("length of fft is ", length(fft_r))
  cat("begin point ", beginPoint, " and end point ", endPoint)
}


r = makeTrace(lowT=0)
plot_freq_dep(r)
r = makeTrace(lowT=1)
plot_freq_dep(r)

plot_freq_dep(makeTrace())
require(audio)
require(tuneR)

w <- noise(kind = c("white"))
par(mfrow=c(2,1))
plot(w,main="white noise")
plot(p,main="pink noise")
writeWave(p,"p.wav")#writes pink noise on your hard drive
require(audio)#loads `audio` package to use `load.wave` function
p.vec <- load.wave("path/to/p.wav")#this will load pink noise as a vector