library(tikzDevice) # load the tikzDevice package
library(ggplot2)

source("SourceCode.R")
sysMinDist <- cppMakeSys("introduction.R", reportVars = 3)
resSimul <- cppRK4(sysMinDist)

View(resSimul[,c("")])


par(mfrow = c(2,3), mar=c(2.1, 5.1, 4.1, 4.1), xpd=T)

plot(x = resSimul$time+2016, y = resSimul$Y, type = 'l', col = 'black', xlab = '', ylab = '', main = "Output")
plot(x = resSimul$time+2016, y = resSimul$Tax, type = 'l', col = 'black', lty = 2,main = "Tax")
plot(x = resSimul$time+2016, y = resSimul$Yd, type = 'l', col = 'black', lty = 2,main = "Yd")
plot(x = resSimul$time+2016, y = resSimul$C, type = 'l', col = 'black', lty = 2,main = "Consumption")
plot(x = resSimul$time+2016, y = resSimul$S, type = 'l', col = 'black', lty = 2, main = "Savings/bonds")

