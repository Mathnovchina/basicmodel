library(tikzDevice) # load the tikzDevice package
library(ggplot2)


source("SourceCode.R")


sysMinDist <- cppMakeSys("introduction.R", reportVars = 3)
resSimul <- cppRK4(sysMinDist)

View(resSimul[,c("")])

par(mfrow = c(2,2), mar=c(2.1, 5.1, 4.1, 4.1), xpd=T)
plot(x = resSimul$time+2016, y = resSimul$lambda, type = 'l', col = 'blue', xlab = '', ylab = '', ylim = c(0,1))
lines(x = resSimul$time+2016, y = resSimul$omega, type = 'l', col = 'royalblue')
lines(x = resSimul$time+2016, y = resSimul$pi, type = 'l', col = 'purple')
#lines(x = resSimul$time+2016, y = resSimul$n, type = 'l', col = 'magenta')
lines(x = resSimul$time+2016, y = resSimul$dam, type = 'l', col = 'red')
lines(x = resSimul$time+2016, y = (resSimul$Div/(resSimul$p*resSimul$Y)), type = 'l', col = 'orange')
legend("topright", inset=c(-.2,0), legend = c("Employment", "Wage", "Profit",  "Damage", "Dividend"), title = "Rates", col=c("blue", "royalblue", "purple", "red", "orange"), lty=1:1, box.lty=0, cex=0.75)
#"GHG reduction", , "magenta"
plot(x = resSimul$time+2016, y = resSimul$Y_0,type = 'l', col = 'blue', xlab = '', ylab="Trillion $")
lines(x = resSimul$time+2016, y = resSimul$Y, type = 'l', col = 'red')
legend("topright", inset=c(-.2,0), legend = c("Y0","Y"), title = "GDP", col=c("blue", "red"), lty=1:1, box.lty=0, cex=0.75, )

plot(x = resSimul$time+2016, y = resSimul$d     , type = 'l', col = 'blue', xlab = '', ylab = '', ylim = c(-2,6))
lines(x = resSimul$time+2016, y = resSimul$i*100     , type = 'l', col = 'purple', xlab = '', ylab = '')
legend("bottomright", inset=c(-.2,0), legend = c("debt/GDP","inflation (%)"), col=c("blue",'purple'), lty=1:1, box.lty=0, cex=0.75)

plot(x = resSimul$time+2016, y = resSimul$E, type = 'l', col = 'brown', xlab = '', ylab = 'Gtc')
legend("bottomright", inset=c(-.2,0), legend = c("Emissions"), col=c("brown"), lty=1:1, box.lty=0, cex=0.75)

# Y_p = resSimul$Y
# emp_p = resSimul$lambda
# w_p = resSimul$omega
# d_p = resSimul$d
# i_p = resSimul$i
# F_tro_p = resSimul$F_tro
# F_bor_p = resSimul$F_bor
# F_tem_p = resSimul$F_tem
# H_tro_p = resSimul$H_tro
# H_tem_p = resSimul$H_tem
# H_bor_p = resSimul$H_bor
# FO_c_p = resSimul$FO_c
# BI_m3_p = resSimul$BI_m3
# tax_p = resSimul$tax
# HB_tro_p = resSimul$HB_tro
# HB_tem_p = resSimul$HB_tem
# HB_bor_p = resSimul$HB_bor
# E_ant_p = resSimul$E_ant
# E_p = resSimul$E
# Temp_p = resSimul$Temp
#----------------------------------------------------

par(mfrow = c(3,3), mar=c(2.1, 5.1, 4.1, 4.1), xpd=T)

plot(x = resSimul$time+2016, y = resSimul$Y, type = 'l', col = 'black', xlab = '', ylab = '',ylim = c(50,450), main = "Output (2015 US$ tril.)")
line(x = resSimul$time+2016, y = Y_p, type = 'l', col = 'black', lty = 2)

plot(x = resSimul$time+2016, y = emp_p, type = 'l', col = 'blue',lty = 2, xlab = '', ylab = '', ylim = c(0.4,0.8), main="Labour")
lines(x = resSimul$time+2016, y = w_p, type = 'l', col = 'olivedrab', lty = 2)
lines(x = resSimul$time+2016, y = resSimul$lambda, type = 'l', col = 'blue', )
lines(x = resSimul$time+2016, y = resSimul$omega, type = 'l',  col = 'olivedrab')
axis(side = 4)  
mtext("Employment rate (blue)", side = 2, line = 3)
mtext("Wage share (%, green)", side = 4, line = 3)

plot(x = resSimul$time+2016, y = d_p    , type = 'l', lty = 2, col = 'blue',ylim=c(1,5), xlab = '', ylab = '', main = "Financial stability")
lines(x = resSimul$time+2016, y = resSimul$d, type = 'l', col = 'blue')
par(new=TRUE)
plot(x = resSimul$time+2016, y = resSimul$i*100, type = 'l', col = 'olivedrab',ylab = '', ylim = c(-1,5),axes = FALSE)
lines(x = resSimul$time+2016, y = i_p*100    , type = 'l',lty = 2, col = 'olivedrab',xlab = '', ylab ='')
axis(side = 4)     
mtext("private debt ratio (blue)", side = 2, line = 3)
mtext("Inflation (%, green)", side = 4, line = 3)

plot(x = resSimul$time+2016, y = resSimul$F_tro, type = 'l', col = 'goldenrod', xlab = '', ylab = '',ylim = c(50,500),main ="Biomass stocks (billion m3)" )
lines(x = resSimul$time+2016, y = resSimul$F_bor, type = 'l', col = 'darkolivegreen4')
lines(x = resSimul$time+2016, y = resSimul$F_tem, type = 'l', col = 'cornflowerblue' )
lines(x = resSimul$time+2016, y = F_tro_p, type = 'l', col = 'goldenrod', lty = 2)
lines(x = resSimul$time+2016, y = F_bor_p, type = 'l', col = 'darkolivegreen4', lty = 2)
lines(x = resSimul$time+2016, y = F_tem_p, type = 'l', col = 'cornflowerblue', lty = 2)

plot(x = resSimul$time+2016, y = resSimul$H_tro, type = 'l', col = 'goldenrod',xlab = '', ylab = '', ylim = c(0.5,10), main = "Total harvest (billion m3)")
lines(x = resSimul$time+2016, y = resSimul$H_bor, type = 'l', col = 'darkolivegreen4')
lines(x = resSimul$time+2016, y = resSimul$H_tem, type = 'l', col = 'cornflowerblue' )
lines(x = resSimul$time+2016, y = H_tro_p, type = 'l', col = 'goldenrod', lty = 2)
lines(x = resSimul$time+2016, y = H_bor_p, type = 'l', col = 'darkolivegreen4', lty = 2)
lines(x = resSimul$time+2016, y = H_tem_p, type = 'l', col = 'cornflowerblue', lty = 2)

plot(x = resSimul$time+2016, y = resSimul$FO_c , type = 'l', col = 'orange4', xlab = '', ylab = '',ylim = c(0,80), main = "Energetic demand")
lines(x = resSimul$time+2016, y = resSimul$BI_m3, type = 'l', col = 'olivedrab',ylab = '',axes = FALSE)
lines(x = resSimul$time+2016, y = FO_c_p, type = 'l', col = 'orange4',ylab = '',axes = FALSE, lty = 2)
lines(x = resSimul$time+2016, y = BI_m3_p, type = 'l', col = 'olivedrab',ylab = '',axes = FALSE, lty = 2)
axis(side=4)
mtext("Fossil (GtC, brown)", side = 2, line = 3, font = 0.5)
mtext("Bioenergy (billion m3, olive)", side = 4, line = 3, font = 0.5, cex = 0.9)

plot(x = resSimul$time+2016, y = resSimul$tax, type = 'l', col = 'black', xlab = '', ylab = '',ylim = c(0,300), main ="Carbon policy ($/tonC)" )
lines(x = resSimul$time+2016, y = tax_p, type = 'l', col = 'black', lty = 2)

plot(x = resSimul$time+2016, y = resSimul$Temp    , type = 'l', col = 'black',xlab = '', ylab ='', ylim=c(0.5,4.5), main = "Temperature anomaly (°C)")
lines(x = resSimul$time+2016, y = Temp_p, type = 'l', lty = 2, col = 'black', )


plot(x = resSimul$time+2016, y = resSimul$E_ant,type = 'l', col = 'orange4', xlab = '', ylab = '',ylim=c(0,90),main = "Emissions (GtCO2-e)")
lines(x = resSimul$time+2016, y = resSimul$E, type = 'l', col = 'olivedrab', )
lines(x = resSimul$time+2016, y = E_ant_p,type = 'l', col = 'orange4', lty = 2 )
lines(x = resSimul$time+2016, y = E_p, type = 'l', col = 'olivedrab', lty = 2)
axis(side=4)
mtext("Gross (brown)", side = 2, line = 3, font = 0.5)
mtext("- sequestration (olive)", side = 4, line = 3, font = 0.5,cex = 0.9)





plot(x = resSimul$time+2016, y = E_n, type = 'l', col = 'blue', xlab = '', ylab = '')
lines(x = resSimul$time+2016, y = resSimul$E, type = 'l', lty = 4, col = 'blue', )
axis(side = 4)     
mtext("Emissions (GtCO2e)", side = 4, line = 3, font = 2)

plot(x = resSimul$time+2016, y = T_n    , type = 'l', col = 'blue',xlab = '', ylab ='', ylim=c(0,5))
lines(x = resSimul$time+2016, y = resSimul$Temp, type = 'l', lty = 4, col = 'blue', )
axis(side = 4)     
mtext("Temperature anomaly (°C)", side = 4, line = 3, font = 2, cex = 0.8)
View(resSimul)

plot(x = resSimul$time+2016, y = resSimul$FO_c, type = 'l', col = 'brown', xlab = '', ylab = '',ylim=c(0,100))
lines(x = resSimul$time+2016, y = resSimul$BI_m3, type = 'l', col = 'olivedrab', )
mtext("Bioenergy (Billion m3, green)", side = 4, line = 1, font = 0.5)


mtext("Fossil fuel (GtC, brown)", side = 2, line = 3, font = 0.5)
mtext("Bioenergy (Billion m3, green)", side = 4, line = 1, font = 0.5)
