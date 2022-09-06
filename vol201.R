
parameters <- c(kt = 1/2.32, qc = 160, Vc = 469, q1 = 746, V1 = 3986, q2 = 190, V2 = 15557)
state <- c(D = 0, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times <- seq(1,191,by = 1)
out <- ode(state, times, odefunction, parameters)
drug <- data.frame(out)
state1 <- c(D = 480, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times1 <- seq(192,551,by = 1)
out1 <- ode(state1, times1, odefunction, parameters)
drug1 <- data.frame(out1)
dim <- dim(drug1)[1]
state2 <- c(D=drug1[dim,2]+960, T1=drug1[dim,3],T2=drug1[dim,4],C=drug1[dim,5],P1=drug1[dim,6],P2=drug1[dim,7])
times2 <- seq(552,720,by=1)
out2 <- ode(state2, times2, odefunction, parameters)
drug2 <- data.frame(out2)
pk201 <- rbind(drug,drug1,drug2)

vol201 <- functional_prob(data201, pk201, 720)

p <- vector()
for (i in 1:1000) {
  p[i] <- sum(vol201[i,]<Inf)/100
}
hist(p, breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), main = "Volunteer 201",cex.main=2.5,cex.axis=2,cex.lab=2,ylim=c(0,1000))

probmat201 <- matrix(0,1000,720)
for (i in 1:1000){
  for (j in 1:720){
    probmat201[i,j] <- length(which(vol201[i,]<=j))/100
  }
}

med <-vector() 

for (i in 1:720){
  med[i] <- quantile(probmat201[,i],0.5)
}

top95 <-vector() 
for (i in 1:720){
  top95[i] <- quantile(probmat201[,i],0.975)
}

bot95 <-vector() 
for (i in 1:720){
  bot95[i] <- quantile(probmat201[,i],0.025)
}
top25 <-vector() 
for (i in 1:720){
  top25[i] <- quantile(probmat201[,i],0.625)
}

bot25 <-vector() 
for (i in 1:720){
  bot25[i] <- quantile(probmat201[,i],0.375)
}

top50 <-vector() 
for (i in 1:720){
  top50[i] <- quantile(probmat201[,i],0.75)
}

bot50 <-vector() 
for (i in 1:720){
  bot50[i] <- quantile(probmat201[,i],0.25)
}

top75 <-vector() 
for (i in 1:720){
  top75[i] <- quantile(probmat201[,i],0.875)
}

bot75 <-vector() 
for (i in 1:720){
  bot75[i] <- quantile(probmat201[,i],0.125)
}

t <- 1:720/24

plot(t,med,ylim=c(0,1),xlim = c(5,30), type="n",xlab="Time(Day)",ylab="Extinction Probability",main="Volunteer 201",cex.main=2.5,cex.axis=2,cex.lab=2)

polygon(c(t, rev(t)), c(top95, rev(bot95)),
        col = b1, border = NA)
polygon(c(t, rev(t)), c(top75, rev(bot75)),
        col = b2, border = NA)
polygon(c(t, rev(t)), c(top50, rev(bot50)),
        col = b3, border = NA)
polygon(c(t, rev(t)), c(top25, rev(bot25)),
        col = b4, border = NA)


lines(t,med,type="l")
lines(t,top95, type="l",lty=2,col=b1)
lines(t,bot95, type="l",lty=2,col=b1)
lines(t,top75, type="l",lty=2,col=b2)
lines(t,bot75, type="l",lty=2,col=b2)
lines(t,top50, type="l",lty=2,col=b3)
lines(t,bot50, type="l",lty=2,col=b3)
lines(t,top25, type="l",lty=2,col=b4)
lines(t,bot25, type="l",lty=2,col=b4)

arrows(8,0.2, 8, 0, col = 'black',lwd=3,length=0.1)
text(8, 0.25, "1st")
arrows(23,0.2, 23, 0, col = 'black',lwd=3, length = 0.1)
text(23, 0.25, "2nd")

