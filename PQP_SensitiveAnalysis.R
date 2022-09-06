#PK model used for simulations
parameters <- c(kt = 1/2.79, qc = 51.4, Vc = 804, q1 = 2020, V1 = 3010, q2 = 149, V2 = 13300)
state <- c(D = 0, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times <- seq(1,167,by = 1)
out <- ode(state, times, odefunction, parameters)
drug <- data.frame(out)
state1 <- c(D = 480, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times1 <- seq(168,720,by = 1)
out1 <- ode(state1, times1, odefunction, parameters)
drug1 <- data.frame(out1)
dr <- rbind(drug,drug1)

#choose range of parameters
kmax <- 10:30/100
ec50 <- 40:180/10
rec <- matrix(0,1,3)
for (i in 1:length(kmax)) {
  for (j in 1:length(ec50)) {
    killing <- kd(dr,kmax[i],ec50[j],12)
    record <- vector()
    for (k in 1:100){
      simulation <- try(cycle(350,7.5,5,80,killing,0.004,0.03,al=42,as=25,720),silent = TRUE)
      record[k] <- ifelse(rowSums(simulation)[which.min(rowSums(simulation))]==0,which.min(rowSums(simulation)),Inf)
    }
    probability <- sum(record!=Inf)/100
    k.e.p <- as.vector(c(kmax[i],ec50[j],probability))
    rec <- rbind(rec,k.e.p)
    }
  }
rec <- data.frame(rec[-1,])

sensitivity <- rec
attach(sensitivity)
colnames(sensitivity) <- c("kmax","EC50","Probability")


library(ggplot2)
a <- ggplot(sensitivity, aes(x=kmax,y=EC50))
a <- a + geom_tile(aes(fill=Probability)) +labs(x=expression(k[max]),y=expression(EC[50]),cex=4)
a
