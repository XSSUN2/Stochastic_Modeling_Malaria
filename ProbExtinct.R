library(truncnorm)
record_all <- list()
set.seed(1071174)
#for each parameter, run stochastic model 100 times and record the time that extinction occurs or Inf if recrudescence
for (k in 1:12){
  record <- matrix(NA, 1000,100)
  for (i in 1:1000){
    for(j in 1:100){
      simulation <- try(cycle(k,data[i,2],data[i,3],data[i,4]
                              ,matrix(0,720,3),data[i,8],data[i,9],al=42,as=25,50),silent = TRUE)
      record[i,j] <- ifelse(rowSums(simulation)[which.min(rowSums(simulation))]==0,which.min(rowSums(simulation)),Inf)
    }
  }
  record_all[[k]] <- record
}

#calculate extinction probability at each time point
p <-matrix(0,12,1000)
for (i in 1:12) {
  for (j in 1:1000) {
    p[i,j] <- sum(record_all[[i]][j,]!=Inf)/100
  }
}

top <-vector() 
for (i in 1:12){
  top[i] <- quantile(p[i,],0.975)
}

bot <-vector() 
for (i in 1:12){
  bot[i] <- quantile(p[i,],0.025)
}

med <- vector()
for (i in 1:12){
  med[i] <- median(p[i,])
}

plot(1:12,med,xlab = "Inoculation size",ylab="Extinction probability",type="l",ylim=c(0,0.8))
b1<- rgb(171,214,255,max = 255, alpha = 80)
polygon(c(1:12, rev(1:12)), c(top, rev(bot)),
        col = b1, border = NA)
legend("topright",lty=c(1,NA),pch = c(NA,15), col = c(1,b1),legend = c("Median","95%PI"))
