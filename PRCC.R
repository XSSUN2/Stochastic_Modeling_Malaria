#vol102, vol104 and vol304 are the results from file "volxxx.csv"
#An example of vol102, changing all "102" to "104" or "304" can get the corresponding PRCC result
p <- vector()
for (i in 1:1000) {
  p[i] <- sum(vol102[i,]<Inf)/100
}

dataset <- as.data.frame(cbind(p,data102[,1:9]))


rank_trans_data <- matrix(0,1000,10)
for (i in 1:10){
  rank_trans_data[,i] <- rank(dataset[,i])
}
rank_trans_data <- as.data.frame(rank_trans_data)
names(rank_trans_data) <- c("P","Inoculation size"
                            ,'mu',"sd","Replication number", "k_max","EC50","gamma","f","delta_p")

prcc_result <- vector()
for (i in 1:9){
  lm_p <- lm(P~. , data = rank_trans_data[,-(i+1)])
  lm_para <- lm(rank_trans_data[,(i+1)]~. ,data = rank_trans_data[,-c(1,i+1)])
  prcc_result[i] <- cor(lm_p$residuals,lm_para$residuals)
}

#the order used are from largest to smallest in vol102
col_names <- c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]), expression(k[max]),expression(EC[50]),expression(gamma),"f",expression(delta[p]))
order102 <- c(6,2,4,9,3,8,1,7,5)
prcc_result_order <- prcc_result[order102, drop=TRUE]
col_names_order <- col_names[order102, drop=TRUE]

dotchart(prcc_result_order, labels=col_names_order,
         xlab="Volunteer 102", pch=19, xlim=c(-0.8,0.8),
         main="PRCC between parameters and extinction probability",cex=1.25)
abline(v=0,lty=2)

