cycle <- function(N,mu,sd,rp,kd,f,detp,al=42,as=25,t){
  pint <- rtruncnorm(round(N,0),1,al,mu,sd) #initial innoculation from truncated normal distribution
  pint <- round(pint,digits = 0)
  rec <- matrix(0,t,al)
  for (i in 1:al) { #record initial number of parasites
    rec[1,i] <- length(which(pint==i)) 
  }
  for (i in 2:t) { 
    for (j in 1:al) { #asexual parasites
      if (j==1){ #state 1 (age = 1 hour)
        p.a.t <- rec[i-1,al]-rpois(1,(detp+kd[i,3])*rec[i-1,al])
        # make sure it's not negative
        p.a.t <- max(0,round(p.a.t*rp,0))
        rec[i,j] <- p.a.t
      }
        else{if (j == as+1){
        p.a.t <- rec[i-1,as]-rpois(1,(detp+kd[i,3])*rec[i-1,as])
        p.a.t <- max(0,p.a.t)
        rec[i,j] <- p.a.t-rbinom(1,p.a.t,f)
      }
        else{ #states except state 1 and sequestration
          p.a.t <- rec[i-1,j-1]-rpois(1,(detp+kd[i,3])*rec[i-1,j-1])
          p.a.t <- max(0,p.a.t)
          rec[i,j] <- p.a.t
        }
      }
    }
    }
  rec <- data.frame(rec)
  return(rec)
}

