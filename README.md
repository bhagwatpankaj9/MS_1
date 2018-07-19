# MS_1

n=100
m=100
alpha1 <- 1
alpha2 <- 1
delta1 <- 1
delta2 <- 1


Est_p_gamma_modified <- function(l,data){
  p_i <- rbeta(n,alpha1,alpha2)
  gama_j <- rgamma(m,delta1,scale = delta2)
  r = max(data)
  N_k <- n_k(data)
  F <- matrix(0,nrow = n,ncol = m)
  for(i in 1:n){
    for(j in 1:m){
      F[i,j] <-  f1(r,N,N_k,p_i[i],gama_j[j])
    }
  }
  
  K_l <- sum(F)
  rho_l <- sum(F%*%gama_j**l)
  theta_l <- sum(t(F)%*%p_i**l)
  
  return(list("alpha1" =  alpha1,"alpha2" = alpha2,"delta1"=delta1,"delta2" = delta2,"gamma_l" = rho_l/K_l, "p_l" = theta_l/K_l))
}

f1 <-function(r,N,N_k,p,gama){
  f <- 1
  for(t in 1:r){
    f <- f*((1-p)**(t*N_k[t]))*((gama+t-1)**(sum(N_k[t:r])))
  }
  return(f*p**(N*gama))
}

n_k <- function(d){
  count_k <- list()
  for(i in 1:max(d)){
    count_k[i] <- length(which(d==i))
  }
  return(as.numeric(count_k))
}


#Simulation
gama <- 1 
p <- 0.5
data <-rnbinom(N,gama,p)
Sums <- c()
for(s in 0:max(data)){
  Sums <- c(Sums,length(which(data==s)))
}

plot(hist(data))

Est_p_gamma_modified(1,data)

#######################
### Replicates
set.seed(123)
N <- length(data)
results_p <- c()
results_gamma <- c()
replics <- 10
for(replic in 1:replics){
  result <- Est_p_gamma_modified(1,data)
  results_p <- c(results_p,result$p_l)
  results_gamma <- c(results_gamma,result$gamma_l)
}

results_p <- apply(matrix(results_p,nrow = 1), 1, function(x) x[!is.nan(x)] )
results_gamma <- apply(matrix(results_gamma,nrow = 1), 1, function(x) x[!is.nan(x)] )
plot(density(results_p))
p_est <- mean(results_p)
plot(density(results_gamma))
gama_est <- mean(results_gamma)

Sums_est <- c()
for(x in 0:max(data)){
  Sums_est <- c(Sums_est,N*dnbinom(x,gama_est,p_est))
}
Sums_est
sum((Sums - Sums_est)**2)

#Comparison With MethodOfMomentsEstimators
data <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
          0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
          2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 
          3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 
          5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 8, 8, 8)

Sums <- c()
for(s in 0:max(data)){
  Sums <- c(Sums,length(which(data==s)))
}

Sums

Sums_est_JM <- c()
for(x in 0:14){
  Sums_est_JM <- c(Sums_est_JM,N*dnbinom(x,a,theta))  ###### a, theta are from MethodOfMomentsin ShuntrAccidents
}
Sums_est_JM
sum((Sums - Sums_est_JM)**2)



