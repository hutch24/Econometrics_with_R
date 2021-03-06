library('quadprog')
data1 = read.csv('/Users/hutch/OneDrive - 중앙대학교/RinEcometrics/Econ_R_Final_2018/data_port_FINAL.csv', header = T)

samplemean = NULL
return_minvar = NULL
return_resample = NULL
return_naive = NULL
phi_minvar = NULL
phi_b_matrix = NULL
phi_naive = c(rep(1/5,5))


#B = 2000 #resample 횟수
Blist = c(100,300,500,700,900,1100,1300,1500,2000,2500,3000,3500,4000)
s_mat=NULL
sharpe=NULL
#for (B in Blist){

for (i in 1:(dim(data1)[1]-120)){ #i = 1,2,...,225
  insample = data1[(i:(i+120-1)),]#(120x6)크기 in sample정의
  outsample = data1[(i+120),]#(1x6)크기 out of sample정의
  
  #표본 평균matrix
  for (j in 2:dim(insample)[2]){
    samplemean[j-1] = mean(insample[,j])
  }
  samplemean = as.vector(samplemean)
  #분산공분산matrix
  samplecov = cov(insample[,2:6]) 

  #최소분산 포트폴리오
  Q_minvar = solve.QP(Dmat=2*samplecov, dvec=c(0,0,0,0,0),
                      Amat=cbind(rep(1,5),diag(1,5)), bvec=c(1,rep(0,5)) )
  #수익률 변수 return_minvar
  return_minvar[i] = Q_minvar$solution %*% t(outsample[2:6])
  #phi를 담는 행렬 정의
  phi_minvar = cbind(phi_minvar, Q_minvar$solution)
  
  #Resampled 포트폴리오
  phi_sum = c(0,0,0,0,0) #phi_b들의 합
  phi_RS = c(0,0,0,0,0) #phi_RS=(1/B)*phi_sum
  for (k in B){       #B번 무작위 추출
    resample = mvrnorm(n=120, mu=samplemean, Sigma=samplecov) #(120x5) Resampled matrix
    Q_resample = solve.QP(Dmat = 2*cov(resample), dvec=c(0,0,0,0,0),
                          Amat=cbind(rep(1,5),diag(1,5)), bvec=c(1,rep(0,5)) )
    phi_b = Q_resample$solution
    phi_sum = (phi_sum + phi_b) #phi_b들의 합
    phi_b_matrix = cbind(phi_b_matrix, phi_b) #그래프를 그리기 위한phi_b 모음 행렬
  }
  phi_RS = phi_sum / B #phi_RS=(1/B)*phi_sum
  return_resample[i] = phi_RS %*% t(outsample[2:6]) #수익률
  
  #Naive 포르트폴리오
  #phi_naive = c( rep(1/5,5) )
  return_naive[i] = phi_naive %*% t(outsample[2:6])
}


#shape 비율 구하는 함수
r = rbind(return_minvar, return_resample, return_naive)
#r = return_resample
r
for (i in 1:dim(r)[1]) {
  m = mean(r[i,]); std = var(r[i,])^0.5; 
  sharpe[i] = m/std }
sharpe

names(sharpe) = c('minvar', 'resample', 'naive')
match(min(sharpe), sharpe) #sharpe비율이 가장 작은 case => naive
match(max(sharpe),sharpe) #sharpe비율이 가장 큰 case => resample or min_var


sharpe_RS = s_mat[1,]
names(sharpe_RS) = Blist

par(mfrow=c(1,1))
plot(sharpe_RS,xlab='', type='l',main='Sharpe by B times', xaxt="n", ylim=c(0.174, max(s_mat)))
axis(side = 1, at = 1:13, labels = Blist, line = 0, cex=5)
abline(h=0.193261, col='blue')
abline(h=0.1743525, col='green')
legend("bottomright",  c("Min variance", "Resample", "Naive"), 
       lty = 1, col = c("black", "blue" ,"green"))

barplot(sharpe, main='Shape ratio', col=c(4,0,0))


#수익률 플롯
par(mfrow=c(3,1))
plot(return_minvar, type='l')
plot(return_resample, type='l')
plot(return_naive, type='l')




#phi of portfolio - minimum vaariance
phi11 = ts(phi_minvar[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi12 = ts(phi_minvar[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi13 = ts(phi_minvar[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi14 = ts(phi_minvar[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi15 = ts(phi_minvar[5,],start = c(2000,01),end = c(2018,09), frequency = 12)


#phi of total assets - min-variance
par(mfrow=c(1,1))
plot(phi11, type='l', main='phi Portfolio 1', xlab='', ylab='',
     ylim=c(min(phi_minvar),max(phi_minvar)))
lines(phi12, type='l',col='blue')
lines(phi13, type='l',col='yellow')
lines(phi14, type='l',col='green')
lines(phi15, type='l',col='red')
legend("center",  c("Cnsmr", "Manuf", "HiTec", "Hlth", "Other"), 
       lty = 1, col = c("black", "blue" ,"yellow", "green", "red"))


#phi of portfolio - resample
phi21 = ts(phi_b_matrix[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi22 = ts(phi_b_matrix[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi23 = ts(phi_b_matrix[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi24 = ts(phi_b_matrix[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi25 = ts(phi_b_matrix[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

#phi of total assets - resample
par(mfrow=c(1,1))
plot(phi21, type='l', main='phi Portfolio 2', xlab='', ylab='', 
     ylim=c(min(phi_b_matrix),max(phi_b_matrix)))
lines(phi22, type='l',col='blue')
lines(phi23, type='l',col='yellow')
lines(phi24, type='l',col='green')
lines(phi25, type='l',col='red')
legend("center",  c("Cnsmr", "Manuf", "HiTec", "Hlth", "Other"), 
      lty = 1, col = c("black", "blue" ,"yellow", "green", "red"))



#phi of portfolio - Naive
phi31 = ts(phi_b_matrix[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi32 = ts(phi_b_matrix[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi33 = ts(phi_b_matrix[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi34 = ts(phi_b_matrix[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi35 = ts(phi_b_matrix[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

#phi of total assets - Naive
par(mfrow=c(1,1))
plot(phi31, type='l', main='phi Portfolio 1', xlab='', ylab='',
     ylim=c(min(phi_minvar),max(phi_minvar)))
lines(phi32, type='l',col='blue')
lines(phi33, type='l',col='yellow')
lines(phi34, type='l',col='green')
lines(phi35, type='l',col='red')
legend("center",  c("Cnsmr", "Manuf", "HiTec", "Hlth", "Other"), 
       lty = 1, col = c("black", "blue" ,"yellow", "green", "red"))
