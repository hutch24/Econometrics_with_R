library('quadprog')
data = read.csv('/Users/hutch/OneDrive - 중앙대학교/RinEcometrics/PS2/data_port_PS2.csv', header=T)
data = ts(data,start = c(1990,01),end = c(2018,09), frequency = 12)
head(data)
meanlist = NULL#in sample 자산별 평균을 담는 변수
phi1 = NULL
phi2 = NULL
phi3 = NULL
phi4 = NULL
phi5 = c(rep(1/5,5))
r1 = NULL
r2 = NULL
r3 = NULL
r4 = NULL
r5 = NULL
c = mean(data[,2:6])#제약식(2)의 c는 모든 자산의 global mean으로 설정
data[121,]
#함수 구현
for (i in 1:(dim(data)[1]-120)){#i=1,2,...,225까지 반복
  #in sample, out of sample 정의
  insample = data[(i:(i+120-1)),]#120개행의 in sample정의
  outsample = data[(i+120),]#t+W번째 행인 out of sample정의
  
  #분산공분산행렬 정의
  covmatrix = cov(insample[,2:6])
  #meanlist정의 - 자산별(열별) 평균값 구함
  for (j in 2:dim(insample)[2]){
    meanlist[j-1] = mean(insample[,j])}
  
  #1번 전략(mean-variance 포트폴리오-숏셀가능)
  Q1=solve.QP( Dmat=2*covmatrix, dvec=c(0,0,0,0,0),
               Amat=cbind(meanlist,rep(1,5)), bvec=c(c,1) )
  
  phi1 = cbind(phi1,Q1$solution)
  r1[i] = Q1$solution%*%outsample[2:6]
  
  #2번 전략(minimun variance 포트폴리오)
  Q2=solve.QP( Dmat=2*covmatrix, dvec=c(0,0,0,0,0),
               Amat=cbind(rep(1,5)), bvec=c(1) )
  
  phi2 = cbind(phi2,Q2$solution)
  r2[i] = Q2$solution%*%outsample[2:6]
  
  #3번 전략(mean-variance 포트폴리오-숏셀불가)
  Q3=solve.QP( Dmat=2*covmatrix, dvec=c(0,0,0,0,0),
               Amat=cbind(meanlist,rep(1,5),diag(1,5)), bvec=c(c,1,rep(0,5)) )
  
  phi3 = cbind(phi3, Q3$solution)
  r3[i] = Q3$solution%*%outsample[2:6]
  
  #4번 전략(minimun variance 포트폴리오-숏셀불가)
  Q4= solve.QP( Dmat=2*covmatrix, dvec=c(0,0,0,0,0),
                Amat=cbind(rep(1,5),diag(1,5)), bvec=c(1,rep(0,5)) )
  
  phi4 = cbind(phi4,Q4$solution)
  r4[i] = Q4$solution%*%outsample[2:6]
  
  #5번 전략(균등비율 포트폴리오)
  r5[i] = phi5%*%outsample[2:6] }

r = rbind(r1,r2,r3,r4,r5)

#shape 비율 구하는 함수
sharpe=NULL
for (i in 1:dim(r)[1]) {
  m = mean(r[i,]); std = var(r[i,])^0.5; sharpe[i] = m/std }
sharpe
names(sharpe) = c(1,2,3,4,5)
match(min(sharpe), sharpe) #sharpe비율이 가장 작은 case
match(max(sharpe),sharpe) #sharpe비율이 가장 큰 case

barplot(sharpe, main='Shape ratio', col=c(5,0,0,4,0))

#Q1 120개 이용해서 추정하고, 121번째 값을 이용해서 표본외 수익률 계산하는 거면, 
#전체 데이터가 345개인데, 롤링윈도우가 344번째까지만 실행되어야 하는가?
#Q2 phi'*mu = c 인 2번제약식에서 c값을 무얼로 두냐에 따라 수렴하는 경우와 
#수렴하지 않는 경우가 있는데 c를 어떤 수로 두어야 하는가?


#t(A)%*%b = bvec 이므로
#Amat정의
t(cbind(meanlist,rep(1,5),rep(-1,5),diag(1,5)))
#I(10*10)과 [1]과 [-1]을 세로 결합
#bvec정의
constant=1
c(constant,1, -1,rep(0, 5)) # 0  0  0  0  0  0  0  0  0  0  1 -1

meanlist = c(1,1,1,1,1)
#완성 코드
# case1
2*covmatrix
cbind(meanlist,rep(1,5))
c(constant,1)
Q1=solve.QP(Dmat=2*covmatrix,dvec=c(0,0,0,0,0),Amat=cbind(meanlist,rep(1,5),rep(-1,5)),bvec=c(constant,1,-1))
sol1[i] = Q1$solution%*%x[2:6]
# case2
cbind(rep(1,5))
c(1)
Q2=solve.QP(Dmat=2*covmatrix,dvec=c(0,0,0,0,0),Amat=cbind(rep(1,5),rep(-1,5)),bvec=c(1,-1))
sol2[i] = Q2$solution%*%x[2:6]
# case3
cbind(meanlist,rep(1,5),diag(1,5))
c(constant,1,rep(0,5))
Q3=solve.QP(Dmat=2*covmatrix,dvec=c(0,0,0,0,0),Amat=cbind(meanlist,rep(1,5),rep(-1,5),diag(1,5)),bvec=c(constant,1,-1,rep(0,5)))
sol3[i] = Q3$solution%*%x[2:6]
# case4
cbind(rep(1,5),diag(1,5))
c(1,rep(0,5))
Q4= solve.QP(Dmat=2*covmatrix,dvec=c(0,0,0,0,0),Amat=cbind(rep(1,5),rep(-1,5),diag(1,5)),bvec=c(1,-1,rep(0,5)))
sol4[i] = Q4$solution%*%x[2:6]
# case5
phi = c(rep(1/5,5))
sol5[i] = phi%*%x[2:6]

#백업
QP1 = solve.QP(Dmat=2*covmatrix,dvec=c(0,0,0,0,0),Amat=cbind(matrix(c(1,1,1,1,1),5,1), meanlist), bvec=c(1,1) )
sol1[i] = QP$solution%*%x[2:6]
QP2 = solve.QP(Dmat=2*covmatrix,dvec=c(0,0,0,0,0),Amat=matrix(c(1,1,1,1,1),5,1), bvec=c(1) )
sol2[i] = QP2$solution%*%x[2:6]


mmat  = NULL
stdmat = NULL
for (j in 1:dim(data)[1]){
  mmat[j] = mean(data[j,2:6])
  stdmat[j] = var(data[j,2:6])^0.5}
plot(y=mmat, x=stdmat, xlab='std.err',ylab='mean')

for (k in 1:5){
  lines(phi1[k,],col=k)
}
lines(phi1[,2])
#std.err 구할 수 있나?
phi1[,2]


phi41 = ts(phi4[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi42 = ts(phi4[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi43 = ts(phi4[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi44 = ts(phi4[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi45 = ts(phi4[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

par(mfrow=c(3,2))
plot(phi41, type='l', main='phi1 of Portfolio 4', xlab='', ylab='')
abline(h=0)
plot(phi42, type='l', main='phi2 of Portfolio 4', xlab='', ylab='')
abline(h=0)
plot(phi43, type='l', main='phi3 of Portfolio 4', xlab='', ylab='')
abline(h=0)
plot(phi44, type='l', main='phi4 of Portfolio 4', xlab='', ylab='')
abline(h=0)
plot(phi45, type='l', main='phi5 of Portfolio 4', xlab='', ylab='')
abline(h=0)

mr=NULL
sr=NULL
for (i in 1:5){
  mr[i] = mean(r[i,])
  sr[i] = var(r[i,])^0.5
}
totmat = t(rbind(mr,sr))
plot(totmat)
r1 
plot(density(r1))

plot(density(phi1[1,]))
plot(density(phi1[2,]))
plot(density(phi2))
plot(density(phi3))
plot(density(phi4))
plot(density(phi5))





phi11 = ts(phi1[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi12 = ts(phi1[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi13 = ts(phi1[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi14 = ts(phi1[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi15 = ts(phi1[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

par(mfrow=c(2,3))
plot(phi11, type='l', main='phi of Cnsmr - Portfolio 1', xlab='', ylab='')
abline(h=0)
plot(phi12, type='l', main='phi of Manuf - Portfolio 1', xlab='', ylab='')
abline(h=0)
plot(phi13, type='l', main='phi of HiTec - Portfolio 1', xlab='', ylab='')
abline(h=0)
plot(phi14, type='l', main='phi of Hlth - Portfolio 1', xlab='', ylab='')
abline(h=0)
plot(phi15, type='l', main='phi of Other - Portfolio 1', xlab='', ylab='')
abline(h=0)

phi21 = ts(phi2[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi22 = ts(phi2[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi23 = ts(phi2[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi24 = ts(phi2[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi25 = ts(phi2[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

par(mfrow=c(2,3))
plot(phi21, type='l', main='phi of Cnsmr - Portfolio 2', xlab='', ylab='')
abline(h=0)
plot(phi22, type='l', main='phi of Manuf - Portfolio 2', xlab='', ylab='')
abline(h=0)
plot(phi23, type='l', main='phi of HiTec - Portfolio 2', xlab='', ylab='')
abline(h=0)
plot(phi24, type='l', main='phi of Hlth - Portfolio 2', xlab='', ylab='')
abline(h=0)
plot(phi25, type='l', main='phi of Other - Portfolio 2', xlab='', ylab='')
abline(h=0)


phi31 = ts(phi3[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi32 = ts(phi3[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi33 = ts(phi3[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi34 = ts(phi3[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi35 = ts(phi3[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

par(mfrow=c(2,3))
plot(phi31, type='l', main='phi of Cnsmr - Portfolio 3', xlab='', ylab='')
plot(phi32, type='l', main='phi of Manuf - Portfolio 3', xlab='', ylab='')
plot(phi33, type='l', main='phi of HiTec - Portfolio 3', xlab='', ylab='')
plot(phi34, type='l', main='phi of Hlth - Portfolio 3', xlab='', ylab='')
plot(phi35, type='l', main='phi of Other - Portfolio 3', xlab='', ylab='')

phi41 = ts(phi4[1,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi42 = ts(phi4[2,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi43 = ts(phi4[3,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi44 = ts(phi4[4,],start = c(2000,01),end = c(2018,09), frequency = 12)
phi45 = ts(phi4[5,],start = c(2000,01),end = c(2018,09), frequency = 12)

par(mfrow=c(2,3))
plot(phi41, type='l', main='phi of Cnsmr - Portfolio 4', xlab='', ylab='')
plot(phi42, type='l', main='phi of Manuf - Portfolio 4', xlab='', ylab='')
plot(phi43, type='l', main='phi of HiTec - Portfolio 4', xlab='', ylab='')
plot(phi44, type='l', main='phi of Hlth - Portfolio 4', xlab='', ylab='')
plot(phi45, type='l', main='phi of Other - Portfolio 4', xlab='', ylab='')
