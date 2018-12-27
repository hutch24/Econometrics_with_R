
spot = read.csv('/Users/hutch/OneDrive - 중앙대학교/RinEcometrics/PS1/spot_future.csv', header=T)
spot=ts(spot, start = c(1990,3,23), end = c(2010,3,24), frequency = 365.25)
head(diff(spot))


lnspot = log(spot[,-1]) #시간자료인 date변수만 빼고 로그취해줌
diffspot = diff(lnspot, lag=1, differences = 1) #변화율을 구하기 위해 로그화된 각 변수를 차분


head(diffspot)
dim(diffspot)

par(mfrow=c(2,2), mex=0.7)

#(5) 최적헷지비율 구하는 함수
coeffmatrix = NULL 
resimatrix = NULL
for (i in 1:(dim(diffspot)[2]/2) ) { 
  model = lm( diffspot[,(2*i-1)] ~ diffspot[,(2*i)] ) #홀수번째 변수인 ln(현물) 에 대해 
  #짝수번째 변수인 ln(선물)을 적합시킨다 
  coeff = summary(model)$coefficients[2,] 
  coeffmatrix = rbind(coeffmatrix, coeff) #최적헷지비율 추정
  resimatrix = cbind(resimatrix, summary(model)$residuals) #모델별 잔차 matrix 
}

coeffmatrix #최적헷지비율 추정

collist =c('soybean','soyoil', 'wheat', 'corn', 'coffee' )
#최적헷지비율의 추정치 matrix
coef.LSE = coeffmatrix[,1]
names(coef.LSE) = collist

#최적헷지비율의 Standard Error matrix
std.coef.LSE = coeffmatrix[,2]
names(std.coef.LSE) = collist

coef.LSE #모델별 최적헷지비율 추정치
std.coef.LSE #최적헷지비율의 Standard Error

plot(coef.LSE, type='l', main="Optimal Hedge Ratio(LSE)", 
     xlab='',ylab='',ylim=c(0.5, 1), axes=FALSE)
axis(side = 1, at = 1:5, labels = collist[1:5], line = 0)
axis(side = 2, lty = 1, line=1)
lines(coef.LSE+1.96*as.numeric(coeffmatrix[,2]), col="blue", lty="dashed")
lines(coef.LSE-1.96*as.numeric(coeffmatrix[,2]), col="blue", lty="dashed")

par(mfrow=c(3,2))
#(6) 르장드레 beta(t)
t = (1:(dim(B_t_1)[1])) #t=1:7305
t.norm <- (t-min(t))/(max(t)-min(t)) #정규화된 t
t.norm.in <- as.matrix(t.norm)

#자산별 information행렬 구축
P_selector = NULL
for (k in 1:(dim(diffspot)[2]/2)) {#총 5개 자산
  ft= diffspot[1:7305,(2*k)] #독립변수 ft
  st = matrix(diffspot[,((2*k)-1)]) #종속변수 st
  #information행렬 정의
  info <- matrix(NA, nrow=2, ncol=6)
  for (j in 1:6){
    beta <- indep(t.norm.in, j) # beta(t)
    model_6 <- lm(st ~cbind(ft*beta))
    info[1,j] <- extractAIC(model_6)[2]  				                	#AIC
    info[2,j] <- extractAIC(model_6, 
                            k=log(length(diffspot[,((2*k)-1)])))[2]#SIC
  }
  info <- rbind(c(1:6),info)
  AICn <- info[1,order(info[2,])[1]]   #AIC
  SICn <- info[1,order(info[3,])[1]]   #SIC
  P_selector = rbind(P_selector, cbind(AICn, SICn)) #자산별 P차수 행렬

  P_selector
  
  nn = P_selector[k,2] #SIC
  MM <- indep(t.norm.in,nn)
  b <- as.matrix(coef(model_6))
  cov.c <- vcov(model_6, 'robust' = TRUE)
  cov.1 <- cov.c[2:(2+nn), 2:(2+nn)]
  sd.1 <- rep(0, length(t.norm.in))
  #cov matrix정의
   for (i in 1:length(t.norm.in)){ 
     ttrd <- rbind(1,legendre(1,t.norm.in[i]),legendre(2,t.norm.in[i]),
                    legendre(3,t.norm.in[i]),legendre(4,t.norm.in[i]),legendre(5,t.norm.in[i]),
                    legendre(6,t.norm.in[i]),legendre(7,t.norm.in[i]), legendre(8,t.norm.in[i]),
                    legendre(9,t.norm.in[i]),legendre(10,t.norm.in[i]))
     ttrd <- as.matrix(ttrd[(1:(nn+1)),1])
     sd.1[i] <- sqrt(t(ttrd)%*%cov.1%*%ttrd) }
  
  Le.M <- cbind(1,legendre(1,t.norm.in),legendre(2,t.norm.in),legendre(3,t.norm.in),
              legendre(4,t.norm.in),legendre(5,t.norm.in),legendre(6,t.norm.in),
              legendre(7,t.norm.in),legendre(8,t.norm.in),legendre(9,t.norm.in),
              legendre(10,t.norm.in))
  Le.M <- Le.M[,1:(1+nn)]
  bb.6 <- Le.M%*%b[2:(nn+2),1]
  bb.11 <- bb.6[7305] #out-of-sample betas
  #신뢰구간 구축
  bb.1.lo <- bb.6 - 1.96*sd.1 
  bb.1.up <- bb.6 + 1.96*sd.1
  
  
  label = c('β(t) of','Asset',k)
  plot(t.norm.in,bb.6, main=label)
  
  seq <- order(t.norm.in)
lines(t.norm.in[seq],bb.6[seq])
lines(t.norm.in[seq],bb.1.lo[seq],col="red")
lines(t.norm.in[seq],bb.1.up[seq],col="red")
rug(jitter(t.norm.in),col="blue")

}




par(mfrow=c(3,2))
#(7) beta( B_t-1 )
#자산별 basis 정의
basis_bean = matrix(lnspot[,1]- lnspot[,2])
basis_soyoil = matrix(lnspot[,3]- lnspot[,4])
basis_wheat = matrix(lnspot[,5]- lnspot[,6])
basis_corn = matrix(lnspot[,7]- lnspot[,8])
basis_coffee = matrix(lnspot[,9]- lnspot[,10])
#각 자산별 basis행렬=> B_t-1 (7305x5)
B_t_1 = cbind(basis_bean, basis_soyoil, basis_wheat, basis_corn, 
              basis_coffee)[1:7305,] 

#정규화된 basis행렬 만들기
norm_basis_matrix = NULL
for (i in 1:dim(B_t_1)[2]){ #자산 갯수 i=1:5
  norm_basis <- (B_t_1[,i] - min(B_t_1[,i]))/(max(B_t_1[,i])-min(B_t_1[,i]))
  norm_basis_matrix = cbind(norm_basis_matrix, norm_basis) #정규화된 basis행렬
}
#자산별 정규화된 basis(7305x1)
norm_basis_matrix = matrix(norm_basis_matrix[,i])

#자산별 information행렬
P_selector = NULL
for (k in 1:(dim(diffspot)[2]/2) ) { 
  ft= diffspot[1:7305,(2*k)]
  st = matrix(diffspot[,((2*k)-1)])
  #information행렬 정의
  info <- matrix(NA, nrow=2, ncol=6)
  for (j in 1:6){
    beta = indep(norm_basis_matrix, j) # beta(B_t-1)
    model_6 = lm(st ~cbind(ft*beta))
    info[1,j] <- extractAIC(model_6)[2]  				                	# AIC
    info[2,j] <- extractAIC(model_6, k=log(length(diffspot[,((2*k)-1)])))[2]#SIC
  }
  info <- rbind(c(1:6),info)
  AICn <- info[1,order(info[2,])[1]]   #AIC
  SICn <- info[1,order(info[3,])[1]]   #SIC
  P_selector = rbind(P_selector, cbind(AICn, SICn))

P_selector
nn = P_selector[k,1] #AIC
MM <- indep(norm_basis_matrix,nn)
b <- as.matrix(coef(model_6))
cov.c <- vcov(model_6, 'robust' = TRUE)
cov.1 <- cov.c[2:(2+nn), 2:(2+nn)]
sd.1 <- rep(0, length(norm_basis_matrix))
#cov matrix구현
for (i in 1:length(norm_basis_matrix)){ 
  ttrd <- rbind(1,legendre(1,norm_basis_matrix[i]),legendre(2,norm_basis_matrix[i]),
                legendre(3,norm_basis_matrix[i]),legendre(4,norm_basis_matrix[i]),
                legendre(5,norm_basis_matrix[i]),
                legendre(6,norm_basis_matrix[i]),legendre(7,norm_basis_matrix[i]), 
                legendre(8,norm_basis_matrix[i]),legendre(9,norm_basis_matrix[i]),
                legendre(10,norm_basis_matrix[i]))
  ttrd <- as.matrix(ttrd[(1:(nn+1)),1])
  sd.1[i] <- sqrt(t(ttrd)%*%cov.1%*%ttrd) }

Le.M <- cbind(1,legendre(1,norm_basis_matrix),legendre(2,norm_basis_matrix),legendre(3,norm_basis_matrix),
              legendre(4,norm_basis_matrix),legendre(5,norm_basis_matrix),
              legendre(6,norm_basis_matrix),legendre(7,norm_basis_matrix),legendre(8,norm_basis_matrix),
              legendre(9,norm_basis_matrix),legendre(10,norm_basis_matrix))
Le.M <- Le.M[,1:(1+nn)]
bb.1 <- Le.M%*%b[2:(nn+2),1]
bb.11 <- bb.1[7305] #out-of-sample betas
#신뢰구간 지정
bb.1.lo <- bb.1 - 1.96*sd.1 
bb.1.up <- bb.1 + 1.96*sd.1

label = c('β(B_t-1) of','Asset',k)
plot(norm_basis_matrix,bb.1, main=label)
seq <- order(norm_basis_matrix)
lines(norm_basis_matrix[seq],bb.1[seq])
lines(norm_basis_matrix[seq],bb.1.lo[seq],col="red")
lines(norm_basis_matrix[seq],bb.1.up[seq],col="red")
rug(jitter(norm_basis_matrix),col="blue")

}
plot(norm_basis_matrix[,5],bb.1, main=label)

var(bb.1)
#legendre() 정의
legendre = function(n,x){
  x = 2*x - 1
  if (n==0) lp <- 1
  if (n==1) lp <- x
  if (n==2) lp <- 0.5*(3*x^2 - 1)
  if (n==3) lp <- 0.5*(5*(x^3) - 3*x)
  if (n==4) lp <- (1/8)*(35*x^4 - 30*x^2 + 3)
  if (n==5) lp <- (1/8)*(63*x^5 - 70*x^3 + 15*x)
  if (n==6) lp <- (1/16)*(231*x^6 - 315*x^4 + 105*x^2 - 5)
  if (n==7) lp <- (1/16)*(429*x^7 - 693*x^5 + 315*x^3 -35*x)
  if (n==8) lp <- (1/128)*(6435*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35)
  if (n==9) lp <- (1/128)*(12155*x^9 - 25730*x^7 + 18018*x^5 - 4620*x^3 + 315*x)
  if (n==10) lp <- (1/256)*(46189*x^(10) - 109395*x^8 + 90090*x^6 - 30030*x^4 + 3465*x^2 - 63)
  return(lp)
}
# Time-varying regressor
indep <- function(x,n){
  MM <- matrix(1, nrow=length(x),ncol=(n+1))
  for (i in 1:n){
    MM[,i+1] <- legendre(i,x)
  }
  return(MM)
}



#(6)+(8)
ft = diffspot[1:7305,2] # Δft
beta <- indep(t.norm.in, 3) #SIC = 3
soybeans = matrix(diffspot[,1]) #종속변수 Δst

soybean.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                                submodel = NULL, external.regressors = NULL, 
                                                variance.targeting = FALSE ),
                          mean.model = list(armaOrder=c(3,1), include.mean = TRUE, 
                                            external.regressors = cbind(ft*beta)),
                          distribution.model = "norm", start.pars = list(), fixed.pars=list())
soybean.fit = ugarchfit(spec=soybean.spec, data=soybeans, solver.control = list(tol = 1e-5))
soybean.fit
coef(soybean.fit) #9.236609e-01
std.soybean.fit = 0.004197


ft = diffspot[1:7305,4]
beta <- indep(t.norm.in, 5)
soyoilf  = cbind(ft*beta)
soyoils = matrix(diffspot[,3]) #종속변수 st

soyoil.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                               submodel = NULL, external.regressors = NULL, variance.targeting = FALSE ),
                         mean.model = list(armaOrder=c(4,1), include.mean = TRUE, 
                                           external.regressors = soyoilf ), 
                         distribution.model = "norm", start.pars = list(), fixed.pars=list())
soyoil.fit = ugarchfit(spec=soyoil.spec, data=soyoils, solver.control = list(tol = 1e-4))
soyoil.fit
coef(soyoil.fit) # 0.849096 
std.soyoil.fit = 0.004895


ft = diffspot[1:7305,6]
beta <- indep(t.norm.in, 5)
wheatf  = cbind(ft*beta)
wheats = matrix(diffspot[,5]) #종속변수 st

wheat.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                              submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                        mean.model = list(armaOrder=c(2,1), include.mean = TRUE, 
                                          external.regressors = wheatf), 
                        distribution.model = "norm", start.pars = list(), fixed.pars=list())
wheat.fit = ugarchfit(spec=wheat.spec, data=wheats, solver.control = list(tol = 1e-12))
wheat.fit
coef(wheat.fit) #0.795488
std.wheat.fit =   0.006146



ft = diffspot[1:7305,8]
beta <- indep(t.norm.in, 6)
cornf  = cbind(ft*beta)
corns = matrix(diffspot[,7]) #종속변수 st

corn.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                             submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                       mean.model = list(armaOrder=c(0,1), include.mean = TRUE, 
                                         external.regressors = cornf), 
                       distribution.model = "norm", start.pars = list(), fixed.pars=list())
corn.fit = ugarchfit(spec=corn.spec, data=corns, solver.control = list(tol = 1e-12))
corn.fit
coef(corn.fit) #0.904010
std.corn.fit = 0.006055


ft = diffspot[1:7305,10]
beta <- indep(t.norm.in, 6)
coffeef  = cbind(ft*beta)
coffees = matrix(diffspot[,9]) #종속변수 st

coffee.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                               submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                         mean.model = list(armaOrder=c(1,1), include.mean = TRUE, 
                                           external.regressors = coffeef), 
                         distribution.model = "norm", start.pars = list(), fixed.pars=list())
coffee.fit = ugarchfit(spec=coffee.spec, data=coffees, solver.control = list(tol = 1e-12))
coffee.fit
coef(coffee.fit) #0.664715
std.coffee.fit = 0.008801
coef(soybean.fit)[6]
coef(soyoil.fit)[7]
coef(wheat.fit)[5]
coef(corn.fit)[3]
coef(coffee.fit)[4]
coef.Legendre_GARCH = c(coef(soybean.fit)[6],coef(soyoil.fit)[7],coef(wheat.fit)[5],coef(corn.fit)[3],coef(coffee.fit)[4])
std.GARCH= c(std.soybean.fit,std.soyoil.fit ,std.wheat.fit,std.corn.fit,std.coffee.fit)
names(coef.Legendre_GARCH) = collist
names(std.GARCH) = collist

coef.Legendre_GARCH #최적헷지비율의 추정치
std.GARCH #최적헷지비율의 std.err의 추정치

par(mfrow=c(1,1))
plot(coef.Legendre_GARCH, type='l', main="Optimal Hedge Ratio(Legendre, GARCH)", xlab='',ylab='',
     ylim=c(0.7,1),axes=FALSE)
lines(coef.Legendre_GARCH+1.96*std.GARCH, col="blue", lty="dashed")
lines(coef.Legendre_GARCH-1.96*std.GARCH, col="blue", lty="dashed")
axis(side = 1, at = 1:5, labels = collist[1:5], line = 0)
axis(side = 2, lty = 1, line=1)


#(7)+(8)
#자산별 정규화된 basis
norm_basis_matrix = NULL
for (i in 1:dim(B_t_1)[2]){ #자산 갯수 i=1:5
  norm_basis <- (B_t_1[,i] - min(B_t_1[,i]))/(max(B_t_1[,i])-min(B_t_1[,i]))
  norm_basis_matrix = cbind(norm_basis_matrix, norm_basis) #정규화된 basis행렬
}

#AIC = 5,5,6,6,6
ft = diffspot[1:7305,2] # Δft
beta <- indep(matrix(norm_basis_matrix[,1]), 5) #AIC =5
soybeans = matrix(diffspot[,1]) #종속변수 Δst
soybean.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                                submodel = NULL, external.regressors = NULL, 
                                                variance.targeting = FALSE ),
                          mean.model = list(armaOrder=c(3,1), include.mean = TRUE, 
                                            external.regressors = cbind(ft*beta)),
                          distribution.model = "norm", start.pars = list(), fixed.pars=list())
soybean.fit = ugarchfit(spec=soybean.spec, data=soybeans, solver.control = list(tol = 1e-5))
soybean.fit
coef(soybean.fit) #1.050150
std.soybean.fit = 0.000338


ft = diffspot[1:7305,4]
beta <- indep(matrix(norm_basis_matrix[,2]), 5)
soyoilf  = cbind(ft*beta)
soyoils = matrix(diffspot[,3]) #종속변수 st
soyoil.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                               submodel = NULL, external.regressors = NULL, variance.targeting = FALSE ),
                         mean.model = list(armaOrder=c(4,1), include.mean = TRUE, 
                                           external.regressors = soyoilf ), 
                         distribution.model = "norm", start.pars = list(), fixed.pars=list())
soyoil.fit = ugarchfit(spec=soyoil.spec, data=soyoils, solver.control = list(tol = 1e-4))
soyoil.fit
coef(soyoil.fit) #0.980299
std.soyoil.fit = 0.016652


ft = diffspot[1:7305,6]
beta <- indep(matrix(norm_basis_matrix[,3]), 5)
wheatf  = cbind(ft*beta)
wheats = matrix(diffspot[,5]) #종속변수 st
wheat.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                              submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                        mean.model = list(armaOrder=c(2,1), include.mean = TRUE, 
                                          external.regressors = wheatf), 
                        distribution.model = "norm", start.pars = list(), fixed.pars=list())
wheat.fit = ugarchfit(spec=wheat.spec, data=wheats, solver.control = list(tol = 1e-12))
wheat.fit
coef(wheat.fit) #1.102288
std.wheat.fit =   0.037372



ft = diffspot[1:7305,8]
beta <- indep(matrix(norm_basis_matrix[,4]), 6)
cornf  = cbind(ft*beta)
corns = matrix(diffspot[,7]) #종속변수 st
corn.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                             submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                       mean.model = list(armaOrder=c(0,1), include.mean = TRUE, 
                                         external.regressors = cornf), 
                       distribution.model = "norm", start.pars = list(), fixed.pars=list())
corn.fit = ugarchfit(spec=corn.spec, data=corns, solver.control = list(tol = 1e-12))
corn.fit
coef(corn.fit) #1.075283
std.corn.fit = 0.026598


ft = diffspot[1:7305,10]
beta <- indep(matrix(norm_basis_matrix[,5]), 6)
coffeef  = cbind(ft*beta)
coffees = matrix(diffspot[,9]) #종속변수 st
coffee.spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
                                               submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                         mean.model = list(armaOrder=c(1,1), include.mean = TRUE, 
                                           external.regressors = coffeef), 
                         distribution.model = "norm", start.pars = list(), fixed.pars=list())
coffee.fit = ugarchfit(spec=coffee.spec, data=coffees, solver.control = list(tol = 1e-4))
coffee.fit
coef(coffee.fit) #0.861147
std.coffee.fit = 0.021387


coef.LegendreB_GARCH = c(coef(soybean.fit)[6],coef(soyoil.fit)[7],coef(wheat.fit)[5],coef(corn.fit)[3],coef(coffee.fit)[4])
std.GARCH= c(std.soybean.fit,std.soyoil.fit ,std.wheat.fit,std.corn.fit,std.coffee.fit)
names(coef.LegendreB_GARCH) = collist
names(std.GARCH) = collist

coef.GARCH #최적헷지비율 추정치
std.GARCH #최적헷지비율 std.err의 추정치

par(mfrow=c(1,1))
plot(coef.LegendreB_GARCH, type='l', main="Optimal Hedge Ratio(Legendre, GARCH)", xlab='',ylab='',
     ylim=c(0.75,1.2),axes=FALSE)
lines(coef.LegendreB_GARCH+1.96*std.GARCH, col="blue", lty="dashed")
lines(coef.LegendreB_GARCH-1.96*std.GARCH, col="blue", lty="dashed")
axis(side = 1, at = 1:5, labels = collist[1:5], line = 0)
axis(side = 2, lty = 1, line=1)



#각 변수들의 분산을 구하는 함수
vars = NULL
for ( i in 1:dim(diffspot)[2] )  {
  vars[i] = var(diffspot[,i])
}
# Cov(st, ft)를 구하는 함수
covs = NULL
for ( j in 1:dim(diffspot)[2]/2 ) {
  covs[j] = cov( diffspot[,(2*j-1)], diffspot[,(2*j)] )
}

names(vars) = c("soybeans","soybeanf","soyaoils","soyaoilf","wheats","wheatf",
                "corns","cornf","coffees","coffeef")
names(covs) = collist
vars
covs

#Var(Pt) matrix구하기
coefs = rbind(coef.LSE, #coef.Legendre자리 ,coef.LegendreB자리,
              coef.LSE_GARCH,coef.Legendre_GARCH,coef.LegendreB_GARCH)
varPt=NULL
varPtmatrix=NULL
for (j in 1:dim(coefs)[1]) { #j = 1:6
  for (i in 1:dim(diffspot)[2]/2) { #총 5개 자산
    varPt[i]= vars[(2*i-1)] + (coefs[j,i]^2)*vars[(2*i)]-2*(coefs[j,i])*covs[i]
  } #Var(Δst)+β^2×Var(Δft)-2×β×Cov(Δst,Δft)
  varPtmatrix = rbind(varPtmatrix,varPt) #6개 모델, 3개 자산에 대해 분산산출
}

names(varPt) = collist
varPt
#르장드레 끼워넣기
coef.Legendre = NULL
coef.LegendreB = NULL
for (i in 1:5){
  coef.Legendre = cbind(coef.Legendre,var(diffspot[,((2*i)-1)]) + var(bb.6*diffspot[,2*i])-2*cov(diffspot[,((2*i)-1)], bb.6*diffspot[,2]))
  coef.LegendreB = cbind(coef.LegendreB,var(diffspot[,((2*i)-1)]) + var(bb.1*diffspot[,2*i])-2*cov(diffspot[,((2*i)-1)], bb.1*diffspot[,2]))
}
#합치기
varPtmatrix = rbind(varPtmatrix[1,], coef.Legendre,coef.LegendreB,varPtmatrix[2:4,])
colnames(varPtmatrix) = collist
rownames(varPtmatrix) = c("LSE","Legendre β(t)","Legendre β(B_t-1)",
                          "GARCH",
                          "Legendre β(t) with GARCH","Legendre β(B_t-1) with GARCH")

varPtmatrix

 
plot(varPtmatrix[1,], type='l', xlab='',ylab='',xaxt="n",
     ylim=c(min(varPtmatrix),max(varPtmatrix)), main="Var(Pt)")
lines(varPtmatrix[2,], type='l', col='red', lty=2)
lines(varPtmatrix[3,], type='l', col='blue')
lines(varPtmatrix[4,], type='l', col='green')
lines(varPtmatrix[5,], type='l', col='orange')
lines(varPtmatrix[6,], type='l', col='pink')
legend(x = 'topleft', legend=c("LSE","Legendre β(t)","Legendre β(B_t-1)","GARCH",
                               "Legendre β(t) with GARCH","Legendre β(B_t-1) with GARCH"),
       col=c("black","red", "blue","green","orange","pink"), lty=c(1,2,1,1,1,1), bg = "white", cex = 1)
axis(side = 1, at = 1:5, labels =collist, line = 0, cex=5)
