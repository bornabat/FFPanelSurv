

install.packages('frailtySurv') #for this I had to change my CRAN mirror. Try  chooseCRANmirror(ind=68) if it did not work.
install.packages('frailtyHL')
install.packages('frailtyEM')

library(FFPanelSurv) #our model
library(frailtySurv)
library(frailtyHL)
library(frailtyEM)

data = read.csv('ffdatatest.csv') #Whatever data you generate. It has to have an "id" column, a "response" column for the outcomes, a "censor" column that equals 1 if outcome is censored, and some other features (X).

#or

data = FFSurv_gen(50,4, baseline_hazard = list(formula = "y = log(1+x/125)", max = 300))

data = data[, !colnames(data) %in% 'intercept']

Model1 = FFPanelSurv::FFSurv_est(data)

Model2 = frailtySurv::fitfrail(Surv(response,I(1-censor)) ~ x1+ x2 + x3 + cluster(id), data, frailty = 'gamma')

Model3 = frailtyHL::frailtyHL(Surv(response,I(1-censor)) ~ x1+ x2 + x3 + (1|id), data, RandDist = "Gamma")

Model4 = frailtyEM::emfrail(Surv(response,I(1-censor)) ~ x1+ x2 + x3 + cluster(id), data)

#an array of coefficients

beta_table = matrix(nrow = 4, ncol = 3)

rownames(beta_table) = c("FFpanelsurv", "Coxph shared", "H-likelihood", "EM")

beta_table[1,] = as.numeric(Model1$Coefficients[2:4,1])
beta_table[2,] = Model2$beta
beta_table[3,] = Model3$FixCoef[,1]
beta_table[4,] = Model4$coefficients

#estimated baseline hazard

delta1 = Model1$`Baseline Hazard`

delta1 = cbind(c(1:max(data$response)), delta1)

delta2 = Model2$Lambda

delta2 = delta2[2:nrow(delta2),]

delta4 = cbind(delta2[,1],cumsum(Model4$hazard))

delta_oracle = cbind(delta1[,1], cumsum(log(1+delta1[,1]/125)))

#normalizing the scale

nonzero = delta1[,2] - c(0,delta1[1:(nrow(delta1)-1),2]) != 0
delta1_nonzero = delta1[nonzero,]
d1 = lm(delta1_nonzero[,2]~ delta1_nonzero[,1])$coefficients[2]
d2 = lm(delta2[,2]~ delta2[,1])$coefficients[2]
d4 = lm(delta4[,2]~ delta4[,1])$coefficients[2]
do = lm(delta_oracle[nonzero,2]~delta_oracle[nonzero,1])$coefficients[2]

delta1[,2] = delta1[,2]/d1
delta2[,2] = delta2[,2]/d2
delta4[,2] = delta4[,2]/d4
delta_oracle[,2] = delta_oracle[,2]/do

plot(delta1, type = 'l', lwd =2 , col ='darkgreen')
lines(delta2, lwd =2 , col ='steelblue')
lines(delta4, lwd = 2, col = 'brown1')
lines(delta_oracle, lwd = 3, lty = 2, col = 'pink')
grid()
legend('topleft',c('Feed-Forward', "Coxph shared", "EM"), col = c('darkgreen', 'steelblue', 'brown1'), lty = c(1,1,1), lwd = c(2,2,2))
