library(readr)#ładuje wszystkie biblioteki
library(car)
library(carData)
library(nortest)
library(lmtest)
library(sur)
library(MASS)


Mydata=read.csv("imo_data.csv")

summary(Mydata)

model=lm(imo2015_tasks_done ~ GDP + population + 
           + high_tech_exports + migration + gov_educ_expenditure + 
           + area + prob_of_death + internet_user_percentage + 
           + gross_enrollment_ratio + unemployment_rate, data=Mydata)


summary(model)



X=Mydata[,-c(1, 2)]#tworzę macierz z zmiennymi objaśniającymi

#współliniowość
vif(model)
M=t(as.matrix(X))%*%as.matrix(X)
sqrt(max(eigen(M)$values)/eigen(M)$values)
cor(X, use="pairwise.complete.obs")

#normalność błędów
par(mfrow=c(1, 2))
hist(model$residuals, main='Histogram rozkładu residułów', xlab='Wartości residułów', ylab='Częstość występowania', prob=TRUE)
lines(density(model$residuals))
qqPlot(model$residuals, xlab='Teoretyczne residuła', ylab='Rzeczywistę residuała', main='Q-Q Plot')
cat(ad.test(model$residuals)$p.value, shapiro.test(model$residuals)$p.value, lillie.test(model$residuals)$p.value)

#stała wariancja(problem)
bptest(model)
gqtest(model)
hmctest(model)
par(mfrow=c(1, 1))
plot(model, 1)

#autokorelacja
dwtest(model)
acf(model$residuals, main='Wykres autokorelacji residuów')
Box.test(model$residuals, type='Ljung-Box')

#liniowa struktura(problem)
raintest(model)
resettest(model)
harvtest(model)


#---------------------------------------------------
#2 częsć

#identyfikacja nietypowych obserwacji

#zmienne wpływowe oraz zmienne z wysoką dźwignią
leverage(model)
cooks.distance(model)

high_leverage=which(leverage(model) > mean(leverage(model)) + 2*sd(leverage(model)))
high_CD=as.vector(which(cooks.distance(model) > mean(cooks.distance(model)) + 2*sd(cooks.distance(model))))
high_leverage
high_CD

#redukowanie nieistotnych zmiennych
backward_elimination=step(model, direction = "backward")
aic=step(model, direction = "both", k = 2)
bic=step(model, direction = "both", k = log(nrow(Mydata)), criteria = "BIC")
#na podstawie tej informacji dochodzimy do wniosku, że możemy zrobić model używając tylko zmiennych statystycznie istotnych:
#high_tech_exports, area, prob_of_death, internet_user_percentage




#Przejście do drugiego modelu
NewModel=lm(imo2015_tasks_done ~ high_tech_exports + area + prob_of_death +internet_user_percentage, data=Mydata)
summary(NewModel)
cat(ad.test(NewModel$residuals)$p.value, shapiro.test(NewModel$residuals)$p.value)



#Przejście do trzeciego modelu
NewModel3=lm(imo2015_tasks_done ~ high_tech_exports + I(area^(1/2)) + prob_of_death + I(exp(internet_user_percentage)), data=Mydata)
summary(NewModel3)

#Sprawdzanie założeń trzeciego modelu

#współliniowość
vif(NewModel3)

#normalność błędów
par(mfrow=c(1, 2))
hist(NewModel3$residuals, main='Histogram rozkładu residułów', xlab='Wartości residułów', ylab='Częstość występowania', prob=TRUE)
lines(density(NewModel3$residuals))
qqPlot(NewModel3$residuals, xlab='Teoretyczne residuła', ylab='Rzeczywistę residuała', main='Q-Q Plot')
cat(ad.test(NewModel3$residuals)$p.value, shapiro.test(NewModel3$residuals)$p.value, lillie.test(NewModel3$residuals)$p.value)

#stała wariancja
bptest(NewModel3)
gqtest(NewModel3)
hmctest(NewModel3)

par(mfrow=c(1, 1))
plot(NewModel3, 1)

#autokorelacja
dwtest(NewModel3)
acf(NewModel3$residuals, main='Wykres autokorelacji residuów')
Box.test(NewModel3$residuals, type='Ljung-Box')

#liniowa struktura
raintest(NewModel3)
resettest(NewModel3)
harvtest(NewModel3)



#Próba przejścia do czwartego modelu
#transformacja Boxa-Coxa
NewModel4=lm(imo2015_tasks_done ~ high_tech_exports + I(area^(1/2)) + prob_of_death + I(exp(internet_user_percentage)), data=Mydata[-c(56),])
#odrzucamy 56 obserwację, bo imo2015_tasks_done[56]=0, a my będziemy chcieli zastosować transformację Boxa-Coxa
boxcox(NewModel4)#1 należy do przedziału ufności, więc nie zmieniamy nasz model




#Porównywanie modeli
cat(summary(model)$r.squared, summary(NewModel)$r.squared, summary(NewModel3)$r.squared)

cat(summary(model)$adj.r.squared, summary(NewModel)$adj.r.squared, summary(NewModel3)$adj.r.squared)

cat(AIC(model), AIC(NewModel), AIC(NewModel3))

cat(BIC(model), BIC(NewModel), BIC(NewModel3))

#walidacja krzyżowa
n=nrow(Mydata)
#będziemy za chwilę robić histogram z wartościami MSE(mean squared error) dla tych modeli
wartościMSE1=c()#mse dla model
wartościMSENew=c()#mse dla NewModel
wartościMSENew3=c()#mse dla NewModel3

for (i in 1:1000) {
  S=sample(n, 0.8*n)#losowa próbka z gala o długości 80%(to jest rozmiar training set)
  training=Mydata[S,]#training data set
  testing=Mydata[-S,]#testing data set(-S oznacza wszystko oprócz S)
  
  model_1=lm(imo2015_tasks_done ~ GDP + population + 
             + high_tech_exports + migration + gov_educ_expenditure + 
             + area + prob_of_death + internet_user_percentage + 
             + gross_enrollment_ratio + unemployment_rate, data=training)
  NewModel_1=lm(imo2015_tasks_done ~ high_tech_exports + area + prob_of_death +internet_user_percentage, data=training)
  
  NewModel3_1=lm(imo2015_tasks_done ~ high_tech_exports + I(area^(1/2)) + prob_of_death + I(exp(internet_user_percentage)), data=Mydata)

  MSE_1=mean((predict(model_1, testing)-Mydata$imo2015_tasks_done)^2)
  MSE_New=mean((predict(NewModel_1, testing)-Mydata$imo2015_tasks_done)^2)
  MSE_New3=mean((predict(NewModel3_1, testing)-Mydata$imo2015_tasks_done)^2)

  wartościMSE1=append(wartościMSE1, MSE_1)
  wartościMSENew=append(wartościMSENew, MSE_New)
  wartościMSENew3=append(wartościMSENew3, MSE_New3)
}
#histogramy dla MSE
par(mfrow=c(1, 3))
hist(wartościMSE1) 
hist(wartościMSENew)
hist(wartościMSENew3)
cat(mean(wartościMSE1), mean(wartościMSENew), mean(wartościMSENew3))
# na podstawie tego histogramu można stwierdzić, że trzeci model jest lepszy






