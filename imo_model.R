# Ładowanie potrzebnych bibliotek
library(readr)
library(car)
library(carData)
library(nortest)
library(lmtest)
library(sur)
library(MASS)

# Wczytanie danych z pliku CSV do zmiennej Mydata
Mydata=read.csv("imo_data.csv")

# Wyświetlenie podsumowania danych w Mydata
summary(Mydata)

# Tworzenie modelu regresji liniowej (imo2015_tasks_done jako zmienna objaśniana)
model=lm(imo2015_tasks_done ~ GDP + population + high_tech_exports + migration + 
           gov_educ_expenditure + area + prob_of_death + internet_user_percentage + 
           gross_enrollment_ratio + unemployment_rate, data=Mydata)

# Wyświetlenie podsumowania modelu
summary(model)

# Tworzenie macierzy X z zmiennymi objaśniającymi (bez pierwszych dwóch kolumn)
X=Mydata[,-c(1, 2)]

# Sprawdzanie współliniowości zmiennych objaśniających modelu
vif(model)
# Obliczanie pierwiastka z największej wartości własnej macierzy X'X
sqrt(max(eigen(t(as.matrix(X))%*%as.matrix(X))$values) / eigen(t(as.matrix(X))%*%as.matrix(X))$values)
# Sprawdzanie korelacji między zmiennymi objaśniającymi
cor(X, use="pairwise.complete.obs")

# Testowanie normalności reszt modelu
par(mfrow=c(1, 2))
hist(model$residuals, main='Histogram rozkładu residułów', xlab='Wartości residułów', ylab='Częstość występowania', prob=TRUE)
lines(density(model$residuals))
qqPlot(model$residuals, xlab='Teoretyczne residuła', ylab='Rzeczywistę residuała', main='Q-Q Plot')
# Testowanie normalności reszt przy użyciu testów
cat(ad.test(model$residuals)$p.value, shapiro.test(model$residuals)$p.value, lillie.test(model$residuals)$p.value)

# Homoskedastyczność
bptest(model)
gqtest(model)
hmctest(model)
par(mfrow=c(1, 1))
plot(model, 1)

# Testowanie autokorelacji w resztach modelu
dwtest(model)
acf(model$residuals, main='Wykres autokorelacji residuów')
Box.test(model$residuals, type='Ljung-Box')

# Testowanie liniowej struktury modelu
raintest(model)
resettest(model)
harvtest(model)

#---------------------------------------------------
#2 częsć


# Identyfikacja nietypowych obserwacji na podstawie dźwigni i wartości Cooka
leverage(model)
cooks.distance(model)

# Wyszukanie obserwacji o wysokiej dźwigni (high_leverage) i wartości Cooka (high_CD)
high_leverage=which(leverage(model) > mean(leverage(model)) + 2*sd(leverage(model)))
high_CD=as.vector(which(cooks.distance(model) > mean(cooks.distance(model)) + 2*sd(cooks.distance(model))))
high_leverage
high_CD

# Redukcja nieistotnych zmiennych za pomocą metody "backward elimination"
backward_elimination=step(model, direction = "backward")

# Redukcja nieistotnych zmiennych za pomocą metod AIC (Akaike Information Criterion) i BIC (Bayesian Information Criterion)
aic=step(model, direction = "both", k = 2)
bic=step(model, direction = "both", k = log(nrow(Mydata)), criteria = "BIC")

# Tworzenie nowego modelu zawierającego tylko zmienne statystycznie istotne z wcześniejszego modelu
# high_tech_exports, area, prob_of_death, internet_user_percentage
NewModel=lm(imo2015_tasks_done ~ high_tech_exports + area + prob_of_death +internet_user_percentage, data=Mydata)
summary(NewModel)
cat(ad.test(NewModel$residuals)$p.value, shapiro.test(NewModel$residuals)$p.value)

# Tworzenie trzeciego modelu ze zmodyfikowanymi zmiennymi
NewModel3=lm(imo2015_tasks_done ~ high_tech_exports + I(area^(1/2)) + prob_of_death + I(exp(internet_user_percentage)), data=Mydata)
summary(NewModel3)

# Testowanie założeń trzeciego modelu
# Współliniowość
vif(NewModel3)

# Testowanie normalności reszt modelu
par(mfrow=c(1, 2))
hist(NewModel3$residuals, main='Histogram rozkładu residułów', xlab='Wartości residułów', ylab='Częstość występowania', prob=TRUE)
lines(density(NewModel3$residuals))
qqPlot(NewModel3$residuals, xlab='Teoretyczne residuła', ylab='Rzeczywistę residuała', main='Q-Q Plot')
# Testowanie normalności reszt przy testów
cat(ad.test(NewModel3$residuals)$p.value, shapiro.test(NewModel3$residuals)$p.value, lillie.test(NewModel3$residuals)$p.value)

# Homoskedastyczność
bptest(NewModel3)
gqtest(NewModel3)
hmctest(NewModel3)
par(mfrow=c(1, 1))
plot(NewModel3, 1)

# Testowanie autokorelacji w resztach modelu
dwtest(NewModel3)
acf(NewModel3$residuals, main='Wykres autokorelacji residuów')
Box.test(NewModel3$residuals, type='Ljung-Box')

# Testowanie liniowej struktury modelu
raintest(NewModel3)
resettest(NewModel3)
harvtest(NewModel3)

# Próba przejścia do czwartego modelu z transformacją Boxa-Coxa
NewModel4=lm(imo2015_tasks_done ~ high_tech_exports + I(area^(1/2)) + prob_of_death + I(exp(internet_user_percentage)), data=Mydata[-c(56),])
# Odrzucamy 56 obserwację, bo imo2015_tasks_done[56]=0, a zamierzamy zastosować transformację Boxa-Coxa
boxcox(NewModel4) # Wartość 1 należy do przedziału ufności, więc nie zmieniamy naszego modelu


# Porównywanie modeli na podstawie różnych wskaźników
cat(summary(model)$r.squared, summary(NewModel)$r.squared, summary(NewModel3)$r.squared)
cat(summary(model)$adj.r.squared, summary(NewModel)$adj.r.squared, summary(NewModel3)$adj.r.squared)
cat(AIC(model), AIC(NewModel), AIC(NewModel3))
cat(BIC(model), BIC(NewModel), BIC(NewModel3))

# Walidacja krzyżowa
n=nrow(Mydata)
# Będziemy robić histogramy wartości MSE (Mean Squared Error) dla tych modeli
wartościMSE1=c() # MSE dla modelu
wartościMSENew=c() # MSE dla NewModel
wartościMSENew3=c() # MSE dla NewModel3

for (i in 1:1000) {
  S=sample(n, 0.8*n) # Losowa próbka o długości 80% jako training set
  training=Mydata[S,] # Training data set
  testing=Mydata[-S,] # Testing data set (wszystko oprócz S)
  
  # Tworzenie modeli na podstawie danych treningowych
  model_1=lm(imo2015_tasks_done ~ GDP + population + high_tech_exports + migration + gov_educ_expenditure + area + prob_of_death + internet_user_percentage + gross_enrollment_ratio + unemployment_rate, data=training)
  NewModel_1=lm(imo2015_tasks_done ~ high_tech_exports + area + prob_of_death +internet_user_percentage, data=training)
  NewModel3_1=lm(imo2015_tasks_done ~ high_tech_exports + I(area^(1/2)) + prob_of_death + I(exp(internet_user_percentage)), data=Mydata)
  
  # Obliczanie wartości MSE dla każdego modelu
  MSE_1=mean((predict(model_1, testing)-Mydata$imo2015_tasks_done)^2)
  MSE_New=mean((predict(NewModel_1, testing)-Mydata$imo2015_tasks_done)^2)
  MSE_New3=mean((predict(NewModel3_1, testing)-Mydata$imo2015_tasks_done)^2)
  
  wartościMSE1=append(wartościMSE1, MSE_1)
  wartościMSENew=append(wartościMSENew, MSE_New)
  wartościMSENew3=append(wartościMSENew3, MSE_New3)
}

# Histogramy dla wartości MSE
par(mfrow=c(1, 3))
hist(wartościMSE1) 
hist(wartościMSENew)
hist(wartościMSENew3)
cat(mean(wartościMSE1), mean(wartościMSENew), mean(wartościMSENew3))
# Na podstawie histogramu możemy stwierdzić, że trzeci model jest lepszy
