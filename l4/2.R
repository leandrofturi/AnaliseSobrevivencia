library(dplyr)
library(survival)
library(survminer)
library(maxstat)


dados = read.csv("HOD_NHL.csv", sep=";")
dados = dados[order(dados$Time), ]


################################################################################
# Kaplan-Meier 

dados$graft_cat = ifelse(dados$Graft == 1, "alogenico", "autologo")
ajuste_km_graft = survfit(Surv(Time, D_R) ~ graft_cat, data=dados)
summary(ajuste_km_graft)
ggsurvplot(
  ajuste_km_graft,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier para Tipo de transplante"
)


dados$karnofsky_cat = ifelse(dados$Karnofsky < 80, "<80", ">=80")
ajuste_km_karnofsky = survfit(Surv(Time, D_R) ~ karnofsky_cat, data=dados)
summary(ajuste_km_karnofsky)
ggsurvplot(
  ajuste_km_karnofsky,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier para Escore de Karnofsk"
)

dados$surv = with(dados, Surv(Time, D_R))
teste_maxstat = maxstat.test(surv ~ Karnofsky, data = dados, smethod = "LogRank", pmethod = "exactGauss")
teste_maxstat$estimate


################################################################################
# ajustes via survreg

fit_weib = survreg(Surv(Time, D_R) ~ graft_cat + karnofsky_cat, data=dados, dist="weibull")
summary(fit_weib)


fit_llogis = survreg(Surv(Time, D_R) ~ graft_cat + karnofsky_cat, data=dados, dist="loglogistic")
summary(fit_llogis)


################################################################################
# razao de taxas de falha
# modelo: Y = log(T) = x beta + sigma epsilon
# HR = exp(-beta/sigma)

coefs = coef(fit_weib)
scale = fit_weib$scale
HR_graft = exp(-coefs["graft_catautologo"] / scale)
HR_karnofsky = exp(-coefs["karnofsky_cat>=80"] / scale)
HR_graft
HR_karnofsky

coefs = coef(fit_llogis)
scale = fit_llogis$scale
HR_graft = exp(-coefs["graft_catautologo"] / scale)
HR_karnofsky = exp(-coefs["karnofsky_cat>=80"] / scale)
HR_graft
HR_karnofsky


################################################################################
# estatistica de Wald e p-valor

wald_test = summary(fit_weib)$table
wald_test
# compare com alfa = 0.10
# se p < 0.10, rejeitamos H0 (diferenca significativa)
# se p >= 0.10, nao rejeitamos H0


################################################################################
# Nelsol-Aalen

ggsurvplot(
  ajuste_km_graft,
  fun = "cumhaz",
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "ajuste de Nelson-Aalen para Tipo de transplante"
)

ggsurvplot(
  ajuste_km_karnofsky,
  fun = "cumhaz",
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "ajuste de Nelson-Aalen para Escore de Karnofsk"
)
