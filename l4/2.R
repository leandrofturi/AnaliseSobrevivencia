library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(flexsurv)
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


################################################################################
# razao de taxas de falha
# weibull: exp{(xk - xj)beta}

beta = coef(fit_weib)["graft_catautologo"]
# x_j = 0 (alogenico), x_k = 1 (autologo)
HR = exp((1 - 0) * beta); HR


beta = coef(fit_weib)["karnofsky_cat>=80"]
# x_j = 0 (<80), x_k = 1 (>=80)
HR = exp((1 - 0) * beta); HR


# OU
shape = 1 / fit_weib$scale
beta = -coef(fit_weib) * shape
HR = exp(beta["graft_catautologo"]); HR


################################################################################
# estatistica de Wald e p-valor

wald_test = summary(fit_weib)$table
wald_test
# compare com alfa = 0.10
# se p < 0.10, rejeitamos H0 (diferenca significativa)
# se p >= 0.10, nao rejeitamos H0


# IGUAL
se_beta = sqrt(diag(fit_weib$var))[c("(Intercept)", "graft_catautologo", "karnofsky_cat>=80")]
wald_z = fit_weib$coefficients[c("(Intercept)", "graft_catautologo", "karnofsky_cat>=80")] / se_beta
p_values = 2 * (1 - pnorm(abs(wald_z))); p_values


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


# estimativas da funcao de taxa de falha acumulada?
fit_graft = survreg(Surv(Time, D_R) ~ graft_cat, data=dados, dist="weibull")

beta0 = fit_graft$coefficients["(Intercept)"]
beta = fit_graft$coefficients["graft_catautologo"]
sigma = fit_graft$scale
shape_w = 1 / sigma; shape_w
scale_w = exp(beta0); scale_w


################################################################################
# residuos de Cox-Snell

ei = (dados$Time*exp(-beta0 -beta * dados$Graft)) ^ shape_w


# OU
fit_flex = flexsurvreg(Surv(Time, D_R) ~ graft_cat, data=dados, dist="weibull")
cs = coxsnell_flexsurvreg(fit_flex)

qy = qexp(ppoints(nrow(cs),0))
qqplot(qy, cs$est)
abline(a=0, b=1, col="red", lwd=2)


# OU https://stats.stackexchange.com/questions/246812/cox-snell-residuals-in-r
fit = survreg(Surv(Time, D_R) ~ graft_cat, data=dados, dist="weibull")
rC = exp(((fit$y[,1]) - log(predict(fit, dados))) / fit$scale)
rC = rC + (1 - fit$y[,2]) * 1


################################################################################
# deviance

mi = dados$D_R - cs$est
di = sign(mi) * sqrt(-2*(mi + dados$D_R*log(dados$D_R - mi)))


# OU
cd = residuals(fit, type = "deviance") # martingale


# plot

dados$qy = qy
dados$cs = cs$est
dados$cd = cd

g = ggplot(dados) +
  geom_line(aes(x=qy, y=cs, colour="cox-snell"), linewidth=1) +
  geom_line(aes(x=qy, y=cd, colour="deviance"), linewidth=1) +
  labs(title="residuos") +
  theme_minimal()
print(g)
