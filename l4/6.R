library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(flexsurv)
library(maxstat)

dados = readxl::read_excel("viarapida.xls")
dados = dados[order(dados$T_CC), ]


################################################################################
# Kaplan-Meier 

dados$sexo_cat = ifelse(dados$SEXO == 0, "fem", "masc")
ajuste_km_sexo = survfit(Surv(T_CC, DELTA) ~ sexo_cat, data=dados)
summary(ajuste_km_sexo)
ggsurvplot(
  ajuste_km_sexo,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier para SEXO"
)


dados$raca_cat = ifelse(dados$RACA == 1, "branco", ifelse(dados$RACA == 2, "negro", "amarelo"))
ajuste_km_raca = survfit(Surv(T_CC, DELTA) ~ raca_cat, data=dados)
summary(ajuste_km_raca)
ggsurvplot(
  ajuste_km_raca,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier para RACA"
)


dados$tipac_cat = ifelse(dados$TIPAC == 0, "congenito", "coronariano")
ajuste_km_tipac = survfit(Surv(T_CC, DELTA) ~ tipac_cat, data=dados)
summary(ajuste_km_tipac)
ggsurvplot(
  ajuste_km_tipac,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier para TIPAC"
)


dados$protocolo_cat = ifelse(dados$Protocolo == 0, "convencional", "via-rapida")
ajuste_km_protocolo = survfit(Surv(T_CC, DELTA) ~ protocolo_cat, data=dados)
summary(ajuste_km_protocolo)
ggsurvplot(
  ajuste_km_protocolo,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier para Protocolo"
)


################################################################################
# log-rank

logrank_sexo = survdiff(Surv(T_CC, DELTA) ~ sexo_cat, data = dados)
logrank_raca = survdiff(Surv(T_CC, DELTA) ~ raca_cat, data = dados)
logrank_tipac = survdiff(Surv(T_CC, DELTA) ~ tipac_cat, data = dados)
logrank_protocolo = survdiff(Surv(T_CC, DELTA) ~ protocolo_cat, data = dados)

logrank_sexo
logrank_raca
logrank_tipac
logrank_protocolo


################################################################################
# taxa de falha acumulada

# como testar exponencial, Weibull ou log-normal?
ajuste_km_protocolo = survfit(Surv(T_CC, DELTA) ~ protocolo_cat, data=dados)
ggsurvplot(
  ajuste_km_protocolo,
  fun = "cumhaz",
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "ajuste de Nelson-Aalen para Tipo de transplante"
)


################################################################################
# cox-snell

fit_flex_ex = flexsurvreg(Surv(T_CC, DELTA) ~ protocolo_cat, data=dados, dist="exponential")
cs_ex = coxsnell_flexsurvreg(fit_flex_ex)

fit_flex_wb = flexsurvreg(Surv(T_CC, DELTA) ~ protocolo_cat, data=dados, dist="weibull")
cs_wb = coxsnell_flexsurvreg(fit_flex_wb)

fit_flex_ln = flexsurvreg(Surv(T_CC, DELTA) ~ protocolo_cat, data=dados, dist="lognormal")
cs_ln = coxsnell_flexsurvreg(fit_flex_ln)


dados$cs_ex = cs_ex$est
dados$cs_wb = cs_wb$est
dados$cs_ln = cs_ln$est

g = ggplot(dados) +
  geom_line(aes(x=T_CC, y=cs_ex, colour="exponential"), linewidth=1) +
  geom_line(aes(x=T_CC, y=cs_wb, colour="weibull"), linewidth=1) +
  geom_line(aes(x=T_CC, y=cs_ln, colour="lognormal"), linewidth=1) +
  labs(title="residuos") +
  theme_minimal()
print(g)
