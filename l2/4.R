library(dplyr)
library(ggplot2)
library(survival)
library(qqplotr)
library(gridExtra)

dados = read.table("~/Documentos/Survival/l2/Lista2_cancerLingua.txt")
colnames(dados) = c("perfil", "tempo", "status")
dados = dados[order(dados$tempo), ]

dados_an = subset(dados, perfil == 1)
ajuste_km = survfit(Surv(tempo, status) ~ perfil, data=dados_an)
plot(ajuste_km, xlab = "tempo (meses)", ylab = "s(t)", main = "kaplan-meier para perfil aneuploidia")

summary(ajuste_km, times = 60)

# taxa de falha acumulada
tabela = data.frame(
  tempo = ajuste_km$time,
  d = ajuste_km$n.event,
  n = ajuste_km$n.risk,
  d_n = ajuste_km$n.event / ajuste_km$n.risk,
  lambda_acumulada = cumsum(ajuste_km$n.event / ajuste_km$n.risk)
)

lambda_60 = sum(tabela$d_n[tabela$tempo <= 60])
erro_60 = sqrt(sum(tabela$d[tabela$tempo <= 60] / tabela$n[tabela$tempo <= 60]^2))
S_60_na = exp(-lambda_60)

summary_km = summary(ajuste_km, times = 60)

cat("S(60): (", round(summary_km$surv, 3), ",", round(summary_km$std.err, 3), ")\n")
cat("Ŝ(60): (", round(S_60_na, 3), ",", round(erro_60, 3), ")\n")


summary_km = summary(ajuste_km, times = 60)
S_hat = summary_km$surv
se_S = summary_km$std.err

loglog = log(-log(S_hat))
se_loglog = se_S / (S_hat * log(S_hat))

lower_loglog = loglog - qnorm(0.975) * se_loglog
upper_loglog = loglog + qnorm(0.975) * se_loglog

IC_inf = exp(-exp(upper_loglog))
IC_sup = exp(-exp(lower_loglog))

cat("IC 95% para S(60): (", round(IC_inf, 3), ",", round(IC_sup, 3), ")\n")

