library(survival)
library(survminer)
library(flexsurv)
library(dplyr)
library(ggplot2)


dados <- read.table("Lista2_cancerLingua.txt")
colnames(dados) <- c("perfil", "time", "status")
df <- dados[order(dados$time), ]
df <- df %>%
  mutate(perfil = ifelse(perfil == 1, "aneuploidia", "normal")) %>%
  mutate(perfil = factor(perfil, levels = c("aneuploidia", "normal"))) %>%
  mutate(perfil = relevel(perfil, "normal"))

################################################################################
# Kaplan-Meier e Teste de Log-rank
################################################################################

surv_obj <- Surv(time = df$time, event = df$status)

fit_km <- survfit(surv_obj ~ perfil, data = df)

# curvas KM
gg <- ggsurvplot(
  fit_km,
  data = df,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = TRUE,  
  legend.title = "Perfil",
  legend.labs = c("Diploide (normal)", "Aneuploidia"),
  xlab = "Tempo (dias)",
  ylab = "Probabilidade de sobrevivência",
  ggtheme = theme_minimal(base_size = 13)
)
print(gg)

# log-rank
# se o p-valor do log-rank < 0.05, rejeitamos a hipótese nula de que as curvas são iguais, indicando efeito do tratamento
logrank <- survdiff(surv_obj ~ perfil, data = df)
logrank

# Chisq= 2.8  on 1 degrees of freedom, p= 0.09
# como o p-valor é 0.09 (> 0.05), não rejeitamos H0 ao nível de 5%
# ou seja, não há evidência estatistica de diferença entre os dois grupos testados


################################################################################
# modelo de Cox com Z
################################################################################
# Z = 1 se aneuploidia; referência = normal
df <- df %>%
  mutate(
    Z = as.integer(perfil == "aneuploidia")
  )

# métodos de empates: Breslow, Efron, Exato
fit_breslow <- coxph(surv_obj ~ Z, data = df, ties = "breslow")
fit_efron <- coxph(surv_obj ~ Z, data = df, ties = "efron")
fit_exact <- coxph(surv_obj ~ Z, data = df, ties = "exact")

summary(fit_breslow)
summary(fit_efron)
summary(fit_exact)

# HR para Z: efeito da aneuploidia vs diploide
# HR = exp((1-0)*beta) = exp(coef)
# Breslow: HR = 0.6307, p = 0.1
# Efron: HR = 0.6273 p = 0.0963
# Exato: HR = 0.6260, p = 0.0977
# o HR está em torno de 0.6 em todos os métodos, ou seja, a aneuploidia está associada a 40% menor risco em comparação com diploidia
# se olharmos só o HR: pacientes com aneuploidia parecem ter risco menor
# porém, o p-valor ~ 0.09 é maior que 0.05: não há evidência estatisticamente significativa ao nível de 5%
# isso significa que não podemos concluir formalmente que a aneuploidia altera o prognóstico


################################################################################
# resíduos de Cox-Snell
################################################################################

fit <- fit_breslow
resm <- resid(fit, type="martingale")
res <- df$status - resm # resíduos de Cox-Snell

fit_res <- survfit(Surv(res, df$status) ~ 1)
res <- sort(res)
exp1 <- exp(-res)

gg <- ggplot() +
  geom_step(aes(x = fit_res$time, y = fit_res$surv, color = "Kaplan–Meier"), size = 1) +
  geom_line(aes(x = res, y = exp1, color = "Exponencial padrão"), linetype = "dashed", size = 1) +
  labs(
    x = "Resíduos de Cox–Snell",
    y = "S(e)",
    title = "Kaplan–Meier dos resíduos de Cox–Snell",
    color = "Curva"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1))
print(gg)

sexp <- exp(-fit_res$time)

gg <- ggplot() +
  geom_point(aes(x = fit_res$surv, y = sexp), size = 1) +
  geom_abline() +
  labs(
    x = "S(e): Kaplan-Meier",
    y = "S(e): exponencial padrão",
    title = "Kaplan–Meier dos resíduos de Cox–Snell",
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1))
print(gg)

# se modelo de Cox está bem ajustado, então os resíduos de Cox–Snell devem seguir uma Exponencial(1)
# a curva do Kaplan–Meier nos resíduos acompanha bem a linha vermelha tracejada da Exp(1)
# conclusão: o modelo de Cox ajustado parece adequado para esses dados



################################################################################
# resíduos Deviance
################################################################################

fit <- fit_breslow
res <- resid(fit, type="deviance")
pl <- fit$linear.predictors

gg <- ggplot() +
  geom_point(aes(x = pl, y = res), color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = c(-3, 3), linetype = 3, color = "red") + # 3 do livro do Colosimo
  labs(
    x = "Preditor linear",
    y = "Resíduo deviance",
    title = "Resíduos deviance do modelo de Cox"
  ) +
  theme_minimal(base_size = 14)
plot(gg)

# a maioria dos pontos está concentrada entre -2 e +2
# alguns chegam perto de +=2,5, mas nenhum passa de 3
# isso indica que não há outliers


################################################################################
# avaliação da suposição de proporcionalidade
################################################################################

# pl <- fit$linear.predictors
# risk <- exp(pl)
# t_events <- sort(unique(df$time[df$status == 1]))
# d_at_tj <- function(tj) sum(df$status == 1 & df$time == tj) # d_j: n de falhas em t_j
# denom_at_tj <- function(tj) sum(risk[df$time >= tj]) # risco R(t_j)

# breslow_steps <- lapply(t_events, function(tj) {
#   dj <- d_at_tj(tj)
#   denom <- denom_at_tj(tj)
#   dLambda0 <- dj / denom
#   data.frame(tj = tj, dj = dj, denom = denom, dLambda0 = dLambda0)
# }) %>% bind_rows() %>%
#   mutate(Lambda0 = cumsum(dLambda0)) # soma cumulativa

fit <- coxph(Surv(time, status) ~ strata(perfil), data = df, ties = "breslow")
bh <- basehaz(fit, centered = FALSE)

gg <- ggplot() +
  geom_step(aes(x = log(bh$time), y = log(bh$hazard), color = bh$strata), linewidth = 1) +
  labs(x = "log(t)", y = "log(Lambda0j(t))",
       title = "Método gráfico descritivo",
       color = "Grupo") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
print(gg)

# suposição de riscos proporcionais: as curvas devem ser aproximadamente paralelas ao longo do tempo
# as duas curvas sobem de forma semelhante, sem cruzamentos grandes
# porém, há momentos em que a curva normal se afasta mais da aneuploidia, principalmente em tempos intermediários, indicando que a diferença de riscos não é exatamente constante


# teste de proporcionalidade dos riscos (Schoenfeld)
fit <- coxph(Surv(time, status) ~ perfil, data = df, ties = "breslow")
ph <- cox.zph(fit, transform = "identity")
print(ph) 

plot(ph, xlab = "tempo", ylab = "Beta(t) para Aneuploidia", main = "Schoenfeld")

# tanto o teste global quanto o teste para a covariável indicam evidências a favor da hipótese nula de taxas de falha proporcionais (p > 0.05)
# o grafico corrobora com este fato, visto que não há tendências marcantes das curvas suavizadas ao longo do tempo, o que implica ser razoável considerar que a inclinação das curvas é nula
# logo, há evidências a favor da suposição de taxas de falhas proporcionais


################################################################################
# Weibull vs Cox
################################################################################

fit_weib <- survreg(Surv(time, status) ~ perfil, data = df, dist = "weibull")
summary(fit_weib)

beta <- fit_weib$coefficients[2]
HR <- beta / fit_weib$scale; HR


# resíduos de Cox-Snell
fit_flex = flexsurvreg(Surv(time, status) ~ perfil, data = df, dist = "weibull")
res = coxsnell_flexsurvreg(fit_flex) # resíduos de Cox-Snell

fit_res <- survfit(Surv(res$est, res$status) ~ 1)
res <- sort(res$est)
exp1 <- exp(-res)

gg <- ggplot() +
  geom_step(aes(x = fit_res$time, y = fit_res$surv, color = "Kaplan–Meier"), size = 1) +
  geom_line(aes(x = res, y = exp1, color = "Exponencial padrão"), linetype = "dashed", size = 1) +
  labs(
    x = "Resíduos de Cox–Snell",
    y = "S(e)",
    title = "Kaplan–Meier dos resíduos de Cox–Snell",
    color = "Curva"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1))
print(gg)

sexp <- exp(-fit_res$time)

gg <- ggplot() +
  geom_point(aes(x = fit_res$surv, y = sexp), size = 1) +
  geom_abline() +
  labs(
    x = "S(e): Kaplan-Meier",
    y = "S(e): exponencial padrão",
    title = "Kaplan–Meier dos resíduos de Cox–Snell",
    color = "Curva"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1))
print(gg)

# a curva do Kaplan–Meier nos resíduos acompanha bem a linha vermelha tracejada da Exp(1)
# conclusão: o modelo de Cox ajustado parece adequado para esses dados

# comparação
# modelo de Cox (semi paramétrico):
# - não assume forma funcional para o risco basal
# - mais robusto, mas menos eficiente se o verdadeiro modelo for de fato Weibull
# - HR é constante por construção

# modelo Weibull (paramétrico):
# - supõe explicitamente distribuição
# - se os dados se ajustarem bem, o modelo é mais eficiente (intervalos de confiança mais estreitos, permite extrapolar sobrevivência)

fit_cox <- coxph(Surv(time, status) ~ perfil, data = df, ties = "breslow")
surv_cox = survfit(fit_cox)

fit_weib <-  survreg(Surv(time, status) ~ perfil, data = df, dist = "weibull")
alpha <- exp(fit_weib$coefficients[1])
gama <- 1/fit_weib$scale
surv_weib <- exp(-(df$time/alpha)^gama)

gg <- ggplot() +
  geom_line(aes(surv_cox$time, surv_cox$surv, color = "Cox"), linewidth = 1) +
  geom_line(aes(df$time, surv_weib, color = "Weibull"), linewidth = 1) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "tempo", y = "S(t)",
       title = "Cox vs. Weibull",
       color = "") +
  theme_minimal(base_size = 14)
print(gg)

# como os dados se ajustaram bem, o modelo Weibull é mais eficiente (intervalos de confiança mais estreitos, permite extrapolar sobrevivência)
