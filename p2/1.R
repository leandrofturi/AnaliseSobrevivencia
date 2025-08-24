library(survival)
library(survminer)
library(dplyr)
library(tibble)


df <- tribble(
  ~time, ~status, ~group,
  # CONTROLE
  20, 1, "controle",
  21, 1, "controle",
  23, 1, "controle",
  24, 1, "controle",
  24, 1, "controle",
  26, 1, "controle",
  26, 1, "controle",
  27, 1, "controle",
  28, 1, "controle",
  30, 1, "controle",
  # RADIOTERAPIA
  26, 1, "radio",
  28, 1, "radio",
  29, 1, "radio",
  29, 1, "radio",
  30, 1, "radio",
  30, 1, "radio",
  31, 1, "radio",
  31, 1, "radio",
  32, 1, "radio",
  35, 0, "radio",
  # RADIO + BPA
  31, 1, "radio_bpa",
  32, 1, "radio_bpa",
  34, 1, "radio_bpa",
  35, 1, "radio_bpa",
  36, 1, "radio_bpa",
  38, 1, "radio_bpa",
  38, 1, "radio_bpa",
  39, 1, "radio_bpa",
  42, 0, "radio_bpa",
  42, 0, "radio_bpa",
) %>%
  mutate(group = factor(group, levels = c("controle", "radio", "radio_bpa")))

################################################################################
# Kaplan-Meier e Teste de Log-rank
################################################################################

# interpretação esperada:
# controle deve apresentar menor sobrevida
# radioterapia tende a prolongar um pouco a vida dos ratos
# radioterapia + BPA deve apresentar a maior sobrevida
surv_obj <- Surv(time = df$time, event = df$status)

fit_km <- survfit(surv_obj ~ group, data = df)

# curvas KM
gg <- ggsurvplot(
  fit_km,
  data = df,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = TRUE,  
  legend.title = "Grupo",
  legend.labs = c("Controle", "Radioterapia", "Radio+BPA"),
  xlab = "Tempo (dias)",
  ylab = "Probabilidade de sobrevivência",
  ggtheme = theme_minimal(base_size = 13)
)
print(gg)

# log-rank
# se o p-valor do log-rank < 0.05, rejeitamos a hipótese nula de que as curvas são iguais, indicando efeito do tratamento
logrank <- survdiff(surv_obj ~ group, data = df)
logrank

# Chisq= 33.4  on 2 degrees of freedom, p= 6e-08
# indica uma diferença significativa entre as curvas de sobrevivência
# rejeitamos a hipótese nula de que todos os grupos têm a mesma função de sobrevivência
# grupo radioterapia + BPA apresenta menos óbitos, indicando benefício adicional da BNCT (com BPA)


################################################################################
# modelo de Cox com Z1 e Z2
################################################################################
# Z1 = 1 se radio; Z2 = 1 se radio_bpa; referência = controle
df <- df %>%
  mutate(
    Z1 = as.integer(group == "radio"),
    Z2 = as.integer(group == "radio_bpa")
  )

# métodos de empates: Breslow, Efron, Exato
fit_breslow <- coxph(surv_obj ~ Z1 + Z2, data = df, ties = "breslow")
fit_efron <- coxph(surv_obj ~ Z1 + Z2, data = df, ties = "efron")
fit_exact <- coxph(surv_obj ~ Z1 + Z2, data = df, ties = "exact")

summary(fit_breslow)
summary(fit_efron)
summary(fit_exact)

# HR para Z1: efeito da radioterapia vs controle
# HR = exp((1-0)*beta) = exp(coef)
# Breslow: HR = 0.163, p = 0.001
# Efron: HR = 0.146, p = 0.0007
# Exato: HR = 0.102, p = 0.0014
# risco de morte no grupo radioterapia é 1-0.163 = 83% a 1-0.102 = 89% menor do que no grupo controle
# efeito é estatisticamente significativo (p < 0.01 em todos os métodos)

# HR para Z2: efeito de radioterapia + BPA vs controle
# Breslow: HR = 0.0285, p = 0.000002
# Efron: HR = 0.0232, p = 0.0000009
# Exato: HR = 0.0145, p = 0.000003
# O risco de morte no grupo radioterapia + BPA é cerca de 97% menor do que no controle
# efeito é estatisticamente significativo (p < 0.01 em todos os métodos)


################################################################################
# razão de verossimilhanças
################################################################################

# coxph fit$loglik = initial values and with the final values of the coefficients
fit_full <- fit_breslow
ll0 <- fit_full$loglik[1]   # loglik do modelo sem covariáveis
ll1 <- fit_full$loglik[2]   # loglik do modelo com Z1 e Z2

lrt <- 2*(ll1 - ll0)
pval <- pchisq(lrt, df=2, lower.tail=FALSE)
c(LRT=lrt, df=2, p=pval)

# valor da estatística (27.37) é maior do que o esperado sob H0
# p-valor é pequeno (< 0.05)
# rejeitamos H0: há evidência de que os tratamentos (radioterapia e radioterapia+BPA) afetam a sobrevivência dos ratos


################################################################################
# teste de Wald
################################################################################

# W = beta^T Var(beta)^{-1} beta ~ Chi^2_df, df = 2
# matriz de informação observada de Fisher é o inverso da matriz de variância
fit_full <- fit_breslow
beta_hat <- fit_full$coefficients
V_beta <- vcov(fit_full)
W <- as.numeric(t(beta_hat) %*% solve(V_beta) %*% beta_hat) # fit_full$wald.test
pW <- pchisq(W, df = length(beta_hat), lower.tail = FALSE)

c(W=W, df=length(beta_hat), p=pW)

# valor da estatística é grande, com p < 0.05.
# rejeitamos H0: existe efeito significativo de tratamento
# tanto a radioterapia quanto a radioterapia + BPA estão associados a redução significativa no risco de morte em comparação ao controle


################################################################################
# razão de verossimilhanças
################################################################################

# comparar:
# modelo completo: group (3 níveis)
# modelo restrito: group com radio e radio_bpa no mesmo nível

df <- df %>%
  mutate(trat = ifelse(group == "controle", "controle", "tratamento")) %>%
  mutate(trat = factor(trat, levels = c("controle","tratamento")))

fit_full <- coxph(surv_obj ~ group, data = df, ties = "breslow")
summary(fit_full)

fit_restr <- coxph(surv_obj ~ trat, data = df, ties = "breslow")
summary(fit_restr)

ll_restr <- fit_restr$loglik[2]  # loglik do modelo restrito
ll_full <- fit_full$loglik[2]
lrt <- 2*(ll_full - ll_restr)
pval <- pchisq(lrt, df=1, lower.tail=FALSE)  # df=1, pois só 1 parâmetro difere
c(LRT=lrt, df=1, p=pval)

# p-valor é < 0.01, então rejeitamos H0
# significa que os efeitos da radioterapia e da radioterapia + BPA não são iguais
# o uso adicional do BPA realmente proporciona um benefício extra além da radioterapia isolada


################################################################################
# teste de Wald
################################################################################

C <- matrix(c(1, -1), nrow = 1) # (beta1 - beta2)
beta_vec <- coef(fit_full) # (beta1 - beta2)
V <- vcov(fit_full) # matriz covariância 2x2
W_contrast <- as.numeric(C %*% beta_vec) # (beta1 - beta2)
V_contrast <- as.numeric(C %*% V %*% t(C)) # Var(beta1 - beta2) = Var(beta1) + Var(beta2) − 2Cov(beta1 - beta2)
Wald_stat <- (W_contrast^2) / V_contrast # (beta1 - beta2)^2 / Var(beta1 - beta2)
p_wald_contrast <- pchisq(Wald_stat, df = 1, lower.tail = FALSE)

c(W=Wald_stat, p=p_wald_contrast)

# p-valor é < 0.01, então rejeitamos H0, ou seja, os coeficientes não são iguais
# radioterapia isolada (beta1) e radioterapia + BPA (beta2) têm efeitos significativamente diferentes sobre a sobrevivência
# como os HRs estimados mostraram que radio+BPA tem risco ainda menor que radio sozinho, concluímos que o BPA acrescenta benefício real além da radioterapia

