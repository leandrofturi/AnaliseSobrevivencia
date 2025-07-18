library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(maxstat)
library(qqplotr)
library(gridExtra)

dados = readxl::read_excel("Lista2_Hodgkins.xlsx")
dados = dados[order(dados$survivaltime), ]

################################################################################
# sexo

dados$sex = factor(dados$sex, labels = c("fem", "masc"))

ajuste_km_sexo <- survfit(Surv(survivaltime, dead) ~ sex, data = dados)

plot(ajuste_km_sexo, xlab = "tempo", ylab = "s(t)", main = "kaplan-meier")

summary(ajuste_km_sexo)

ggsurvplot(
  ajuste_km_sexo,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier por sexo"
)

# visualmente, as curvas de sobrevivência para homens e mulheres são muito 
# semelhantes ao longo do tempo. ambas apresentam quedas graduais e estão 
# acompanhadas por intervalos de confiança sobrepostos, indicando que não há 
# diferenças significativas na sobrevivência entre os sexos

logrank_test = survdiff(Surv(survivaltime, dead) ~ sex, data = dados)
logrank_test

wilcoxon_test = survdiff(Surv(survivaltime, dead) ~ sex, data = dados, rho=1)
wilcoxon_test

# os testes tem como hipótese nula que as funções de sobrevivência são iguais 
# para os dois grupos. em ambos os casos, os valores de p (0.9 e 0.8) são muito 
# altos, não sendo possível rejeitar a hipótese nula. ou seja, não há evidência 
# estatística de que a sobrevivência dos pacientes diferencie-se por sexo.


################################################################################
# idade

dados$surv = with(dados, Surv(survivaltime, dead))

# maximiza a estatística do teste log-rank, separando os grupos com maior 
# diferença na sobrevivência
teste_maxstat = maxstat.test(surv ~ age, data = dados, smethod = "LogRank", pmethod = "exactGauss")
teste_maxstat$estimate

plot(teste_maxstat)

ponto_corte = teste_maxstat$estimate
dados$idade_cat = ifelse(dados$age <= ponto_corte, "novos", "velhos")

ajuste_idade = survfit(surv ~ idade_cat, data = dados)
ggsurvplot(
  ajuste_idade,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier por idade"
)

# vemos que os pacientes mais velhos tem uma curva de sobrevivência 
# inferior à dos mais novos. a separação entre as curvas é clara, com 
# intervalos de confiança que não se sobrepõem na maior parte do tempo.

logrank_test = survdiff(Surv(survivaltime, dead) ~ idade_cat, data = dados)
logrank_test

wilcoxon_test = survdiff(Surv(survivaltime, dead) ~ idade_cat, data = dados, rho=1)
wilcoxon_test

# ambos os testes rejeitam a hipótese nula de igualdade entre as curvas de 
# sobrevivência, com alto grau de significância. isso sugere que a idade exerce 
# forte influência na sobrevivência dos pacientes neste conjunto de dados

################################################################################
# idade

dados$idade_cat2 = ifelse(dados$age < 25, "x<25", 
                          ifelse(dados$age < 38, "25<=x<38",
                                 ifelse(dados$age < 53, "38<=x<53", ">=53")))

ajuste_idade2 = survfit(surv ~ idade_cat2, data = dados)
ggsurvplot(
  ajuste_idade2,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier por idade"
)

# o grupo com > 53 anos apresenta sobrevida visivelmente mais baixa. os demais 
# grupos apresentam curvas mais próximas entre si. há sobreposição dos 
# intervalos de confiança em alguns grupos intermediários, indicando semelhanças

logrank_test = survdiff(Surv(survivaltime, dead) ~ idade_cat2, data = dados)
logrank_test

wilcoxon_test = survdiff(Surv(survivaltime, dead) ~ idade_cat2, data = dados, rho=1)
wilcoxon_test

# com o valor de p alto, ambos rejeitam a hipótese nula de igualdade entre os 
# grupos, reforçando a conclusão gráfica, com ênfase em diferenças nos tempos 
# iniciais de falha

################################################################################
# estagio doenca

dados$stage = factor(dados$stage, labels = c("inicial", "avancado"))

ajuste_stage = survfit(surv ~ stage, data = dados)
ggsurvplot(
  ajuste_stage,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier por stage"
)

logrank_test = survdiff(Surv(survivaltime, dead) ~ stage, data = dados)
logrank_test

wilcoxon_test = survdiff(Surv(survivaltime, dead) ~ stage, data = dados, rho=1)
wilcoxon_test

# baseado nas curvas de kaplan-meier e nos testes log-rank e wilcoxon, não há 
# diferença significativa na sobrevivência entre pacientes em estágio inicial 
# e avançado da doença nessa amostra. os resultados visuais (curvas sobrepostas)
# e numéricos (alto p valor) reforçam

################################################################################
# histologia

dados$hist = factor(dados$hist, labels = c("escl nod", "misto cel", "depl de linf"))

ajuste_hist = survfit(surv ~ hist, data = dados)
ggsurvplot(
  ajuste_hist,
  data = dados,
  conf.int = TRUE,
  pval = TRUE,
  title = "curvas de Kaplan-Meier por hist"
)

logrank_test = survdiff(Surv(survivaltime, dead) ~ hist, data = dados)
logrank_test

wilcoxon_test = survdiff(Surv(survivaltime, dead) ~ hist, data = dados, rho=1)
wilcoxon_test

# a curva do grupo depleção de linfócitos se destaca por apresentar menor 
# probabilidade de sobrevivência ao longo do tempo. os grupos esclerose nodular 
# e celularidade mista têm curvas mais próximas entre si e apresentam melhor 
# sobrevida. ambos os testes rejeitam a hipótese nula de igualdade entre as 
# curvas de sobrevivência, ou seja, para a variável histologia, há forte 
# evidência de diferença entre os grupos, com pior sobrevida para pacientes com 
# o subtipo de depleção de linfócitos
