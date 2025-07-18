library(ggplot2)
library(survival)
library(qqplotr)
library(gridExtra)


lambda = 1
rho = 0.5

t = seq(0, 10, length.out=100)

S_t = function(t, lambda, rho)
{
  return(exp(-lambda * t^rho))
}

lambda_t = function(t, lambda, rho)
{
  return(lambda * rho * t^{rho - 1})
}

df = expand.grid(t=t, lambda=lambda, rho=rho)
df$S_t = apply(df, 1, function(x) S_t(x["t"], x["lambda"], x["rho"]))
df$lambda_t = apply(df, 1, function(x) lambda_t(x["t"], x["lambda"], x["rho"]))

df$params = apply(df, 1, function(x) paste0("lambda = ", x["lambda"], " rho = ", x["rho"]))

g = ggplot(df, aes(x=t, y = S_t, color=params)) +
  geom_line(linewidth=1) +
  labs(title="funcao de sobrevivencia") +
  theme_minimal()
print(g)

g2 = ggplot(df, aes(x=t, y=lambda_t, color=params)) +
  geom_line(linewidth=1) +
  labs(title="taxa de falha") +
  theme_minimal()
print(g2)

n = 100
shape = 0.5
scale = 1.0
dfwb = data.frame(wb=rweibull(n, shape=shape, scale=scale))

g3 = ggplot(dfwb, aes(x="", y=wb)) +
  geom_boxplot() +
  labs(title="weibull") +
  theme_minimal()
print(g3)

g4 = ggplot(dfwb, aes(x=wb)) +
  geom_histogram(aes(y=..density..), bins=30, fill="skyblue", color="black") +
  stat_function(fun=rweibull, args=list(shape=shape, scale=scale), color = "red") +
  labs(title="histograma com densidade teorica")
print(g4)

km = survfit(Surv(dfwb$wb) ~ 1)
surv_df = data.frame(
  time=km$time,
  surv=km$surv,
  theor=1 - rweibull(km$time, shape=shape, scale=scale)
)

g5 = ggplot(surv_df) +
  geom_step(aes(x=time, y=surv), color="black", linetype="dashed") +
  geom_line(aes(x=time, y=theor), color="red") +
  labs(title="curva de sobrevivencia empirica vs teorica")
print(g5)

g6 = ggplot(dfwb, aes(sample=wb)) +
  stat_qq_point(distribution="weibull", dparams=c(shape, scale)) +
  stat_qq_line(distribution="weibull", dparams=c(shape, scale), color="red") +
  labs(title="qqplot: weibull")
print(g6)

dfwb$wbpad = scale(dfwb$wb)
g7 = ggplot(dfwb, aes(sample=wbpad)) +
  stat_qq_point(distribution="weibull", dparams=c(shape, scale)) +
  stat_qq_line(distribution="weibull", dparams=c(shape, scale), color="red") +
  labs(title="qqplot: weibull padronizado")
print(g7)

