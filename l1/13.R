library(ggplot2)
library(survival)


lambdas = c(0.5, 1)
rhos = c(0.5, 1, 2)

t = seq(0, 10, length.out=100)

S_t = function(t, lambda, rho)
{
  return(exp(-lambda * t^rho))
}

lambda_t = function(t, lambda, rho)
{
  return(lambda * rho * t^{rho - 1})
}

df = expand.grid(t=t, lambda=lambdas, rho=rhos)
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

