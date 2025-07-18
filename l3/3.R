library(survival)


tempos <- c(0.19, 0.78, 0.96, 1.31, 2.78, 3.16, 4.67, 4.85,
            6.50, 7.35, 8.27, 12.07, 32.52, 33.91, 36.71)
# censura (10 não falharam)
# até que um número fixo de falhas, censura a direita tipo II
tempos_cens <- rep(36.71, 10)

dados <- data.frame(
  tempo = c(tempos, tempos_cens),
  status = c(rep(1, 15), rep(0, 10))  # 1 = falha, 0 = censurado
)

################################################################################
# ajustes

# log-verossimilhanças sob censura a direita:
# f para observei o evento
# S para censura

loglike_exp <- function(par) {
  mu <- par[1]
  lambda <- exp(-mu)
  with(dados, {
    sum(status * (log(lambda) - lambda * tempo) +
          (1 - status) * (-lambda * tempo))
  })
}

loglike_weib <- function(par) {
  mu <- par[1]
  sigma <- par[2]
  shape <- 1 / sigma # shape
  scale <- exp(mu) # scale
  with(dados, {
    logf <- log(shape) - shape * log(scale) + (shape - 1) * log(tempo) - (tempo / scale)^shape
    logS <- - (tempo / scale)^shape
    sum(status * logf + (1 - status) * logS)
  })
}


loglike_lnorm <- function(par) {
  mu <- par[1]
  sigma <- par[2]
  with(dados, {
    logf <- dlnorm(tempo, meanlog = mu, sdlog = sigma, log = TRUE)
    logS <- plnorm(tempo, meanlog = mu, sdlog = sigma, lower.tail = FALSE, log.p = TRUE)
    sum(status * logf + (1 - status) * logS)
  })
}


loglike_llogis <- function(par) {
  mu <- par[1]
  sigma <- par[2]
  alpha <- 1 / sigma # shape
  beta <- exp(mu) # scale
  with(dados, {
    logf <- log(alpha) - log(beta) + (alpha - 1) * log(tempo / beta) -
      2 * log(1 + (tempo / beta)^alpha)
    logS <- -log(1 + (tempo / beta)^alpha)
    sum(status * logf + (1 - status) * logS)
  })
}


start <- log(mean(dados$tempo)) # ponto inicial: peguei da net

fit_exp_optim <- optim(par = c(start), fn = function(par) -loglike_exp(par), method = "BFGS")
fit_weib_optim <- optim(par = c(start, 1), fn = function(par) -loglike_weib(par), method = "BFGS")
fit_lnorm_optim <- optim(par = c(start, 1), fn = function(par) -loglike_lnorm(par), method = "BFGS")
fit_llogis_optim <- optim(par = c(start, 1), fn = function(par) -loglike_llogis(par), method = "BFGS")

loglik_values <- -c(
  exp   = fit_exp_optim$value,
  weib  = fit_weib_optim$value,
  lnorm = fit_lnorm_optim$value,
  llog  = fit_llogis_optim$value
)
loglik_values


# AIC = -2 * loglik + 2 * k, k = número de parâmetros
aic_values <- -2 * loglik_values + c(2, 4, 4, 4)
aic_values

lambda_exp <- fit_exp_optim$par[1]

shape_weib <- fit_weib_optim$par[1]
scale_weib <- fit_weib_optim$par[2]

mu_lnorm <- fit_lnorm_optim$par[1]
sigma_lnorm <- fit_lnorm_optim$par[2]

alpha_llog <- fit_llogis_optim$par[1]
beta_llog <- fit_llogis_optim$par[2]


################################################################################
# ajustes via survreg

fit_exp <- survreg(Surv(tempo, status) ~ 1, data=dados, dist="exponential")
fit_weib <- survreg(Surv(tempo, status) ~ 1, data=dados, dist="weibull")
fit_lnorm <- survreg(Surv(tempo, status) ~ 1, data=dados, dist="lognormal")
fit_llogis <- survreg(Surv(tempo, status) ~ 1, data=dados, dist="loglogistic")


AIC(fit_exp, fit_weib, fit_lnorm, fit_llogis)
# melhor: fit_lnorm

lambda_expS <- fit_exp[[2]][1]

shape_weibS <- fit_weib[[2]][1]
scale_weibS <- fit_weib[[2]][2]

mu_lnormS <- fit_lnorm[[2]][1]
sigma_lnormS <- fit_lnorm[[2]][2]

alpha_llogS <- fit_llogis[[2]][1]
beta_llogS <- fit_llogis[[2]][2]


################################################################################
# comparação

params <- data.frame(
  expS = as.numeric(c(lambda_expS, lambda_expS)),
  exp = c(lambda_exp, lambda_exp),
  weibS = as.numeric(c(shape_weibS, exp(scale_weibS))),
  weib = c(shape_weib, scale_weib),
  lnormS = as.numeric(c(mu_lnormS, exp(sigma_lnormS))),
  lnorm = c(mu_lnorm, sigma_lnorm),
  llogisS = as.numeric(c(alpha_llogS, exp(beta_llogS))),
  llogis = c(alpha_llog, beta_llog)
)
params


################################################################################
# medianas

log(2) * exp(fit_exp$coefficients[1])
exp(fit_weib$coefficients[1]) * (log(2))^fit_weib$scale
exp(fit_lnorm$coefficients[1])
exp(fit_llogis$coefficients[1])


################################################################################
# medias

exp(fit_exp$coefficients[1])
exp(fit_weib$coefficients[1]) * gamma(1 + fit_weib$scale)
exp(fit_lnorm$coefficients[1] + 0.5 * fit_lnorm$scale^2)
exp(fit_llogis$coefficients[1]) * (pi * fit_llogis$scale) / sin(pi * fit_llogis$scale)


################################################################################
# quantis

-log(1 - 0.2) * exp(fit_exp$coefficients[1])
exp(fit_weib$coefficients[1]) * (-log(1 - 0.2))^fit_weib$scale
exp(fit_lnorm$coefficients[1] + fit_lnorm$scale * qnorm(0.2))
exp(fit_llogis$coefficients[1]) * (0.2 / 0.8)^fit_llogis$scale

