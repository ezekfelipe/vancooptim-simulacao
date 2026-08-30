# =============================================================================
# VancoOptim — Estudo de simulação (seções 2.8 e 3.6 do manuscrito)
#
# Avalia a recuperação da AUC24 pela estimativa bayesiana, sob mecanismo gerador
# EXTERNO à biblioteca de priors implementada na ferramenta.
#
# Estruturado pelo arcabouço ADEMP (Morris, White & Crowther, Stat Med 2019).
#
# REPRODUTIBILIDADE: semente fixa registrada abaixo. Rodar este script produz
# exatamente os mesmos números, em qualquer máquina com a mesma versão de R.
#
# Uso:   source("simulacao_vancooptim.R")
#        res <- rodar_tudo()
#
# Requisitos: R base. Nenhum pacote externo é necessário.
# =============================================================================

# ── SEMENTE ──────────────────────────────────────────────────────────────────
# Fixada para reprodutibilidade. NÃO alterar entre execuções do mesmo estudo:
# mudar a semente muda os pacientes sorteados e, portanto, os resultados.
SEMENTE <- 20260812L

# ── PARÂMETROS DO ESTUDO ─────────────────────────────────────────────────────
# Justificativa do número de réplicas pelo erro de Monte Carlo:
#   SE(estimativa) = DP(erros) / sqrt(N)
# Com DP empírico ~18 pontos percentuais para o MAPE, N = 2000 dá SE ~0,40 pp.
# Resolução suficiente para distinguir diferenças de 1 ponto entre cenários.
N_REPLICAS <- 2000L

# =============================================================================
# A — MECANISMO GERADOR: Roberts et al. 2011
#     Antimicrob Agents Chemother 55:2704-2709, Tabela 2.
#
# Externo à biblioteca do VancoOptim, o que evita avaliar o otimizador contra o
# próprio modelo que o alimenta. Totalmente especificado: efeitos fixos,
# variabilidade entre sujeitos e erro residual.
#
# RESSALVA: o modelo deriva de infusão CONTÍNUA, condição em que o volume de
# distribuição é fracamente identificável — no equilíbrio a concentração é a
# razão entre taxa e clearance, e o volume só é informado pelas transições.
# O valor de 1,53 L/kg está acima da faixa usual citada na própria publicação.
# Por isso `vd_kg` é fator do desenho, e não constante: o cenário de Roberts é
# tratado como estresse, não como caso primário.
# =============================================================================

ROBERTS <- list(
  cl_tipico   = 4.58,   # L/h, para ClCr de 100 mL/min/1,73 m²
  vd_kg       = 1.53,   # L/kg
  cv_cl       = 0.389,  # variabilidade entre sujeitos, CL
  cv_vd       = 0.374,  # variabilidade entre sujeitos, Vd
  ruv_prop    = 0.199,  # erro residual proporcional
  ruv_adit    = 2.4,    # erro residual aditivo, mg/L
  loq         = 0.6     # limite de quantificação do ensaio, mg/L
)

# Distribuição de covariáveis, conforme a coorte original (n = 206, Tabela 1).
COORTE <- list(
  peso_media = 74.8, peso_dp = 15.8,
  crcl_media = 90.7, crcl_dp = 60.4
)

#' Sorteio lognormal com média e coeficiente de variação especificados.
#' Usado para variabilidade entre sujeitos, que Roberts modela como exponencial.
rlnorm_cv <- function(n, media, cv) {
  s <- sqrt(log(1 + cv^2))
  media * exp(rnorm(n, 0, s) - s^2 / 2)
}

# =============================================================================
# SEGUNDO MECANISMO GERADOR: Llopis-Salvia & Jimenez-Torres 2006
#     J Clin Pharm Ther 31:447-454, Tabela 3.
#
# Incluído para responder à limitação do gerador único. As três razões da
# escolha:
#
#   1. Deriva de administração INTERMITENTE (infusão de no mínimo 60 min, doses
#      duas vezes ao dia), que é a modalidade contemplada pela ferramenta. O
#      modelo de Roberts deriva de infusão contínua.
#   2. Reporta os três componentes exigidos: efeitos fixos, variabilidade entre
#      sujeitos para os três volumes, e erro residual COMBINADO (proporcional
#      mais aditivo), na mesma forma adotada por Roberts.
#   3. População de pacientes críticos, comparável à de Roberts, o que mantém a
#      comparação entre geradores centrada na modalidade de administração e não
#      no perfil dos pacientes.
#
# Estrutura de DOIS compartimentos. Como o estimador da ferramenta opera em um
# compartimento no método de ponto único, este gerador introduz também má
# especificação ESTRUTURAL, e não apenas de parâmetro. Isso é deliberado: é a
# condição que um paciente real apresenta.
#
# NOTA SOBRE O VOLUME: Vss = (0,414 + 1,32) × peso ≈ 1,73 L/kg, portanto ainda
# MAIOR que o 1,53 L/kg de Roberts. Dois modelos independentes, um de infusão
# contínua e outro de intermitente, ambos de população crítica, convergem em
# volume elevado. O volume alto parece característico de pacientes críticos, e
# não artefato da modalidade de administração.
# =============================================================================

LLOPIS <- list(
  # CL (L/h) = 0,034 × ClCr(mL/min) + 0,015 × TBW(kg)
  cl_theta1   = 0.034,
  cl_theta2   = 0.015,
  vc_kg       = 0.414,  # L/kg
  q_lh        = 7.48,   # L/h, sem covariável
  vp_kg       = 1.32,   # L/kg
  cv_cl       = 0.292,  # variabilidade entre sujeitos, CL
  cv_vc       = 0.364,  # variabilidade entre sujeitos, Vc
  cv_vp       = 0.398,  # variabilidade entre sujeitos, Vp
  ruv_prop    = 0.239,  # erro residual proporcional (sigma1)
  ruv_adit    = 1.85,   # erro residual aditivo, mg/L (sigma2 = 18,5%)
  loq         = 0.6
)

# Coorte do Llopis-Salvia (grupo A, n = 50; Tabela 1).
COORTE_LLOPIS <- list(
  peso_media = 71.2, peso_dp = 12.3,
  crcl_media = 69.9, crcl_dp = 24.9,
  crcl_max   = 120    # o modelo trunca o ClCr em 120 mL/min
)

#' Concentração no steady-state para modelo de DOIS compartimentos, infusão
#' intermitente, medida `t` horas após o fim da infusão.
#'
#' Usada apenas na GERAÇÃO dos dados. O estimador continua operando em um
#' compartimento, o que é justamente a má especificação estrutural que este
#' gerador introduz.
conc_ss_2comp <- function(cl, vc, q, vp, dose, tinf, tau, t) {
  k10 <- cl / vc; k12 <- q / vc; k21 <- q / vp
  soma <- k10 + k12 + k21
  disc <- soma^2 - 4 * k10 * k21
  if (disc < 0) return(NA_real_)
  d <- sqrt(disc)
  alpha <- (soma + d) / 2
  beta  <- (soma - d) / 2
  if (alpha <= 0 || beta <= 0) return(NA_real_)
  taxa <- dose / tinf
  A <- (k21 - alpha) / ((beta - alpha) * vc)
  B <- (k21 - beta)  / ((alpha - beta) * vc)
  termo <- function(coef, lam)
    coef / lam * ((1 - exp(-lam * tinf)) * exp(-lam * t)) / (1 - exp(-lam * tau))
  taxa * termo(A, alpha) + taxa * termo(B, beta)
}

#' Sorteia um paciente virtual pelo modelo de Llopis-Salvia.
#' Devolve também a AUC24 verdadeira, que no modelo de 2 compartimentos
#' continua sendo a dose diária dividida pelo clearance.
paciente_llopis <- function(dose, tau) {
  peso <- max(40, rnorm(1, COORTE_LLOPIS$peso_media, COORTE_LLOPIS$peso_dp))
  crcl <- min(COORTE_LLOPIS$crcl_max,
              max(15, rnorm(1, COORTE_LLOPIS$crcl_media, COORTE_LLOPIS$crcl_dp)))
  cl_tipico <- LLOPIS$cl_theta1 * crcl + LLOPIS$cl_theta2 * peso
  list(
    peso = peso, crcl = crcl,
    cl = rlnorm_cv(1, cl_tipico, LLOPIS$cv_cl),
    vc = rlnorm_cv(1, LLOPIS$vc_kg * peso, LLOPIS$cv_vc),
    q  = LLOPIS$q_lh,
    vp = rlnorm_cv(1, LLOPIS$vp_kg * peso, LLOPIS$cv_vp)
  )
}

# =============================================================================
# MODELO ESTRUTURAL — 1 compartimento, infusão intermitente, steady-state
# =============================================================================

#' Concentração no steady-state, `t` horas APÓS o fim da infusão.
#'
#' C(t) = (Dose/(T·CL)) · (1 − e^(−k·T)) · e^(−k·t) / (1 − e^(−k·τ))
#'
#' Note a divisão por T: a fórmula parte da TAXA de infusão, não da dose.
#' As duas formas coincidem quando T = 1 h, o que já mascarou um defeito de
#' implementação neste projeto.
conc_ss <- function(cl, vd, dose, tinf, tau, t) {
  k <- cl / vd
  den <- 1 - exp(-k * tau)
  ifelse(den > 0,
         (dose / (tinf * cl)) * (1 - exp(-k * tinf)) * exp(-k * t) / den,
         NA_real_)
}

# =============================================================================
# M — MÉTODOS AVALIADOS
# =============================================================================

#' Estimativa MAP do clearance a partir de UMA concentração.
#'
#' Minimiza  obj(CL) = ((C_pred − C_obs)/σ_obs)² + ((CL − CL_prior)/ω_CL)²
#' com o volume de distribuição FIXO no valor populacional.
#'
#' A decisão de não estimar o volume decorre de identificabilidade: uma única
#' concentração é compatível com múltiplas combinações de CL e Vd, e deixar
#' ambos variarem consome informação que deveria recair sobre o clearance — do
#' qual a AUC depende diretamente. A Tabela 8 do manuscrito testa essa decisão.
map_1ponto <- function(c_obs, t_obs, dose, tinf, tau,
                       cl_prior, cl_prior_dp, vd_pop, cv_obs = 0.15) {
  sigma_obs <- max(cv_obs * c_obs, 0.5)

  obj <- function(cl) {
    cp <- conc_ss(cl, vd_pop, dose, tinf, tau, t_obs)
    if (!is.finite(cp) || cp <= 0) return(1e12)
    ((cp - c_obs) / sigma_obs)^2 + ((cl - cl_prior) / cl_prior_dp)^2
  }

  # optimize() é determinístico: dada a mesma entrada, devolve o mesmo ótimo.
  # Não introduz aleatoriedade adicional à do sorteio dos pacientes.
  op <- optimize(obj, interval = c(cl_prior * 0.15, cl_prior * 3.5), tol = 1e-6)
  cl_est <- op$minimum

  # Incerteza pela curvatura: obj = −2·log(posterior), logo Cov = 2·H⁻¹.
  # O fator 2 é a armadilha clássica desta derivação; omiti-lo estreita o
  # intervalo por um fator de 1,41.
  h <- cl_est * 1e-3
  d2 <- (obj(cl_est + h) - 2 * obj(cl_est) + obj(cl_est - h)) / h^2
  dp_cond <- if (is.finite(d2) && d2 > 0) sqrt(2 / d2) else NA_real_

  list(cl = cl_est, dp = dp_cond, obj_fn = obj)
}

#' Estimativa MAP permitindo que o volume também varie. Só para a Tabela 8.
map_1ponto_vd_livre <- function(c_obs, t_obs, dose, tinf, tau,
                                cl_prior, cl_prior_dp, vd_prior, vd_prior_dp,
                                cv_obs = 0.15) {
  sigma_obs <- max(cv_obs * c_obs, 0.5)
  obj <- function(par) {
    cl <- par[1]; vd <- par[2]
    if (cl <= 0 || vd <= 0) return(1e12)
    cp <- conc_ss(cl, vd, dose, tinf, tau, t_obs)
    if (!is.finite(cp) || cp <= 0) return(1e12)
    ((cp - c_obs) / sigma_obs)^2 +
      ((cl - cl_prior) / cl_prior_dp)^2 +
      ((vd - vd_prior) / vd_prior_dp)^2
  }
  op <- optim(c(cl_prior, vd_prior), obj, method = "Nelder-Mead",
              control = list(reltol = 1e-8))
  op$par[1]
}

#' Propaga a incerteza do volume POPULACIONAL para o clearance, pelo método delta.
#'
#' Sem esta correção a cobertura do IC95 fica muito abaixo do nominal: o volume
#' entra na estimação como se fosse conhecido, quando é populacional.
propagar_vd <- function(c_obs, t_obs, dose, tinf, tau,
                        cl_prior, cl_prior_dp, vd_pop, vd_pop_dp, cv_obs = 0.15) {
  base <- map_1ponto(c_obs, t_obs, dose, tinf, tau, cl_prior, cl_prior_dp, vd_pop, cv_obs)
  if (is.na(base$dp)) return(base)
  hv <- vd_pop * 0.05
  cl_mais  <- map_1ponto(c_obs, t_obs, dose, tinf, tau, cl_prior, cl_prior_dp, vd_pop + hv, cv_obs)$cl
  cl_menos <- map_1ponto(c_obs, t_obs, dose, tinf, tau, cl_prior, cl_prior_dp, vd_pop - hv, cv_obs)$cl
  dcl_dvd <- (cl_mais - cl_menos) / (2 * hv)
  base$dp <- sqrt(base$dp^2 + (dcl_dvd^2) * (vd_pop_dp^2))
  base
}

#' Estimativa com DOIS pontos: Sawchuk-Zaske seguido de atualização conjugada.
map_2pontos <- function(c_pico, t_pico, c_vale, t_vale, dose, tinf, tau,
                        cl_prior, cl_prior_dp, cv_obs = 0.15) {
  k <- log(c_pico / c_vale) / (t_vale - t_pico)
  if (!is.finite(k) || k <= 0) return(NA_real_)
  c_pico_corr <- c_pico * exp(k * t_pico)
  vd <- ((dose / tinf) * (1 - exp(-k * tinf))) /
        (k * (1 - exp(-k * tau)) * c_pico_corr)
  cl_obs <- k * vd
  if (!is.finite(cl_obs) || cl_obs <= 0) return(NA_real_)
  s_obs <- cv_obs * cl_obs
  (cl_prior / cl_prior_dp^2 + cl_obs / s_obs^2) / (1 / cl_prior_dp^2 + 1 / s_obs^2)
}

# =============================================================================
# PRIOR DE ESTIMAÇÃO — modelo estratificado por ClCr do VancoOptim
# É o prior que a ferramenta aplica; NÃO é o gerador.
# =============================================================================

prior_estratificado <- function(crcl, peso) {
  if (crcl >= 90)      list(cl = 5.5, vd_kg = 0.50)
  else if (crcl >= 60) list(cl = 4.0, vd_kg = 0.55)
  else if (crcl >= 30) list(cl = 2.5, vd_kg = 0.60)
  else if (crcl >= 15) list(cl = 1.2, vd_kg = 0.65)
  else                 list(cl = 0.6, vd_kg = 0.65)
}

# =============================================================================
# P — MEDIDAS DE DESEMPENHO
# =============================================================================

desempenho <- function(erros, dentro_ic, acertos) {
  list(
    n       = length(erros),
    mpe     = mean(erros),
    mdpe    = median(erros),
    mape    = mean(abs(erros)),
    mdape   = median(abs(erros)),
    # Erro de Monte Carlo: quanto a própria simulação é imprecisa.
    # Reportá-lo é o que distingue estudo de simulação de teste ad hoc.
    mc_mpe  = sd(erros) / sqrt(length(erros)),
    mc_mape = sd(abs(erros)) / sqrt(length(erros)),
    cobertura = 100 * mean(dentro_ic, na.rm = TRUE),
    acerto    = 100 * mean(acertos, na.rm = TRUE)
  )
}

faixa_auc <- function(auc) ifelse(auc < 400, 0L, ifelse(auc <= 600, 1L, 2L))

# =============================================================================
# CENÁRIO: gera N pacientes e avalia os métodos
# =============================================================================

#' Cenário com o gerador de Llopis-Salvia (2 compartimentos, intermitente).
#'
#' Diferença central em relação a `rodar_cenario()`: a concentração é gerada por
#' modelo de DOIS compartimentos, enquanto o estimador continua operando em um.
#' Introduz, portanto, má especificação estrutural além da de parâmetro.
rodar_cenario_llopis <- function(t_coleta, n = N_REPLICAS,
                                 dose = 1000, tinf = 2, tau = 12,
                                 metodo = "1ponto") {

  erros <- numeric(0); dentro <- logical(0); acertos <- logical(0)

  for (i in seq_len(n)) {
    pv <- paciente_llopis(dose, tau)
    auc_verdadeira <- dose * (24 / tau) / pv$cl

    c_real <- conc_ss_2comp(pv$cl, pv$vc, pv$q, pv$vp, dose, tinf, tau, t_coleta)
    if (!is.finite(c_real) || c_real <= 0) next
    c_obs <- c_real * (1 + LLOPIS$ruv_prop * rnorm(1)) + LLOPIS$ruv_adit * rnorm(1)
    if (c_obs < LLOPIS$loq) next

    pr <- prior_estratificado(pv$crcl, pv$peso)
    cl_prior_dp <- pr$cl * 0.30
    vd_pop <- pr$vd_kg * pv$peso
    vd_pop_dp <- 0.15 * pv$peso

    if (metodo == "1ponto") {
      est <- propagar_vd(c_obs, t_coleta, dose, tinf, tau,
                         pr$cl, cl_prior_dp, vd_pop, vd_pop_dp)
      cl_est <- est$cl
      auc_est <- dose * (24 / tau) / cl_est
      if (!is.na(est$dp)) {
        lo <- dose * (24 / tau) / (cl_est + 1.96 * est$dp)
        hi <- dose * (24 / tau) / max(1e-6, cl_est - 1.96 * est$dp)
        dentro <- c(dentro, auc_verdadeira >= lo && auc_verdadeira <= hi)
      } else dentro <- c(dentro, NA)
    } else {   # dois pontos: pico 1 h após o fim da infusão
      t_pico <- 1
      cp_real <- conc_ss_2comp(pv$cl, pv$vc, pv$q, pv$vp, dose, tinf, tau, t_pico)
      cp <- cp_real * (1 + LLOPIS$ruv_prop * rnorm(1)) + LLOPIS$ruv_adit * rnorm(1)
      if (!is.finite(cp) || cp <= c_obs) next
      cl_est <- map_2pontos(cp, t_pico, c_obs, t_coleta, dose, tinf, tau,
                            pr$cl, cl_prior_dp)
      if (is.na(cl_est) || cl_est <= 0) next
      auc_est <- dose * (24 / tau) / cl_est
      dentro <- c(dentro, NA)
    }

    if (!is.finite(auc_est)) next
    erros <- c(erros, 100 * (auc_est - auc_verdadeira) / auc_verdadeira)
    acertos <- c(acertos, faixa_auc(auc_est) == faixa_auc(auc_verdadeira))
  }

  desempenho(erros, dentro, acertos)
}

rodar_cenario <- function(vd_kg_gerador, t_coleta, n = N_REPLICAS,
                          dose = 1000, tinf = 2, tau = 12, metodo = "1ponto") {

  erros <- numeric(0); dentro <- logical(0); acertos <- logical(0)

  for (i in seq_len(n)) {
    # ── paciente virtual
    peso <- max(35, rnorm(1, COORTE$peso_media, COORTE$peso_dp))
    crcl <- max(20, rnorm(1, COORTE$crcl_media, COORTE$crcl_dp))
    cl_v <- rlnorm_cv(1, ROBERTS$cl_tipico * (crcl / 100), ROBERTS$cv_cl)
    vd_v <- rlnorm_cv(1, vd_kg_gerador * peso, ROBERTS$cv_vd)

    auc_verdadeira <- dose * (24 / tau) / cl_v

    # ── concentração observada, com erro residual combinado
    c_real <- conc_ss(cl_v, vd_v, dose, tinf, tau, t_coleta)
    if (!is.finite(c_real)) next
    c_obs <- c_real * (1 + ROBERTS$ruv_prop * rnorm(1)) + ROBERTS$ruv_adit * rnorm(1)
    if (c_obs < ROBERTS$loq) next   # abaixo do limite de quantificação

    pr <- prior_estratificado(crcl, peso)
    cl_prior_dp <- pr$cl * 0.30
    vd_pop <- pr$vd_kg * peso
    vd_pop_dp <- 0.15 * peso

    if (metodo == "1ponto") {
      est <- propagar_vd(c_obs, t_coleta, dose, tinf, tau,
                         pr$cl, cl_prior_dp, vd_pop, vd_pop_dp)
      cl_est <- est$cl
      auc_est <- dose * (24 / tau) / cl_est
      if (!is.na(est$dp)) {
        lo <- dose * (24 / tau) / (cl_est + 1.96 * est$dp)
        hi <- dose * (24 / tau) / max(1e-6, cl_est - 1.96 * est$dp)
        dentro <- c(dentro, auc_verdadeira >= lo && auc_verdadeira <= hi)
      } else dentro <- c(dentro, NA)

    } else if (metodo == "1ponto_vd_livre") {
      cl_est <- map_1ponto_vd_livre(c_obs, t_coleta, dose, tinf, tau,
                                    pr$cl, cl_prior_dp, vd_pop, vd_pop_dp)
      auc_est <- dose * (24 / tau) / cl_est
      dentro <- c(dentro, NA)

    } else {  # dois pontos: pico 1 h após o fim da infusão
      t_pico <- 1
      cp_real <- conc_ss(cl_v, vd_v, dose, tinf, tau, t_pico)
      cp <- cp_real * (1 + ROBERTS$ruv_prop * rnorm(1)) + ROBERTS$ruv_adit * rnorm(1)
      if (!is.finite(cp) || cp <= c_obs) next
      cl_est <- map_2pontos(cp, t_pico, c_obs, t_coleta, dose, tinf, tau,
                            pr$cl, cl_prior_dp)
      if (is.na(cl_est) || cl_est <= 0) next
      auc_est <- dose * (24 / tau) / cl_est
      dentro <- c(dentro, NA)
    }

    if (!is.finite(auc_est)) next
    erros <- c(erros, 100 * (auc_est - auc_verdadeira) / auc_verdadeira)
    acertos <- c(acertos, faixa_auc(auc_est) == faixa_auc(auc_verdadeira))
  }

  desempenho(erros, dentro, acertos)
}

# =============================================================================
# EXECUÇÃO
# =============================================================================

fmt <- function(x, d = 1) formatC(x, format = "f", digits = d)

rodar_tudo <- function() {

  # A semente é reposicionada no início. Toda a sequência de sorteios que segue
  # é determinada por ela — reexecutar produz exatamente os mesmos números.
  set.seed(SEMENTE)

  cat("=======================================================================\n")
  cat("VancoOptim — estudo de simulação\n")
  cat("Gerador: Roberts et al. 2011 (externo à biblioteca de priors)\n")
  cat("Semente:", SEMENTE, "| Réplicas por cenário:", N_REPLICAS, "\n")
  cat("R:", R.version.string, "\n")
  cat("=======================================================================\n\n")

  # ── Tabela 7: desempenho do método com concentração única
  cat("TABELA 7 — estimativa com concentração única\n")
  cat(sprintf("%-12s %-8s %8s %8s %8s %8s %10s %8s %12s\n",
              "Vd gerador", "Coleta", "MPE", "MdPE", "MAPE", "MdAPE",
              "Cobertura", "Acerto", "MC(MAPE)"))
  cat(strrep("-", 92), "\n")

  t7 <- list()
  for (vd in c(1.53, 1.10, 0.90, 0.70, 0.55)) {
    for (tc in c(10, 5)) {
      r <- rodar_cenario(vd, tc)
      t7[[paste(vd, tc)]] <- r
      cat(sprintf("%-12s %-8s %7s%% %7s%% %7s%% %7s%% %9s%% %7s%% %11s\n",
                  paste0(fmt(vd, 2), " L/kg"), paste0(tc, " h"),
                  fmt(r$mpe), fmt(r$mdpe), fmt(r$mape), fmt(r$mdape),
                  fmt(r$cobertura), fmt(r$acerto), paste0("±", fmt(r$mc_mape, 2))))
    }
  }

  # ── Tabela 8: volume fixo versus estimado
  cat("\n\nTABELA 8 — volume populacional fixo versus estimado (MAPE)\n")
  cat(sprintf("%-12s %-8s %12s %12s\n", "Vd gerador", "Coleta", "Vd fixo", "Vd estimado"))
  cat(strrep("-", 48), "\n")

  t8 <- list()
  for (vd in c(1.53, 1.10, 0.90, 0.70, 0.55)) {
    for (tc in c(10, 5)) {
      set.seed(SEMENTE + 1000L)   # mesma amostra para os dois métodos
      a <- rodar_cenario(vd, tc, metodo = "1ponto")
      set.seed(SEMENTE + 1000L)
      b <- rodar_cenario(vd, tc, metodo = "1ponto_vd_livre")
      t8[[paste(vd, tc)]] <- list(fixo = a$mape, livre = b$mape)
      cat(sprintf("%-12s %-8s %11s%% %11s%%\n",
                  paste0(fmt(vd, 2), " L/kg"), paste0(tc, " h"),
                  fmt(a$mape), fmt(b$mape)))
    }
  }

  # ── Comparação com o método de dois pontos
  cat("\n\nCOMPARAÇÃO — um ponto versus dois pontos (MAPE)\n")
  cat(sprintf("%-12s %-8s %12s %12s\n", "Vd gerador", "Coleta", "1 ponto", "2 pontos"))
  cat(strrep("-", 48), "\n")

  cmp <- list()
  for (vd in c(1.53, 0.90, 0.55)) {
    set.seed(SEMENTE + 2000L)
    a <- rodar_cenario(vd, 10, metodo = "1ponto")
    set.seed(SEMENTE + 2000L)
    b <- rodar_cenario(vd, 10, metodo = "2pontos")
    cmp[[as.character(vd)]] <- list(um = a$mape, dois = b$mape)
    cat(sprintf("%-12s %-8s %11s%% %11s%%\n",
                paste0(fmt(vd, 2), " L/kg"), "10 h", fmt(a$mape), fmt(b$mape)))
  }

  # ── Tabela 12: segundo gerador, derivado de administração intermitente
  cat("\n\nTABELA 12 — gerador de Llopis-Salvia (intermitente, 2 compartimentos)\n")
  cat(sprintf("%-10s %-10s %8s %8s %8s %8s %10s %12s\n",
              "Método", "Coleta", "MPE", "MdPE", "MAPE", "MdAPE", "Cobertura", "MC(MAPE)"))
  cat(strrep("-", 80), "\n")

  t12 <- list()
  for (met in c("1ponto", "2pontos")) {
    for (tc in c(10, 5)) {
      set.seed(SEMENTE + 4000L)
      r <- rodar_cenario_llopis(tc, metodo = met)
      t12[[paste(met, tc)]] <- r
      cat(sprintf("%-10s %-10s %7s%% %7s%% %7s%% %7s%% %9s %11s\n",
                  ifelse(met == "1ponto", "1 ponto", "2 pontos"), paste0(tc, " h"),
                  fmt(r$mpe), fmt(r$mdpe), fmt(r$mape), fmt(r$mdape),
                  ifelse(is.nan(r$cobertura), "  —  ", paste0(fmt(r$cobertura), "%")),
                  paste0("±", fmt(r$mc_mape, 2))))
    }
  }

  cat("\n\nNOTA: MC(MAPE) é o erro-padrão de Monte Carlo — a imprecisão da\n")
  cat("própria simulação. Valores abaixo de 1 ponto percentual indicam que\n")
  cat("o número de réplicas é suficiente para as comparações pretendidas.\n")

  invisible(list(tabela7 = t7, tabela8 = t8, comparacao = cmp, tabela12 = t12,
                 semente = SEMENTE, replicas = N_REPLICAS,
                 sessao = sessionInfo()))
}

# =============================================================================
# CALIBRAÇÃO DA JANELA DE COLETA — população local
#
# Roda com a distribuição de meia-vida OBSERVADA no Hospital Moriah, e não com
# a do gerador de Roberts. É o que sustenta o achado de que a recomendação de
# janela depende do perfil populacional.
#
# [A DEFINIR] Substituir por dados atualizados quando a coorte crescer.
# =============================================================================

K_OBSERVADOS_MORIAH <- c(0.1638, 0.0783, 0.0654, 0.1770, 0.1251,
                         0.1220, 0.1720, 0.0600, 0.0809, 0.0392, 0.1191)

calibrar_janela <- function(tau = 12, tinf = 2, n = N_REPLICAS, dose = 1000) {
  set.seed(SEMENTE + 3000L)

  t12 <- sort(log(2) / K_OBSERVADOS_MORIAH)
  cat("\n\nCALIBRAÇÃO DA JANELA DE COLETA — população local\n")
  cat("Meia-vida observada (n =", length(t12), "): mediana",
      fmt(median(t12)), "h | quartis", fmt(quantile(t12, 0.25)), "-",
      fmt(quantile(t12, 0.75)), "h\n\n")
  cat(sprintf("Regime: %d mg q%dh, infusão de %d h\n", dose, tau, tinf))
  cat(sprintf("%-22s %12s %12s\n", "h após fim da infusão", "MAPE", "MC(MAPE)"))
  cat(strrep("-", 50), "\n")

  saida <- list()
  for (tc in seq_len(tau - tinf)) {
    erros <- numeric(0)
    for (i in seq_len(n)) {
      peso <- max(35, rnorm(1, COORTE$peso_media, COORTE$peso_dp))
      k_base <- sample(K_OBSERVADOS_MORIAH, 1)
      k_v <- rlnorm_cv(1, k_base, 0.20)
      vd_v <- rlnorm_cv(1, 0.60 * peso, 0.30)
      cl_v <- k_v * vd_v
      auc_v <- dose * (24 / tau) / cl_v

      c_real <- conc_ss(cl_v, vd_v, dose, tinf, tau, tc)
      if (!is.finite(c_real)) next
      c_obs <- c_real * (1 + 0.08 * rnorm(1))
      if (c_obs < ROBERTS$loq) next

      pr <- prior_estratificado(90, peso)   # adulto típico
      est <- propagar_vd(c_obs, tc, dose, tinf, tau,
                         pr$cl, pr$cl * 0.30, pr$vd_kg * peso, 0.15 * peso)
      auc_e <- dose * (24 / tau) / est$cl
      if (is.finite(auc_e)) erros <- c(erros, 100 * (auc_e - auc_v) / auc_v)
    }
    mape <- mean(abs(erros)); mc <- sd(abs(erros)) / sqrt(length(erros))
    saida[[as.character(tc)]] <- list(mape = mape, mc = mc)
    cat(sprintf("%-22s %11s%% %11s\n", paste0(tc, " h"),
                fmt(mape, 2), paste0("±", fmt(mc, 2))))
  }
  invisible(saida)
}

# =============================================================================
# Para executar:
#
#   res  <- rodar_tudo()
#   cal  <- calibrar_janela()
#
# Registrar no manuscrito: a semente, o número de réplicas e a saída de
# sessionInfo(), que identifica a versão de R usada.
#
#   res$sessao
# =============================================================================