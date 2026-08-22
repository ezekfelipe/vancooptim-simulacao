# =============================================================================
# Goti 2018 — implementação de referência em R
#
# Comparador independente para a estimativa bayesiana bicompartimental do
# VancoOptim. Escrito de forma explícita, sem otimizadores de caixa-preta, para
# que cada passo do cálculo seja auditável a partir do código.
#
# MOTIVAÇÃO
# A comparação contra o TDMx é validação externa, mas depende de uma ferramenta
# de terceiro permanecer disponível e inalterada. Esta implementação cumpre
# papel complementar: pode ser depositada junto do estudo, e quem quiser
# replicar não depende de nenhum serviço externo.
#
# Ela também serve como árbitro. Quando VancoOptim e TDMx divergem, uma terceira
# implementação independente indica qual das duas se afasta — foi assim que se
# identificou, em 08/2026, um erro de conversão de tempo no ramo bicompartimental
# do VancoOptim.
#
# LIMITE DESTA VERIFICAÇÃO
# Os parâmetros populacionais abaixo foram transcritos da implementação do
# VancoOptim, e NÃO conferidos contra a publicação original de Goti e
# colaboradores. Uma implementação de referência que herda os parâmetros da
# ferramenta avaliada verifica a ARITMÉTICA, não a fidelidade à fonte. Conferir
# os valores contra o artigo original é passo pendente e está sinalizado abaixo.
#
# Uso:   source("goti_referencia.R")
#        exemplo()
# =============================================================================

# ── PARÂMETROS POPULACIONAIS ────────────────────────────────────────────────
# [A CONFERIR contra Goti V, Chaturvedula A, Fossler MJ, Mok S, Jacob JT.
#  Ther Drug Monit. 2018;40(2):212-221.]

GOTI <- list(
  cl_tipico   = 4.5,    # L/h, para ClCr de 120 mL/min
  cl_expoente = 0.8,    # expoente do ClCr
  cl_ref_crcl = 120,    # mL/min
  vc_70kg     = 58.4,   # L, para 70 kg
  vp          = 38.4,   # L
  q           = 6.5,    # L/h
  omega_cl    = 0.398,  # CV da variabilidade entre sujeitos, CL
  omega_vc    = 0.816,  # CV da variabilidade entre sujeitos, Vc
  sigma_obs   = 3.4,    # desvio-padrão da observação, mg/L (aditivo)
  # Covariáveis de hemodiálise
  fator_cl_hd = 0.7,
  fator_vc_hd = 0.5
)

#' Constantes microscópicas do modelo de dois compartimentos.
#'
#' A partir de CL, Vc, Q e Vp, obtém k21 e os autovalores alpha e beta que
#' descrevem as fases de distribuição e de eliminação terminal.
micro_2comp <- function(cl, vc, q, vp) {
  k10 <- cl / vc
  k12 <- q  / vc
  k21 <- q  / vp
  soma <- k10 + k12 + k21
  disc <- soma^2 - 4 * k10 * k21
  if (disc < 0) return(NULL)
  d <- sqrt(disc)
  list(k10 = k10, k12 = k12, k21 = k21,
       alpha = (soma + d) / 2,
       beta  = (soma - d) / 2)
}

#' Concentração no steady-state, dois compartimentos, infusão intermitente.
#'
#' CONVENÇÃO DE TEMPO: `t` é contado a partir do INÍCIO da infusão.
#' A função trata separadamente o período durante a infusão (t <= tinf) e o
#' período posterior.
#'
#' Esta convenção é a fonte de um erro identificado na implementação do
#' VancoOptim: a chamada convertia o tempo do pico somando a duração da infusão,
#' mas passava o tempo do vale sem a mesma conversão. Como ambos são medidos a
#' partir do FIM da infusão na interface, o vale era tratado como colhido mais
#' cedo do que realmente foi.
conc_ss_2comp <- function(t, dose, tinf, tau, cl, vc, q, vp) {
  m <- micro_2comp(cl, vc, q, vp)
  if (is.null(m) || m$alpha <= 0 || m$beta <= 0) return(NA_real_)

  taxa <- dose / tinf
  A <- (m$k21 - m$alpha) / ((m$beta  - m$alpha) * vc)
  B <- (m$k21 - m$beta)  / ((m$alpha - m$beta)  * vc)

  termo <- function(coef, lam) {
    if (t <= tinf) {
      # Durante a infusão: acúmulo da dose corrente somado ao resíduo das anteriores.
      coef / lam * ((1 - exp(-lam * t)) +
        exp(-lam * tau) * (1 - exp(-lam * tinf)) * exp(-lam * (t - tinf)) /
        (1 - exp(-lam * tau)))
    } else {
      # Após a infusão: decaimento biexponencial.
      coef / lam * ((1 - exp(-lam * tinf)) * exp(-lam * (t - tinf))) /
        (1 - exp(-lam * tau))
    }
  }
  taxa * termo(A, m$alpha) + taxa * termo(B, m$beta)
}

#' Prior populacional do Goti para um paciente.
#'
#' O clearance de creatinina segue Cockcroft-Gault com peso total, com fator de
#' sexo e teto de 150 mL/min, como na implementação avaliada.
prior_goti <- function(idade, peso, scr, sexo = "M", dialise = FALSE) {
  scr_ajust <- if (scr < 1.0 && idade > 65) 1.0 else scr
  crcl <- ((140 - idade) * peso) / (72 * scr_ajust)
  if (sexo == "F") crcl <- crcl * 0.85
  crcl <- min(crcl, 150)

  f_cl <- if (dialise) GOTI$fator_cl_hd else 1
  f_vc <- if (dialise) GOTI$fator_vc_hd else 1

  list(
    crcl = crcl,
    cl = GOTI$cl_tipico * (crcl / GOTI$cl_ref_crcl)^GOTI$cl_expoente * f_cl,
    vc = GOTI$vc_70kg * (peso / 70) * f_vc,
    vp = GOTI$vp,
    q  = GOTI$q
  )
}

#' Estimativa MAP de CL e Vc, com Q e Vp fixos no valor populacional.
#'
#' Minimiza  obj = Σ((pred − obs)/σ_obs)² + ((CL − CL_prior)/ω_CL)²
#'                                        + ((Vc − Vc_prior)/ω_Vc)²
#'
#' que corresponde a −2·log(posterior) para verossimilhança e prior gaussianos.
#'
#' Q e Vp permanecem fixos porque duas concentrações não os identificam. A
#' consequência é que a partição entre compartimento central e periférico é
#' determinada em parte pelo prior, e não pelos dados — ponto relevante ao
#' interpretar divergências entre implementações.
#'
#' @param t_pico,t_vale tempos contados do INÍCIO da infusão.
map_goti <- function(pico, t_pico, vale, t_vale, dose, tinf, tau,
                     pr, sigma_obs = GOTI$sigma_obs) {

  sd_cl <- pr$cl * GOTI$omega_cl
  sd_vc <- pr$vc * GOTI$omega_vc

  obj <- function(par) {
    cl <- par[1]; vc <- par[2]
    if (cl <= 0 || vc <= 0) return(1e12)
    cp <- conc_ss_2comp(t_pico, dose, tinf, tau, cl, vc, pr$q, pr$vp)
    cv <- conc_ss_2comp(t_vale, dose, tinf, tau, cl, vc, pr$q, pr$vp)
    if (!is.finite(cp) || !is.finite(cv)) return(1e12)
    ((cp - pico) / sigma_obs)^2 + ((cv - vale) / sigma_obs)^2 +
      ((cl - pr$cl) / sd_cl)^2 + ((vc - pr$vc) / sd_vc)^2
  }

  # Otimização em duas etapas. A busca em malha localiza a bacia do mínimo, e o
  # Nelder-Mead refina. A malha isolada, empregada na implementação avaliada,
  # limita a resolução ao passo final e afeta a hessiana numérica calculada em
  # seguida.
  melhor <- c(pr$cl, pr$vc); v_melhor <- obj(melhor)
  for (cl in seq(pr$cl * 0.3, pr$cl * 2.5, length.out = 60)) {
    for (vc in seq(pr$vc * 0.3, pr$vc * 2.5, length.out = 60)) {
      v <- obj(c(cl, vc))
      if (v < v_melhor) { v_melhor <- v; melhor <- c(cl, vc) }
    }
  }
  op <- optim(melhor, obj, method = "Nelder-Mead",
              control = list(reltol = 1e-10, maxit = 2000))

  cl_est <- op$par[1]; vc_est <- op$par[2]

  # Covariância posterior pela curvatura. Como obj = −2·log(posterior), a
  # covariância aproximada é 2·H⁻¹. Omitir o fator 2 estreita o intervalo por
  # um fator de 1,41.
  h <- c(cl_est, vc_est) * 1e-4
  hess <- matrix(0, 2, 2)
  for (i in 1:2) for (j in 1:2) {
    ei <- rep(0, 2); ei[i] <- h[i]
    ej <- rep(0, 2); ej[j] <- h[j]
    hess[i, j] <- (obj(op$par + ei + ej) - obj(op$par + ei - ej) -
                   obj(op$par - ei + ej) + obj(op$par - ei - ej)) / (4 * h[i] * h[j])
  }
  cov_post <- tryCatch(2 * solve(hess), error = function(e) NULL)
  se <- if (!is.null(cov_post) && all(diag(cov_post) > 0)) sqrt(diag(cov_post)) else c(NA, NA)

  m <- micro_2comp(cl_est, vc_est, pr$q, pr$vp)
  list(
    cl = cl_est, vc = vc_est, vss = vc_est + pr$vp,
    se_cl = se[1], se_vc = se[2],
    t12_alfa = log(2) / m$alpha,
    t12_beta = log(2) / m$beta,
    obj = op$value, convergiu = op$convergence == 0
  )
}

#' AUC24 no steady-state. Depende apenas do clearance e da dose diária.
auc24 <- function(dose, tau, cl) dose * (24 / tau) / cl

# =============================================================================
# EXEMPLO: caso 3 da comparação com o TDMx
#
# Concentrações geradas por modelo de dois compartimentos com parâmetros
# conhecidos, para que a divergência entre implementações não se confunda com
# ruído de medição.
# =============================================================================

exemplo <- function() {
  # Parâmetros verdadeiros do paciente simulado
  verdadeiro <- list(cl = 3.0, vc = 40.0, q = 6.5, vp = 30.0)
  dose <- 750; tinf <- 2; tau <- 12

  # Concentrações geradas, tempos contados do INÍCIO da infusão
  t_pico <- 1.0 + tinf     # 1 h após o fim
  t_vale <- 9.5 + tinf     # 9,5 h após o fim
  pico <- conc_ss_2comp(t_pico, dose, tinf, tau,
                        verdadeiro$cl, verdadeiro$vc, verdadeiro$q, verdadeiro$vp)
  vale <- conc_ss_2comp(t_vale, dose, tinf, tau,
                        verdadeiro$cl, verdadeiro$vc, verdadeiro$q, verdadeiro$vp)

  pr <- prior_goti(idade = 55, peso = 70, scr = 0.9, sexo = "M")
  est <- map_goti(pico, t_pico, vale, t_vale, dose, tinf, tau, pr)

  cat("=======================================================================\n")
  cat("Goti 2018 — implementação de referência em R\n")
  cat("=======================================================================\n\n")
  cat(sprintf("Paciente: 55 anos, 70 kg, creatinina 0,9 mg/dL, ClCr %.1f mL/min\n",
              pr$crcl))
  cat(sprintf("Esquema : %d mg q%dh, infusão de %d h\n\n", dose, tau, tinf))

  cat(sprintf("Concentrações geradas: pico %.1f | vale %.1f mg/L\n\n", pico, vale))

  cat(sprintf("Prior      : CL %.3f L/h | Vc %.1f L | Vp %.1f L | Q %.1f L/h\n",
              pr$cl, pr$vc, pr$vp, pr$q))
  cat(sprintf("Verdadeiro : CL %.3f L/h | Vc %.1f L | Vss %.1f L\n",
              verdadeiro$cl, verdadeiro$vc, verdadeiro$vc + verdadeiro$vp))
  cat(sprintf("Estimado   : CL %.3f L/h | Vc %.1f L | Vss %.1f L\n\n",
              est$cl, est$vc, est$vss))

  cat(sprintf("Erro do clearance : %+.1f%%\n",
              100 * (est$cl - verdadeiro$cl) / verdadeiro$cl))
  cat(sprintf("Meia-vida alfa    : %.2f h\n", est$t12_alfa))
  cat(sprintf("Meia-vida beta    : %.2f h\n", est$t12_beta))
  cat(sprintf("AUC24 estimada    : %.0f µg×h/mL (verdadeira %.0f)\n\n",
              auc24(dose, tau, est$cl), auc24(dose, tau, verdadeiro$cl)))

  cat("Comparação com as demais implementações:\n")
  cat(sprintf("  %-28s %8s %8s\n", "", "CL", "Vss"))
  cat(sprintf("  %-28s %8.3f %8.2f\n", "verdadeiro", verdadeiro$cl,
              verdadeiro$vc + verdadeiro$vp))
  cat(sprintf("  %-28s %8.3f %8.2f\n", "esta implementação (R)", est$cl, est$vss))
  cat(sprintf("  %-28s %8.3f %8.2f\n", "TDMx", 2.841, 76.016))
  cat(sprintf("  %-28s %8.3f %8.2f\n", "VancoOptim (antes da correção)", 3.150, 76.800))

  invisible(est)
}

#' Demonstração do erro de conversão de tempo.
#'
#' Reproduz o efeito de passar o tempo do vale sem somar a duração da infusão,
#' como ocorria na implementação do VancoOptim até 08/2026.
demo_erro_tempo <- function() {
  verdadeiro <- list(cl = 3.0, vc = 40.0, q = 6.5, vp = 30.0)
  dose <- 750; tinf <- 2; tau <- 12
  t_pico <- 1.0 + tinf
  t_vale <- 9.5 + tinf
  pico <- conc_ss_2comp(t_pico, dose, tinf, tau,
                        verdadeiro$cl, verdadeiro$vc, verdadeiro$q, verdadeiro$vp)
  vale <- conc_ss_2comp(t_vale, dose, tinf, tau,
                        verdadeiro$cl, verdadeiro$vc, verdadeiro$q, verdadeiro$vp)
  pr <- prior_goti(55, 70, 0.9, "M")

  cat("\nEfeito da conversão de tempo do vale\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("%-34s %8s %8s %7s\n", "tempo do vale informado", "CL", "Vss", "erro"))
  for (rot in c("correto", "sem converter")) {
    tv <- if (rot == "correto") t_vale else 9.5
    e <- map_goti(pico, t_pico, vale, tv, dose, tinf, tau, pr)
    cat(sprintf("%-34s %8.3f %8.2f %+6.1f%%\n",
                sprintf("%s (%.1f h do início)", rot, tv),
                e$cl, e$vss, 100 * (e$cl - verdadeiro$cl) / verdadeiro$cl))
  }
  cat("\nO tempo não convertido faz o modelo ajustar a uma concentração medida\n")
  cat("mais cedo do que realmente foi, e compensa com clearance mais alto.\n")
}

# =============================================================================
# Para executar:
#   source("goti_referencia.R"); exemplo(); demo_erro_tempo()
#
# PENDENTE: conferir os parâmetros do bloco GOTI contra a publicação original.
# Enquanto isso não for feito, esta implementação verifica a aritmética da
# estimação, e não a fidelidade dos valores populacionais à fonte.
# =============================================================================
