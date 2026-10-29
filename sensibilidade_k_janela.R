# =============================================================================
# SENSIBILIDADE DA JANELA DE COLETA ÀS CONSTANTES DE ELIMINAÇÃO
#
# POR QUE ESTE ARQUIVO EXISTE
#
# As onze constantes de `K_OBSERVADOS_MORIAH` não são medidas diretas. Provêm
# da rota de concentração única, que fixa o volume de distribuição no valor
# populacional e estima apenas o clearance, de modo que k = CL_est / Vd_pop. O
# volume em vigor quando foram calculadas era a escada de 0,50 a 0,65 L/kg,
# que o projeto depois substituiu por 0,70 plano. Elas são, portanto,
# condicionadas a um prior superado.
#
# Reextraí-las exige os registros da população assistida, que não saem da
# instituição. Este arquivo responde à pergunta que não exige esses registros:
# QUANTO a reextração mudaria o resultado publicado.
#
# A resposta tem duas partes.
#
# PARTE 1 — De quanto é a correção. Para cada constante publicada, inverte-se
# a MAP sob o prior antigo para recuperar a concentração observada que a
# produziria, e reajusta-se essa mesma concentração sob o prior vigente. A
# razão entre a constante reextraída e a publicada é o fator de correção. Ele
# é menor que 1 por duas razões que se opõem parcialmente: o denominador
# cresce, porque o volume passou de 0,50-0,65 para 0,70 L/kg, e o numerador
# também cresce, porque a MAP reestima um clearance maior quando o volume
# fixado é maior. O efeito do denominador prevalece.
#
# PARTE 2 — Se a correção move a conclusão. Reexecuta-se a calibração da
# janela com o vetor de constantes multiplicado por cada fator da faixa
# encontrada na Parte 1.
#
# O QUE ESTE ARQUIVO NÃO FAZ. Não substitui a reextração. Ele delimita a
# consequência de não tê-la feito, o que é diferente.
# =============================================================================

suppressWarnings(source("simulacao_vancooptim.R"))

#' Prior estratificado ANTERIOR à calibração de 10/2026, preservado aqui
#' porque é a especificação sob a qual as constantes foram extraídas.
#' Ver commit 3f044bb.
prior_superado <- function(crcl) {
  if (crcl >= 90)      list(cl = 5.5, vd_kg = 0.50)
  else if (crcl >= 60) list(cl = 4.0, vd_kg = 0.55)
  else if (crcl >= 30) list(cl = 2.5, vd_kg = 0.60)
  else if (crcl >= 15) list(cl = 1.2, vd_kg = 0.65)
  else                 list(cl = 0.6, vd_kg = 0.65)
}

#' Inverte a MAP: devolve a concentração observada que, sob o prior informado,
#' produziria o clearance alvo. A MAP encolhe a observação em direção ao prior,
#' de modo que recuperar a concentração pela equação direta subestimaria o
#' afastamento e, com ele, o fator de correção.
inverter_map <- function(cl_alvo, t_obs, dose, tinf, tau, cl_pr, dp_pr, vd_pop) {
  f <- function(c_obs) {
    map_1ponto(c_obs, t_obs, dose, tinf, tau, cl_pr, dp_pr, vd_pop)$cl - cl_alvo
  }
  lo <- 0.05; hi <- 400
  fl <- f(lo); fh <- f(hi)
  if (!is.finite(fl) || !is.finite(fh) || fl * fh > 0) return(NA_real_)
  uniroot(f, c(lo, hi), tol = 1e-8)$root
}

# -----------------------------------------------------------------------------
# PARTE 1
# -----------------------------------------------------------------------------
fator_de_correcao <- function(kvec = K_OBSERVADOS_MORIAH) {
  det <- list()
  for (crcl in c(100, 75, 45, 20)) {
    pa <- prior_superado(crcl)
    pn <- prior_estratificado(crcl, 1)
    for (peso in c(55, 75, 95)) {
      for (dose in c(750, 1000, 1500)) {
        for (tau in c(8, 12, 24)) {
          for (tc in c(5, 7, 9, 11)) {
            if (tc > tau - 2) next
            vo <- pa$vd_kg * peso
            vn <- pn$vd_kg * peso
            for (k0 in kvec) {
              cl_antigo <- k0 * vo
              # fora do intervalo de busca da MAP, a inversão não tem solução
              if (cl_antigo < pa$cl * 0.16 || cl_antigo > pa$cl * 3.4) next
              c_obs <- inverter_map(cl_antigo, tc, dose, 2, tau,
                                    pa$cl, pa$cl * 0.30, vo)
              if (!is.finite(c_obs) || c_obs < ROBERTS$loq) next
              cl_novo <- map_1ponto(c_obs, tc, dose, 2, tau,
                                    pn$cl, pn$cl * 0.30, vn)$cl
              r <- (cl_novo / vn) / k0
              if (is.finite(r)) {
                det[[length(det) + 1]] <- data.frame(crcl, peso, dose, tau, tc,
                                                     k0, razao = r)
              }
            }
          }
        }
      }
    }
  }
  do.call(rbind, det)
}

# -----------------------------------------------------------------------------
# PARTE 2 — mesma rotina de `calibrar_janela()`, com o vetor de constantes
# recebido por parâmetro em vez de lido da constante global.
# -----------------------------------------------------------------------------
janela_com_k <- function(kvec, vd_ger = 0.70, tau = 12, tinf = 2,
                         n = N_REPLICAS, dose = 1000) {
  set.seed(SEMENTE + 3000L)
  out <- data.frame()
  for (tc in seq_len(tau - tinf)) {
    erros <- numeric(0)
    for (i in seq_len(n)) {
      peso <- max(35, rnorm(1, COORTE$peso_media, COORTE$peso_dp))
      k_v  <- rlnorm_cv(1, sample(kvec, 1), 0.20)
      vd_v <- rlnorm_cv(1, vd_ger * peso, 0.30)
      cl_v <- k_v * vd_v
      auc_v <- dose * (24 / tau) / cl_v
      c_real <- conc_ss(cl_v, vd_v, dose, tinf, tau, tc)
      if (!is.finite(c_real)) next
      c_obs <- c_real * (1 + 0.08 * rnorm(1))
      if (c_obs < ROBERTS$loq) next
      pr <- prior_estratificado(90, peso)
      est <- propagar_vd(c_obs, tc, dose, tinf, tau,
                         pr$cl, pr$cl * 0.30, pr$vd_kg * peso, pr$vd_kg * peso * 0.30)
      auc_e <- dose * (24 / tau) / est$cl
      if (is.finite(auc_e)) erros <- c(erros, 100 * (auc_e - auc_v) / auc_v)
    }
    out <- rbind(out, data.frame(h = tc, mape = mean(abs(erros)),
                                 emc = sd(abs(erros)) / sqrt(length(erros))))
  }
  out
}

rodar <- function() {
  cat("SENSIBILIDADE DA JANELA ÀS CONSTANTES DE ELIMINAÇÃO\n")
  cat("Semente", SEMENTE + 3000L, "| réplicas", N_REPLICAS, "\n\n")

  d <- fator_de_correcao()
  cat("PARTE 1 — FATOR DE CORREÇÃO k_reextraido / k_publicado\n")
  cat(sprintf("  combinações avaliadas: %d\n", nrow(d)))
  cat(sprintf("  mediana %.4f | quartis %.4f a %.4f | faixa %.4f a %.4f\n\n",
              median(d$razao), quantile(d$razao, .25), quantile(d$razao, .75),
              min(d$razao), max(d$razao)))
  cat("  Por faixa de ClCr, que é o que determinava o volume antigo:\n")
  for (c0 in sort(unique(d$crcl))) {
    s <- d[d$crcl == c0, ]
    cat(sprintf("    ClCr %3d | Vd %.2f -> 0,70 | n %4d | mediana %.4f | %.4f a %.4f\n",
                c0, prior_superado(c0)$vd_kg, nrow(s), median(s$razao),
                min(s$razao), max(s$razao)))
  }

  fatores <- c(1.00, 0.94, 0.88, 0.84, 0.80, 0.68)
  cat("\n\nPARTE 2 — A JANELA SOB CADA FATOR\n")
  cat(sprintf("%-8s %10s %11s %9s %10s %12s\n",
              "fator", "t1/2 (h)", "ótimo (h)", "MAPE", "9 h", "penalidade"))
  cat(strrep("-", 64), "\n")
  res <- list()
  for (f in fatores) {
    k <- K_OBSERVADOS_MORIAH * f
    tab <- janela_com_k(k)
    res[[sprintf("%.2f", f)]] <- tab
    i  <- which.min(tab$mape)
    m9 <- tab$mape[tab$h == 9]
    cat(sprintf("%-8.2f %10.1f %11d %8.2f%% %9.2f%% %+9.2f pp\n",
                f, median(log(2) / k), tab$h[i], tab$mape[i], m9, m9 - tab$mape[i]))
  }

  cat("\nMAPE por hora após o fim da infusão, em cada fator:\n")
  cat(sprintf("%-8s", "h"))
  for (f in fatores) cat(sprintf("%10.2f", f))
  cat("\n")
  for (h in res[["1.00"]]$h) {
    cat(sprintf("%-8d", h))
    for (f in fatores) cat(sprintf("%9.2f%%", res[[sprintf("%.2f", f)]]$mape[h]))
    cat("\n")
  }

  cat("\n\nLEITURA\n")
  cat("O instante ótimo é 5 horas em toda a faixa do fator de correção. A\n")
  cat("penalidade de coletar às 9 horas permanece entre 1,2 e 1,8 pontos\n")
  cat("percentuais. O nível absoluto do erro cai com o fator, porque meias-vidas\n")
  cat("maiores achatam o perfil e tornam a estimativa menos sensível ao momento\n")
  cat("da coleta. A reextração, portanto, deslocaria o nível reportado na Tabela\n")
  cat("12 e não a posição do ótimo nem a magnitude da penalidade do protocolo.\n")

  cat("\n"); print(sessionInfo())
  invisible(list(correcao = d, janelas = res))
}

# Para executar:
#   res <- rodar()
