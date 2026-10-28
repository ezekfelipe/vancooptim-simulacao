# =============================================================================
# VancoOptim — sensibilidade da janela de coleta ao volume de distribuição
#              verdadeiro do mecanismo gerador (Tabela 12B do manuscrito)
#
# PERGUNTA: a recomendação de antecipar a coleta depende do volume de
# distribuição assumido na geração dos pacientes?
#
# MOTIVO. A `calibrar_janela()` do `simulacao_vancooptim.R` sorteia o volume
# verdadeiro de uma lognormal centrada em 0,60 L/kg. Esse valor situa-se abaixo
# do prior adotado pela ferramenta, de 0,70 L/kg, e abaixo da faixa relatada
# para adultos de peso normal. Como o achado da janela de coleta é apresentado
# como específico da população assistida, convém separar o que nele é robusto
# do que é condicional a essa escolha.
#
# O que varia aqui é apenas o volume VERDADEIRO do gerador. O prior que o
# estimador assume permanece em 0,70 L/kg nos quatro cenários, como no motor.
#
# REPRODUTIBILIDADE: mesma semente do estudo principal, reposicionada no início
# de cada cenário, de modo que os quatro cenários veem a mesma sequência de
# pacientes e diferem somente pelo volume. Rodar este arquivo produz sempre os
# mesmos números.
#
# Uso:   source("simulacao_vancooptim.R")
#        source("sensibilidade_vd_janela.R")
#        sensibilidade_janela()
#
#   ou, pelo terminal:
#        Rscript -e "source('simulacao_vancooptim.R'); source('sensibilidade_vd_janela.R'); sensibilidade_janela()"
#
# Requisitos: R base, e o `simulacao_vancooptim.R` carregado antes.
# =============================================================================

#' Repete a calibração da janela de coleta para um volume verdadeiro dado.
#'
#' Reproduz `calibrar_janela()` em tudo, exceto por `vd_ger`, que ali está fixo
#' em 0,60 L/kg. Devolve o erro percentual médio absoluto da AUC24 para cada
#' hora da janela útil.
janela_por_vd <- function(vd_ger, tau = 12, tinf = 2,
                          n = N_REPLICAS, dose = 1000) {
  set.seed(SEMENTE + 3000L)
  saida <- numeric(0)

  for (tc in seq_len(tau - tinf)) {
    erros <- numeric(0)
    for (i in seq_len(n)) {
      peso   <- max(35, rnorm(1, COORTE$peso_media, COORTE$peso_dp))
      k_base <- sample(K_OBSERVADOS_MORIAH, 1)
      k_v    <- rlnorm_cv(1, k_base, 0.20)
      vd_v   <- rlnorm_cv(1, vd_ger * peso, 0.30)
      cl_v   <- k_v * vd_v
      auc_v  <- dose * (24 / tau) / cl_v

      c_real <- conc_ss(cl_v, vd_v, dose, tinf, tau, tc)
      if (!is.finite(c_real)) next
      c_obs <- c_real * (1 + 0.08 * rnorm(1))
      if (c_obs < ROBERTS$loq) next

      # Prior do estimador: adulto típico, invariante entre os cenários.
      pr  <- prior_estratificado(90, peso)
      est <- propagar_vd(c_obs, tc, dose, tinf, tau,
                         pr$cl, pr$cl * 0.30,
                         pr$vd_kg * peso, pr$vd_kg * peso * 0.30)
      auc_e <- dose * (24 / tau) / est$cl
      if (is.finite(auc_e)) erros <- c(erros, 100 * (auc_e - auc_v) / auc_v)
    }
    saida <- c(saida, mean(abs(erros)))
  }
  saida
}

#' Imprime a Tabela 12B do manuscrito.
sensibilidade_janela <- function(cenarios = c(0.60, 0.70, 0.90, 1.10),
                                 tau = 12, tinf = 2) {

  res <- sapply(cenarios, janela_por_vd, tau = tau, tinf = tinf)

  cat("\n\nTABELA 12B — janela de coleta por volume verdadeiro do gerador\n")
  cat("Prior do estimador fixo em 0,70 L/kg | semente", SEMENTE,
      "| réplicas por cenário:", N_REPLICAS, "\n")
  cat(sprintf("%-22s", "h após fim da infusão"))
  for (v in cenarios) cat(sprintf("%12s", paste0(formatC(v, format = "f", digits = 2), " L/kg")))
  cat("\n"); cat(strrep("-", 22 + 12 * length(cenarios)), "\n")

  for (j in seq_len(nrow(res))) {
    cat(sprintf("%-22s", paste0(j, " h")))
    for (i in seq_len(ncol(res))) {
      cat(sprintf("%11s%%", formatC(res[j, i], format = "f", digits = 2)))
    }
    cat("\n")
  }

  cat("\nInstante ótimo por cenário:\n")
  for (i in seq_along(cenarios)) {
    cat(sprintf("  Vd verdadeiro %.2f L/kg -> %d h (MAPE %.2f%%)\n",
                cenarios[i], which.min(res[, i]), min(res[, i])))
  }
  cat("\nA coleta às 9 h corresponde ao protocolo institucional vigente.\n")

  invisible(list(cenarios = cenarios, mape = res,
                 semente = SEMENTE, replicas = N_REPLICAS))
}
