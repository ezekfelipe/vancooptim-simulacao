# =============================================================================
# REEXTRAÇÃO DAS CONSTANTES DE ELIMINAÇÃO SOB O PRIOR VIGENTE
#
# RODE ESTE ARQUIVO NA INSTITUIÇÃO, E NÃO AQUI.
#
# As planilhas da população assistida contêm nome completo, prontuário e
# creatininas seriadas. Nada disso entra neste repositório, em artefato, em
# documento de projeto ou em arquivo entregue. Este script foi escrito para que
# a reextração ocorra na máquina de quem tem acesso legítimo aos registros, e
# para que apenas o vetor de onze números saia dali.
#
# ENTRADA. Um CSV desidentificado, preparado a partir da planilha, com uma
# linha por paciente e estas colunas, sem nome, sem prontuário e sem data:
#
#   peso          kg
#   clcr          mL/min, pela mesma fórmula que a ferramenta emprega
#   dose          mg por administração
#   tau           h, intervalo entre doses
#   tinf          h, duração da infusão, derivada de início e término
#   t_coleta      h após o TÉRMINO da infusão
#   concentracao  mg/L
#
# O arquivo de entrada não é versionado. Mantenha-o fora do repositório ou
# acrescente o seu nome ao .gitignore antes de criá-lo.
#
# SAÍDA. O vetor de constantes pronto para substituir `K_OBSERVADOS_MORIAH` em
# `simulacao_vancooptim.R`, mais a mediana e os quartis da meia-vida, que são
# os valores citados na seção 3.6 do manuscrito. Nenhuma linha individual é
# impressa com os dados de entrada.
#
# O QUE O SCRIPT FAZ. Para cada paciente, aplica a MESMA rota de concentração
# única que o motor aplica: fixa o volume de distribuição no valor populacional
# do prior estratificado vigente e estima o clearance por máximo a posteriori.
# A constante é o quociente entre o clearance estimado e o volume assumido,
# k = CL_est / Vd_pop, que é a definição sob a qual as constantes publicadas
# foram obtidas, agora com o volume correto.
#
# USO
#   Rscript reextrair_k.R caminho/para/coorte_desidentificada.csv
# =============================================================================

suppressWarnings(source("simulacao_vancooptim.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("informe o CSV desidentificado: Rscript reextrair_k.R <arquivo.csv>",
       call. = FALSE)
}
caminho <- args[[1]]
if (!file.exists(caminho)) stop("arquivo não encontrado: ", caminho, call. = FALSE)

d <- read.csv(caminho, stringsAsFactors = FALSE)

obrigatorias <- c("peso", "clcr", "dose", "tau", "tinf", "t_coleta", "concentracao")
faltando <- setdiff(obrigatorias, names(d))
if (length(faltando)) {
  stop("colunas ausentes no CSV: ", paste(faltando, collapse = ", "), call. = FALSE)
}

# Recusa colunas que não deveriam existir, em vez de ignorá-las em silêncio.
proibidas <- grep("nome|prontuario|prontuário|registro|paciente|data|nasc|cpf",
                  names(d), ignore.case = TRUE, value = TRUE)
if (length(proibidas)) {
  stop("o CSV traz coluna(s) que não devem sair da instituição: ",
       paste(proibidas, collapse = ", "),
       "\nremova-as do arquivo e repita.", call. = FALSE)
}

k <- numeric(0)
descartados <- 0
for (i in seq_len(nrow(d))) {
  r <- d[i, ]
  pr <- prior_estratificado(r$clcr, r$peso)
  vd_pop <- pr$vd_kg * r$peso
  est <- map_1ponto(r$concentracao, r$t_coleta, r$dose, r$tinf, r$tau,
                    pr$cl, pr$cl * 0.30, vd_pop)
  ki <- est$cl / vd_pop
  if (is.finite(ki) && ki > 0) k <- c(k, ki) else descartados <- descartados + 1
}

if (!length(k)) stop("nenhuma linha produziu constante válida.", call. = FALSE)

t12 <- sort(log(2) / k)
cat("\nREEXTRAÇÃO SOB O PRIOR VIGENTE (Vd 0,70 L/kg)\n")
cat(sprintf("Linhas lidas: %d | constantes válidas: %d | descartadas: %d\n\n",
            nrow(d), length(k), descartados))

cat("Vetor para substituir K_OBSERVADOS_MORIAH em simulacao_vancooptim.R:\n\n")
val <- sprintf("%.4f", k)
linhas <- split(val, ceiling(seq_along(val) / 5))
cat("K_OBSERVADOS_MORIAH <- c(", paste(linhas[[1]], collapse = ", "),
    if (length(linhas) > 1) "," else ")", "\n", sep = "")
if (length(linhas) > 1) {
  for (j in seq(2, length(linhas))) {
    cat("                         ", paste(linhas[[j]], collapse = ", "),
        if (j < length(linhas)) "," else ")", "\n", sep = "")
  }
}

cat(sprintf("\nMeia-vida (n = %d): mediana %.1f h | quartis %.1f a %.1f h\n",
            length(t12), median(t12), quantile(t12, 0.25), quantile(t12, 0.75)))
cat("São os valores citados na seção 3.6 do manuscrito.\n")
cat("\nDepois de substituir o vetor, reexecute calibrar_janela() e confira a\n")
cat("Tabela 12. O ótimo deve permanecer em 5 h; o nível do erro deve cair,\n")
cat("conforme sensibilidade_k_janela.R previu.\n\n")
