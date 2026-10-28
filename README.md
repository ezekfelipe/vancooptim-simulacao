# VancoOptim — código de simulação

Código do estudo de simulação que integra o manuscrito *"VancoOptim: desenvolvimento,
verificação e avaliação por simulação de uma ferramenta web proprietária para dosagem
bayesiana de vancomicina guiada por AUC em ambiente hospitalar"*.

Este repositório contém apenas o código da simulação. Não inclui dados de pacientes nem o
código-fonte da ferramenta clínica, e não deve ser usado para decisão assistencial.

---

## Conteúdo do repositório

| Arquivo | O que é |
|---|---|
| `simulacao_vancooptim.R` | Estudo de simulação completo, estruturado pelo arcabouço ADEMP |
| `sensibilidade_vd_janela.R` | Sensibilidade da janela de coleta ao volume do mecanismo gerador |
| `goti_referencia.R` | Implementação de referência do modelo de Goti 2018, usada como árbitro aritmético |

### Por que existe uma implementação de referência

O `goti_referencia.R` reimplementa o modelo bicompartimental de Goti de forma explícita, sem
otimizadores de caixa-preta, para que cada passo do cálculo seja auditável a partir do código.
Serve como terceira implementação, para confronto com o motor de cálculo e com a verificação
por recálculo manual.

O alcance dessa verificação tem um limite, declarado no próprio arquivo. Os parâmetros
populacionais foram transcritos da implementação do VancoOptim e não conferidos contra a
publicação original de Goti e colaboradores. Uma implementação de referência que herda os
parâmetros da ferramenta avaliada verifica a aritmética, não a fidelidade à fonte. A
conferência dos valores contra o artigo original é passo pendente.

A justificativa é metodológica. O recálculo manual não discrimina erros de convenção, porque
quem recalcula adota as mesmas convenções de quem implementou. Foi exatamente uma divergência
de convenção — o referencial de tempo adotado na predição bicompartimental — que uma
implementação independente revelou, e que nenhuma das outras duas verificações teria
detectado.

---

## O que este código faz

Avalia se a estimativa bayesiana implementada no VancoOptim recupera a exposição verdadeira
de um paciente, sob variabilidade farmacocinética realista.

O procedimento, por réplica:

1. Sorteia um paciente virtual (peso, clearance de creatinina, clearance de vancomicina e
   volume de distribuição) a partir de um modelo populacional externo à biblioteca de priors
   implementada na ferramenta.
2. Calcula a concentração sérica que esse paciente teria no momento da coleta.
3. Adiciona erro de medição, com os componentes proporcional e aditivo do modelo gerador.
4. Submete a concentração ao estimador, que usa o prior da ferramenta, e não o modelo gerador.
5. Compara a AUC₂₄ estimada com a verdadeira, conhecida porque foi sorteada no passo 1.

Repetido 2.000 vezes por cenário, o procedimento permite medir viés, precisão e calibração
dos intervalos.

O gerador precisa ser externo. Se o modelo que gera os pacientes fosse o mesmo que o
estimador assume, a simulação avaliaria o otimizador e não o método. O gerador externo
introduz a discrepância entre modelo e realidade que existe na prática clínica.

---

## Como reproduzir

### Requisitos

R versão 4.0 ou superior, apenas com pacotes de base. Verificado em R 4.6.1 sob Windows e
reproduzido em R 4.3.3 sob Linux, com resultado idêntico em todos os campos.

### Execução

Pelo terminal, sem abrir o R:

```
Rscript -e "source('simulacao_vancooptim.R'); rodar_tudo(); calibrar_janela()"
```

A sensibilidade da janela de coleta ao volume do gerador, que produz a Tabela 12B, roda em
arquivo próprio e depende do principal estar carregado antes:

```
Rscript -e "source('simulacao_vancooptim.R'); source('sensibilidade_vd_janela.R'); sensibilidade_janela()"
```

No Windows, caso o `Rscript` não esteja no `PATH`:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "source('simulacao_vancooptim.R'); rodar_tudo(); calibrar_janela()"
```

Dentro do R ou do RStudio:

```r
source("simulacao_vancooptim.R")
res <- rodar_tudo()
cal <- calibrar_janela()
```

Para gravar a saída em arquivo:

```
Rscript -e "source('simulacao_vancooptim.R'); rodar_tudo(); calibrar_janela()" > saida.txt 2>&1
```

A execução leva alguns minutos. São 2.000 réplicas por cenário, e cada réplica roda uma
otimização numérica.

### Reprodutibilidade

A semente aleatória está fixada no topo do script:

```r
SEMENTE <- 20260812L
```

Executar o script sem alterá-la produz exatamente os mesmos números publicados no
manuscrito. Se algum valor divergir, verifique a versão de R: mudanças no gerador de números
pseudoaleatórios entre versões maiores podem alterar as sequências.

Ao reportar resultados, registre também a saída de `sessionInfo()`:

```r
res <- rodar_tudo()
res$sessao
```

---

## Correspondência com as tabelas do manuscrito

| Saída do script | Tabela no manuscrito |
|---|---|
| `rodar_tudo()`, TABELA 7 | Tabela 7, desempenho da estimativa com concentração única |
| `rodar_tudo()`, TABELA 8 | Tabela 8, volume populacional fixo versus estimado |
| `rodar_tudo()`, COMPARAÇÃO | Tabela 9, uma versus duas concentrações |
| `rodar_tudo()`, TABELA 11 | Tabela 11, desempenho sob o gerador de Llopis-Salvia |
| `calibrar_janela()` | Tabela 12, erro por momento da coleta |
| `sensibilidade_janela()` | Tabela 12B, janela de coleta por volume verdadeiro do gerador |

A Tabela 10, que trata da calibração do σ_obs, foi produzida por rotina separada e não está
incluída nesta versão do repositório.

A análise de comparação entre modelos populacionais, reportada na seção 3.8 do manuscrito,
também não integra esta versão. Ela invoca diretamente o motor de cálculo da ferramenta, e
não uma reimplementação independente, porque exige as definições dos seis priors tal como
implementadas — transcrevê-las criaria uma segunda fonte de verdade para os mesmos
coeficientes. Depositá-la aqui implicaria publicar módulos do código clínico, o que está fora
do escopo declarado deste repositório.

---

## Estrutura do script

O arquivo segue a ordem do arcabouço ADEMP (Morris, White & Crowther, *Stat Med* 2019), que
estrutura o reporte de estudos de simulação em objetivos, mecanismo gerador, estimandos,
métodos e medidas de desempenho.

### Mecanismo gerador: `ROBERTS`, `COORTE`

Parâmetros do modelo de Roberts e colaboradores (2011), reproduzidos da Tabela 2 da
publicação original: clearance típico, volume por quilo, variabilidade entre sujeitos para
ambos os parâmetros, e erro residual proporcional somado a aditivo. As covariáveis são
sorteadas conforme a coorte descrita naquele artigo.

Há uma ressalva importante, documentada também no código. O modelo de Roberts deriva de
infusão contínua, condição em que o volume de distribuição é fracamente identificável: no
equilíbrio, a concentração corresponde à razão entre taxa e clearance, e o volume só é
informado pelas fases de transição. Por isso o volume verdadeiro entra como fator do
desenho, e não como constante. Os cenários varrem de 1,53 a 0,55 L/kg, o que permite separar
o efeito da má especificação do efeito do método.

### Modelo estrutural: `conc_ss()`

Concentração no estado estacionário, um compartimento, infusão intermitente:

```
C(t) = (Dose/(T·CL)) · (1 − e^(−k·T)) · e^(−k·t) / (1 − e^(−k·τ))
```

A divisão por `T` merece atenção. A fórmula parte da taxa de infusão, e não da dose. As duas
formas coincidem quando `T` vale 1 hora, e essa coincidência já mascarou um defeito de
implementação neste projeto. O comentário no código registra o episódio.

### Métodos avaliados

`map_1ponto()` estima o clearance a partir de uma concentração, por maximização da
probabilidade a posteriori, com o volume fixo no valor populacional. A incerteza vem da
curvatura da função objetivo: como `obj = −2·log(posterior)`, a covariância é `2·H⁻¹`. O
fator 2 é a armadilha clássica dessa derivação, e omiti-lo estreita o intervalo por 1,41.

`map_1ponto_vd_livre()` faz a mesma estimativa permitindo que o volume varie. Existe apenas
para testar a decisão de projeto, e alimenta a Tabela 8.

`propagar_vd()` propaga a incerteza do volume populacional para o clearance, pelo método
delta. Sem essa correção a cobertura dos intervalos fica muito abaixo do nominal, porque o
volume entra na estimação como se fosse conhecido.

`map_2pontos()` implementa Sawchuk-Zaske seguido de atualização Normal-Normal conjugada,
como na rota principal da ferramenta.

### Prior de estimação: `prior_estratificado()`

Modelo estratificado por clearance de creatinina, que é o prior padrão da ferramenta. Não
confundir com o gerador: este é o que o estimador assume, e a diferença entre os dois é o
que a simulação avalia.

### Medidas de desempenho: `desempenho()`

Além de viés e precisão, a função devolve o erro de Monte Carlo, que mede a imprecisão da
própria simulação. É calculado como o desvio-padrão dos erros dividido pela raiz do número
de réplicas.

Reportá-lo é o que permite afirmar que uma diferença entre cenários é real. Com 2.000
réplicas, o erro de Monte Carlo fica abaixo de 0,7 ponto percentual em todos os cenários,
o que sustenta comparações da ordem de 1 a 2 pontos.

### Calibração local: `calibrar_janela()`, `K_OBSERVADOS_MORIAH`

Esta rotina roda com a distribuição de constantes de eliminação estimadas na população
assistida, e não com a do gerador. É o que sustenta o achado de que a recomendação de janela
de coleta depende do perfil populacional.

O vetor contém 11 valores de constante de eliminação, derivados de registros de casos reais e
integralmente anonimizados. São constantes farmacocinéticas, sem qualquer identificador, data
ou dado demográfico associado. Por derivarem de registros com concentração única, e não de
pares de concentrações, são estimativas condicionadas ao modelo populacional em vigor no
momento do cálculo, e não medidas diretas.

Além da distribuição de constantes de eliminação, a rotina sorteia o volume de distribuição
verdadeiro de uma lognormal centrada em 0,60 L/kg com coeficiente de variação de 30%. O
resultado depende desse valor. O `sensibilidade_vd_janela.R` repete a análise sobre quatro
volumes verdadeiros, mantido o prior do estimador em 0,70 L/kg: o instante ótimo permanece na
janela de 4 a 5 horas de 0,60 a 1,10 L/kg, mas as magnitudes do erro crescem com ele, e no
cenário de 1,10 L/kg o erro do protocolo institucional ultrapassa o critério de
aceitabilidade adotado.

---

## Limitações do desenho

Estão documentadas no manuscrito e resumidas aqui.

Os dois modelos geradores derivam de populações críticas europeias, e nenhum deles foi
derivado de coorte brasileira. O de Roberts deriva de infusão contínua, condição em que o
volume de distribuição é fracamente identificável. O de Llopis-Salvia deriva de administração
intermitente, que é a modalidade contemplada pela ferramenta, e tem estrutura de dois
compartimentos, o que introduz má especificação estrutural além da de parâmetro.

Os intervalos da estimativa com concentração única cobriram de 69,5% a 82,0%, contra 95%
nominal, em todos os cenários avaliados. A causa predominante é a discrepância sistemática
entre o volume populacional e o individual, que não pode ser capturada a partir de uma única
concentração.

A simulação avalia a recuperação de parâmetros conhecidos. Não substitui a comparação entre
concentrações previstas e observadas em pacientes reais.

---

## Licença

O código é distribuído sob a licença **MIT**, cujo texto integral está em `LICENSE`.

A documentação — este README e demais textos explicativos — é distribuída sob **Creative
Commons Atribuição 4.0 Internacional (CC BY 4.0)**. Ver `LICENSE-docs.md`.

## Citação

Ao usar este código, cite o registro do repositório e o manuscrito associado.

A citação com DOI será inserida aqui após o depósito no Zenodo. O identificador a citar é o
**DOI de versão**, e não o de conceito: o de conceito aponta sempre para a versão mais
recente, o que quebraria a correspondência com os números publicados no manuscrito.

---

## Aviso

Código de pesquisa. Não constitui dispositivo médico, não foi submetido a avaliação
regulatória, e não deve ser empregado em decisão sobre dose de paciente. A ferramenta
clínica a que este estudo se refere está sujeita a validação farmacêutica local em todo uso
assistencial.
