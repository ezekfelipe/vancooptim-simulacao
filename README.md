# VancoOptim — código de simulação

Código do estudo de simulação que integra o manuscrito *"VancoOptim: desenvolvimento,
verificação e avaliação por simulação de uma ferramenta web proprietária para dosagem
bayesiana de vancomicina guiada por AUC em ambiente hospitalar"*.

Este repositório contém apenas o código da simulação. Não inclui dados de pacientes nem o
código-fonte da ferramenta clínica, e não deve ser usado para decisão assistencial.

---

## O que este código faz

Avalia se a estimativa bayesiana implementada no VancoOptim recupera a exposição
verdadeira de um paciente, sob variabilidade farmacocinética realista.

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

R versão 4.0 ou superior, apenas com pacotes de base. Verificado em R 4.6.1 sob Windows.

### Execução

Pelo terminal, sem abrir o R:

```
Rscript -e "source('simulacao_vancooptim.R'); rodar_tudo(); calibrar_janela()"
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
| `calibrar_janela()` | Tabela 11, erro por momento da coleta |

A Tabela 10, que trata da calibração do σ_obs, foi produzida por rotina separada e não está
incluída nesta versão do repositório.

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

Esta rotina roda com a distribuição de constantes de eliminação observada na população
assistida, e não com a do gerador. É o que sustenta o achado de que a recomendação de janela
de coleta depende do perfil populacional: simulações preliminares com parâmetros genéricos
produziram estimativas substancialmente diferentes para o mesmo efeito.

O vetor contém 11 valores de constante de eliminação, derivados de casos reais e
integralmente anonimizados. São constantes farmacocinéticas, sem qualquer identificador,
data ou dado demográfico associado.

---

## Limitações do desenho

Estão documentadas no manuscrito e resumidas aqui.

O estudo empregou um único modelo gerador. O modelo de Roberts deriva de infusão contínua,
e um segundo gerador derivado de administração intermitente fortaleceria o desenho. O
critério para essa escolha é que a publicação de origem reporte efeitos fixos, variabilidade
entre sujeitos e erro residual. Sem os três componentes o gerador fica incompleto.

Os intervalos da estimativa com concentração única cobriram de 59,6% a 79,5%, contra 95%
nominal, em todos os cenários avaliados. A causa predominante é a discrepância sistemática
entre o volume populacional e o individual, que não pode ser capturada a partir de uma única
concentração.

A simulação avalia a recuperação de parâmetros conhecidos. Não substitui a comparação entre
concentrações previstas e observadas em pacientes reais.

---

## Licença e citação

[A DEFINIR. Sugestão: MIT para o código, CC-BY 4.0 para a documentação.]

Ao usar este código, cite o registro do repositório e o manuscrito associado.

```
[A DEFINIR. Inserir a citação com DOI após o depósito no Zenodo.]
```

---

## Aviso

Código de pesquisa. Não constitui dispositivo médico, não foi submetido a avaliação
regulatória, e não deve ser empregado em decisão sobre dose de paciente. A ferramenta
clínica a que este estudo se refere está sujeita a validação farmacêutica local em todo uso
assistencial.
