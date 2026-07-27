# Análise de Séries Temporais

## 1. Componentes de uma Série Temporal

Uma série temporal pode ser decomposta em quatro componentes principais:

### Tendência (Trend)
Movimento de longo prazo da série — direção geral ao longo do tempo (crescimento, queda ou estabilidade). Pode ser linear ou não-linear.

### Sazonalidade (Seasonality)
Padrões que se repetem em intervalos fixos e conhecidos (diário, semanal, mensal, anual). Exemplo: vendas maiores no Natal todo ano.

### Ciclicidade (Cyclicity)
Flutuações de longo prazo sem período fixo, geralmente associadas a ciclos econômicos. Diferente da sazonalidade porque o intervalo não é regular.

### Resíduo / Irregularidade (Residual)
Variação aleatória que não é explicada pelos outros componentes — ruído.

### Modelos de Decomposição

| Modelo | Fórmula | Quando usar |
|--------|---------|-------------|
| Aditivo | `Y = T + S + C + R` | Amplitude da sazonalidade constante |
| Multiplicativo | `Y = T × S × C × R` | Amplitude da sazonalidade cresce com a tendência |

---

## 2. Estacionaridade

### O que é?
Uma série é **estacionária** quando suas propriedades estatísticas não mudam ao longo do tempo:
- Média constante
- Variância constante
- Autocovariância depende apenas do lag, não do tempo

### Por que importa?
A maioria dos modelos clássicos (ARMA, por exemplo) assume estacionaridade. Aplicar esses modelos em séries não-estacionárias gera resultados espúrios.

### Tipos de Não-Estacionaridade

**Não-estacionária com tendência (Trend-Stationary):**
A série tem tendência determinística. Basta remover a tendência (regressão) para torná-la estacionária.

**Não-estacionária com raiz unitária (Difference-Stationary):**
A série precisa ser diferenciada para se tornar estacionária. É o caso mais comum na prática (passeio aleatório).

### Testes de Estacionaridade

| Teste | Hipótese Nula (H₀) | Rejeitar H₀ se |
|-------|-------------------|----------------|
| **ADF** (Augmented Dickey-Fuller) | Série tem raiz unitária (não-estacionária) | p-value < 0.05 → estacionária |
| **KPSS** | Série é estacionária | p-value < 0.05 → não-estacionária |
| **PP** (Phillips-Perron) | Série tem raiz unitária | p-value < 0.05 → estacionária |

> Boa prática: usar ADF e KPSS juntos — conclusões opostas indicam estacionaridade de tendência.

### Como tornar uma série estacionária?

O processo segue uma ordem lógica: primeiro estabilizar a **variância**, depois remover a **tendência**, depois tratar a **sazonalidade**.

#### Passo 1 — Estabilizar a variância (se necessário)

Se a amplitude das oscilações cresce com o tempo (variância não-constante), aplica-se uma transformação antes de qualquer diferenciação:

| Transformação | Quando usar | Efeito |
|---------------|-------------|--------|
| `log(yₜ)` | Crescimento exponencial, dados sempre positivos | Comprime valores grandes, lineariza crescimento |
| `sqrt(yₜ)` | Variância proporcional à média (dados de contagem) | Efeito mais suave que o log |
| `Box-Cox(yₜ, λ)` | Caso geral — encontra o λ ótimo automaticamente | Generaliza log (λ=0) e raiz (λ=0.5) |

> Se o gráfico da série mostrar "funil" (dispersão crescendo ao longo do tempo), aplicar log ou Box-Cox antes de diferenciar.

#### Passo 2 — Remover a tendência

**Opção A — Diferenciação (série com raiz unitária):**

É a abordagem mais comum. Calcula a variação entre observações consecutivas em vez dos valores absolutos:

```
Δyₜ = yₜ - yₜ₋₁          ← 1ª diferença (d=1)
Δ²yₜ = Δyₜ - Δyₜ₋₁       ← 2ª diferença (d=2), raramente necessária
```

- **d=1** resolve a maioria dos casos (tendência linear).
- **d=2** é necessário quando a taxa de crescimento em si ainda tem tendência.
- Diferenciar demais (over-differencing) introduz autocorrelação artificial — pare quando o ADF indicar estacionaridade.

**Opção B — Remoção de tendência determinística:**

Se a tendência é determinística (previsível, não aleatória), basta ajustar uma regressão e trabalhar com os resíduos:

```
yₜ = β₀ + β₁t + εₜ   →   série estacionária = εₜ = yₜ - ŷₜ
```

Quando a tendência é não-linear, usa-se regressão polinomial ou a componente de tendência extraída pela decomposição STL.

#### Passo 3 — Remover a sazonalidade (se houver)

**Diferenciação sazonal:**

Remove padrões que se repetem a cada `s` períodos (ex: s=12 para dados mensais, s=7 para diários):

```
Δₛyₜ = yₜ - yₜ₋ₛ
```

Muitas vezes combina-se diferenciação regular e sazonal:

```
Δ Δ₁₂ yₜ = (yₜ - yₜ₋₁₂) - (yₜ₋₁ - yₜ₋₁₃)   ← usada no SARIMA
```

**Alternativa:** subtrair o índice sazonal calculado pela decomposição (seasonally adjusted series).

#### Resumo do fluxo de transformação

```
Série original
    │
    ├─ Variância instável? ──► log(yₜ) ou Box-Cox
    │
    ├─ Tem tendência? ─────────────────────────────┐
    │       │                                       │
    │  Raiz unitária           Tendência determinística
    │  (ADF não rejeita H₀)   (regressão visível)
    │       │                       │
    │   diferenciar (d=1)    subtrair ŷₜ
    │
    └─ Tem sazonalidade? ──► diferenciar sazonalmente (Δₛ)
    │
    ▼
 Série estacionária → confirmar com ADF/KPSS → modelar
```

> Regra prática: aplicar apenas o mínimo de transformações necessárias para passar nos testes — cada operação descarta informação.

---

## 3. Modelos de Séries Temporais

### Modelos Clássicos (Lineares)

#### AR — Autorregressivo
A variável depende de seus próprios valores passados.
```
yₜ = c + φ₁yₜ₋₁ + φ₂yₜ₋₂ + ... + φₚyₜ₋ₚ + εₜ
```
- Parâmetro: **p** (ordem — quantos lags usar)
- Identifica pelo: PACF (Partial Autocorrelation Function)

#### MA — Médias Móveis
A variável depende dos erros passados.
```
yₜ = μ + εₜ + θ₁εₜ₋₁ + θ₂εₜ₋₂ + ... + θqεₜ₋q
```
- Parâmetro: **q** (ordem)
- Identifica pelo: ACF (Autocorrelation Function)

#### ARMA — AR + MA
Combina os dois. Exige série **estacionária**.
```
yₜ = c + φ₁yₜ₋₁ + ... + θ₁εₜ₋₁ + ... + εₜ
```

#### ARIMA — ARMA + Integração
Adiciona diferenciação para lidar com não-estacionaridade.
- **d**: número de diferenciações necessárias
- Notação: `ARIMA(p, d, q)`

#### SARIMA — ARIMA com Sazonalidade
Estende o ARIMA com termos sazonais.
- Notação: `SARIMA(p, d, q)(P, D, Q)ₛ`
- Parâmetros maiúsculos = componentes sazonais

### Modelos de Suavização Exponencial

| Modelo | Captura |
|--------|---------|
| SES (Simple Exponential Smoothing) | Nível |
| Holt (Double) | Nível + Tendência |
| Holt-Winters (Triple) | Nível + Tendência + Sazonalidade |

### Modelos Modernos / Machine Learning

| Modelo | Características |
|--------|----------------|
| **Prophet** (Meta) | Robusto a outliers e feriados, fácil de usar, decomposição aditiva/multiplicativa |
| **LSTM / GRU** | Redes neurais recorrentes, captura padrões complexos e não-lineares |
| **Transformer (PatchTST, etc.)** | Atenção temporal, estado da arte para séries longas |
| **XGBoost / LightGBM** | Requer feature engineering manual (lags, rolling stats) |
| **N-BEATS / N-HiTS** | Redes neurais especializadas em séries temporais, interpretáveis |

---

## 4. Comparativo Rápido dos Modelos

| Modelo | Estacionaridade necessária | Sazonalidade | Não-linearidade | Interpretabilidade |
|--------|--------------------------|-------------|----------------|-------------------|
| ARMA | Sim | Não | Não | Alta |
| ARIMA | Não (diferencia) | Não | Não | Alta |
| SARIMA | Não (diferencia) | Sim | Não | Alta |
| Holt-Winters | Não | Sim | Não | Média |
| Prophet | Não | Sim | Parcial | Média |
| LSTM | Não | Implícita | Sim | Baixa |
| XGBoost | Não | Com features | Sim | Média |

---

## 5. Fluxo de Análise

```
1. Visualizar a série → tendência, sazonalidade, outliers?
2. Decompor → aditivo ou multiplicativo?
3. Testar estacionaridade → ADF / KPSS
4. Se não-estacionária → diferenciar / transformar
5. Analisar ACF e PACF → definir p e q
6. Ajustar modelo → ARIMA, SARIMA, etc.
7. Verificar resíduos → devem ser ruído branco
8. Avaliar → MAE, RMSE, MAPE no conjunto de teste
```

---

## 6. Métricas de Avaliação

| Métrica | Fórmula | Interpretação |
|---------|---------|--------------|
| MAE | `mean(|yₜ - ŷₜ|)` | Erro médio absoluto |
| RMSE | `sqrt(mean((yₜ - ŷₜ)²))` | Penaliza erros grandes |
| MAPE | `mean(|yₜ - ŷₜ| / yₜ) × 100` | Erro percentual (evitar se yₜ ≈ 0) |
| AIC / BIC | — | Seleção de modelos (menor = melhor) |
