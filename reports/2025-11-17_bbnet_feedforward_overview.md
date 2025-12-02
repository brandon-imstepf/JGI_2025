# CallGenes Feed-Forward Network Overview (2025-11-17)

## Context
- CallGenes multiplies heuristic ORF scores by a BBNet/CellNet feed-forward network before the DP trace so the NN remains advisory (`INIT.md:3-6`).
- Training data uses 356 features per ORF (8 scalar + 348 one-hot bases) with 80/10/10 genome splits drawn from NCBI CDS positives vs. “all genes” negatives (`INIT.md:3-6`).
- Current production jobs (e.g., `slurm/cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_orig_tol3.err:1-6`) instantiate `ml.Trainer` with six dense layers `[358,125,275,75,25,1]`, 400k batches, and four evolutionary cycles.

## Architecture & Implementation
- **Network representation** – `ml.CellNet` builds dense or sparse multilayer perceptrons with configurable per-layer widths, densities, and edge block sizes. Dense mode is default (`bbmap/current/ml/CellNet.java:1118-1125`) and each `Cell` stores its bias, weights, and cached activation/error buffers (`bbmap/current/ml/Cell.java:595-620`).
- **Activations** – Each `Cell` wraps a `Function` instance; defaults are sigmoid hidden layers (`Function.SIG`) and an RSLOG (logistic-style) output layer (`bbmap/current/ml/Cell.java:629-636`, `bbmap/current/ml/Function.java:1-74`). Activation types can be randomized per layer via `Function.TYPE_RATES`.
- **Forward/Backward passes** – Dense inference follows `CellNet.feedForwardDense()` while gradients accumulate per cell so `applyUpdates()` can update weights with per-edge annealing safeguards (edge amplification throttle at 0.98, `bbmap/current/ml/Cell.java:102-120`).
- **Sample management** – `ml.Sample` stores input vectors, goals, cached predictions, and per-sample weights/balancing metadata (pivot, epoch). `ml.Trainer` partitions data into subsets and orchestrates worker threads plus evolutionary “cycles” to keep multiple candidate networks in flight (`bbmap/current/ml/Trainer.java:350-430`, `bbmap/current/ml/Trainer.java:1080-1167`).

### Activation Functions (Latex-Ready)
All eight feed-forward activations live in `ml.Function`/`ml.Functions`. Sigmoid (hidden layers) and RSLOG (output) are the current defaults, but the others are available via type rates or explicit layer configs.

| Function | Formula |
|----------|---------|
| Sigmoid (`SIG`) | \( \sigma(x) = \dfrac{1}{1 + e^{-x}} \) |
| Extended Sigmoid (`ESIG`) | \( \operatorname{eSig}(x) = 2\,\sigma(x) - 1 \) |
| Hyperbolic Tangent (`TANH`) | \( \tanh(x) = \dfrac{e^{x} - e^{-x}}{e^{x} + e^{-x}} \) |
| Rotationally Symmetric Log (`RSLOG`) | \( \operatorname{rslog}(x) = \operatorname{sign}(x)\,\ln\!\big(|x| + 1\big) \) |
| Swish (`SWISH`) | \( \operatorname{swish}(x) = x\,\sigma(x) \) |
| Mirrored Sigmoid (`MSIG`) | \( \operatorname{mSig}(x) = \frac{1}{\sigma(5)} \times \begin{cases} \sigma(2x + 5), & x < 0 \\ \sigma(-2x + 5), & x \ge 0 \end{cases} \) |
| Extended Mirrored Sigmoid (`EMSIG`) | \( \operatorname{emSig}(x) = 2\,\operatorname{mSig}(x) - 1 \) |
| Bell / Gaussian (`BELL`) | \( \operatorname{bell}(x) = e^{-x^{2}} \) |

Derivatives used in `ml.Cell` follow directly from these definitions (e.g., \( \sigma'(x) = \sigma(x)\big(1-\sigma(x)\big) \), \( \tanh'(x) = 1-\tanh^{2}(x) \), \( \frac{d}{dx}\operatorname{rslog}(x) = \frac{1}{|x|+1} \), etc.), matching the helpers in `bbmap/current/ml/Functions.java:1-220`.

**When activations fire:** Each `Cell` in layer \( \ell \) first computes a weighted sum of layer \( \ell-1 \) outputs plus its bias (`Cell.summateDense/summateSparse`). The function \( f_{\ell}(\cdot) \) is immediately applied to that sum to produce the layer output; those values feed the next layer during the same forward pass. Hidden layers all use `Cell.defaultActivationType` (sigmoid unless overridden), and the final layer applies `Cell.finalLayerType` (RSLOG) just before the network emits its probability-like score.

## Training Workflow Highlights
- Command-line knobs (`maxdims`, `mindims`, `density`, etc.) are parsed in `ml.Trainer.processArgs`, enabling Slurm scripts to pin layer widths or let the trainer randomize architectures per seed (`bbmap/current/ml/Trainer.java:350-389`).
- `Trainer.randomNetwork` seeds `CellNet` instances, sets densities derived from `maxEdges` or `maxDensity`, and invokes `CellNet.randomize()` (which fills edges & activations) before handing them to worker threads (`bbmap/current/ml/Trainer.java:1099-1158`, `bbmap/current/ml/CellNet.java:134-214`).
- Training mixes 64 worker threads with 8 trainer threads (see Slurm log) and uses subset triage to emphasize difficult samples over time (see `ml/Subset.java`, not detailed here).

## Loss Functions in Production
Only the BBNet/CellNet stack is in scope (CNN losses are ignored per instructions). Production runs minimize a squared-error objective with classification-aware weighting:

1. **Base squared-error loss** – Every output neuron computes `0.5 * (goal - prediction)^2` (`bbmap/current/ml/Cell.java:52-61`, `bbmap/current/ml/Sample.java:73-105`). This is the sole scalar loss accumulated when Trainer reports `error`.
2. **Weighted classification penalties** – Before backpropagating, `Cell.calcETotalOverOut` scales the gradient `v - target` by multipliers that depend on whether the sample is a false positive, false negative, or an “excess” error on the correct side of the cutoff (`bbmap/current/ml/Cell.java:520-557`). Default multipliers emphasize misclassification near the decision boundary:  
   - False positives: `falsePositiveErrorMult = 10.5`  
   - False negatives: `falseNegativeErrorMult = 10.5`  
   - Correct-side but oversized errors (“excess”): `0.2`  
   - Background positive/negative errors: `1.0` (`bbmap/current/ml/Cell.java:649-657`).  
   These constants, plus a 0.05 “spread” window around the 0.5 cutoff, are what production scripts rely on today.

In short, the BBNet trainer optimizes a weighted mean-squared-error objective tailored to binary ORF classification—no additional loss functions are wired into production feed-forward runs.

**When the loss fires:** After the forward pass produces the RSLOG output for a sample, `ml.Cell.calcError`/`ml.Sample.calcError` compute \( \tfrac{1}{2}(y - \hat{y})^2 \). That scalar loss immediately drives `Cell.calcETotalOverOut`, which multiplies the raw gradient `(prediction - goal)` by the FP/FN weights before the trainer performs backpropagation through the layers. The loss is therefore evaluated once per sample per batch right at the network output, and its weighted gradient is what flows backward to update all hidden activations.

## Slide/Paper Sound Bites
- “BBNet layers are hand-built dense perceptrons (sigmoid hidden layers, RSLOG output) with sparse-aware math primitives; we train 32 candidate nets per Slurm job via evolutionary cycling.”
- “Inputs are 356-dimensional ORF descriptors, and the current stack trains for ~400k batches with aggressive FP/FN penalties baked directly into the squared-error gradient.”
- “Loss = 0.5 (goal − prediction)², scaled x10.5 for boundary FPs/FNs, so the NN learns to preserve recall while tamping down high-score false positives before CallGenes applies its DP trace.” 
