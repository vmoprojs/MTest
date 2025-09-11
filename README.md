# MTest — Bootstrap-based Multicollinearity Diagnostic (Klein & VIF)

Functions to **detect and quantify multicollinearity** via a **nonparametric pairs bootstrap**.

`MTest` reports achieved significance levels (ASL; bootstrap proportions) for two widely used rules:

- **Klein's rule**: flag multicollinearity if $R^2_j > R^2_g$  
- **VIF rule**: flag multicollinearity if $\mathrm{VIF}_j$ is large, with $\mathrm{VIF}_j = \dfrac{1}{1 - R^2_j}$

**Reference:** Morales-Oñate & Morales-Oñate (2023). *MTest: a Bootstrap Test for Multicollinearity*. Revista Politécnica, 51(2), 53–62.  
DOI: https://doi.org/10.33333/rp.vol51n2.05

---

## What MTest does

Given a fitted linear model, `MTest`:

1. Resamples **rows** of the model frame (pairs bootstrap) `nboot` times.  
2. At each bootstrap replicate, recomputes the **global** $R^2_g$ and the **auxiliary** $R^2_j$
   (regressing each predictor on the rest), using the **same** expanded design matrix as the original fit.
   This is robust to `log()`, `I()`, interactions, factors, `poly()`, etc.  
3. Returns bootstrap distributions and **ASL** (bootstrap proportions) for:
   - **VIF rule (threshold on $R^2_j$)**:
     
$$
\mathrm{ASL}_{\mathrm{VIF}}(j) = \mathbb{P}\big(R^2_j > c\big)
$$
     
    Example: `valor_vif = 0.90` implies a VIF cutoff of $1 / (1 - 0.90) = 10$.

   - **Klein's rule**:
     
$$
\mathrm{ASL}_{\mathrm{Klein}}(j) = \mathbb{P}\big(R^2_g < R^2_j\big).
$$

These ASLs are simple **bootstrap proportions** of the corresponding events (no additional parametric assumptions).

---

## Model context

Linear regression model:

$$
Y_i = \beta_0 + \beta_1 X_{1i} + \cdots + \beta_p X_{pi} + u_i, \quad i=1,\ldots,n.
$$

Auxiliary regressions (one per predictor):

$$
X_{ji} = \gamma_0 + \sum_{k \ne j} \gamma_k X_{ki} + e_{ji}, \quad j=1,\ldots,p.
$$

Let $R^2_g$ be the global $R^2$ and $R^2_j$ the $R^2$ of the $j$-th auxiliary regression.

---

## Installation

```r
# From CRAN (if available):
# install.packages("MTest")

# From source (package root):
# R CMD build .
# R CMD INSTALL MTest_x.y.z.tar.gz
```

---

## Quick start

```r
library(MTest)

data(simDataMTest)
fit <- lm(y ~ ., data = simDataMTest)

set.seed(123)
out <- MTest(fit, nboot = 200, valor_vif = 0.90, trace = FALSE)

# ASL (bootstrap p-values)
out$pval_vif      # ASL_VIF(j)   = Pr(R^2_j > valor_vif)
out$pval_klein    # ASL_Klein(j) = Pr(R^2_g < R^2_j)

# Inspect distributions / summary
head(out$Bvals)   # nboot x (p+1): columns = "global" and one per predictor
print(out)        # compact console summary (print.MTest)
```

---

## Returned components

`MTest()` returns an object of class `"MTest"` with:

- `Bvals` — matrix `nboot x (p+1)` with bootstrap $R^2$: first column `"global"` for $R^2_g$, then one column per predictor for $R^2_j$.
- `VIFvals` — matrix `nboot x p` with bootstrap VIF per predictor.
- `pval_vif` — named vector with $\mathbb{P}(R^2_j >$ `valor_vif`$)$.
- `pval_klein` — named vector with $\mathbb{P}(R^2_g < R^2_j)$.
- `vif.tot`, `R.tot` — observed (non-bootstrap) VIF and $R^2$.
- `nsam`, `nboot` — bootstrap sample size and iterations used.

> **Note (factors):** VIF is returned **per design column**. If you need group-wise GVIFs (as in `car::vif`), compute them separately.

---

## Interpreting ASL

- **High `pval_klein[j]`** → $R^2_j$ often exceeds $R^2_g$ (Klein's rule frequently triggered).  
- **High `pval_vif[j]`** → $R^2_j$ often exceeds `valor_vif` (equivalently, $\mathrm{VIF}_j$ exceeds the implied cutoff $1/(1-c)$ ).

Do not rely on a single rule: inspect both and consider the modeling context.

---

## Pairwise KS helper

`pairwiseKStest()` runs pairwise Kolmogorov–Smirnov tests on matrix columns (e.g., the $R^2_j$ columns of `Bvals`), returning a **p-value matrix** and a simple **Suggestion** ranking that mirrors the original behavior.

```r
# Compare only predictors, exclude "global":
ks_res <- pairwiseKStest(out$Bvals[, -1],
                         alternative = "greater",  # "less" or "two.sided"
                         use = "asis",             # same behavior as original
                         exact = NULL)             # let ks.test decide

# Full p-value matrix (rows = x, cols = y)
print(ks_res)             # shows the matrix and a compact summary
ks_res$KSpwMatrix

# Suggestion field:
#   alternative = "greater" -> sort row sums (descending)
#   alternative = "less"    -> sort column sums (descending)
ks_res$Suggestion

# Object class (S3):
class(ks_res)
# [1] "pairwiseKStest" "list"
```

The package also provides a `print.pairwiseKStest()` method that always prints the p-value matrix and, optionally, a row-wise minima summary.

---

## Reproducibility tips

- Set `seed=` in `MTest()` (and/or call `set.seed()` before running).
- Use `trace = TRUE` for a progress bar on long runs.

---

## Minimal API

```r
MTest(object, nboot = 100, nsam = NULL, trace = FALSE, seed = NULL,
      valor_vif = 0.9)

pairwiseKStest(X,
               alternative = c("greater","less","two.sided"),
               use = c("asis","pairwise.complete.obs"),
               exact = NULL)
```

- `valor_vif` is a **threshold on $R^2_j$**. The implied VIF cutoff is $1/(1-$`valor_vif`$)$ (e.g., `0.90` ↔ `10`).  
- `pairwiseKStest()` expects a numeric matrix/data frame; for `Bvals` you typically pass `Bvals[, -1]` to exclude `"global"`.

---

## Citation

> Morales-Oñate, V., & Morales-Oñate, B. (2023).  
> *MTest: a Bootstrap Test for Multicollinearity*. Revista Politécnica, 51(2), 53–62.  
> https://doi.org/10.33333/rp.vol51n2.05

---

## License

MIT (or your package license). Include the corresponding `LICENSE` file in the repo.

---

## Acknowledgements

Thanks to the R community and the maintainers of `stats`, `car`, and related tooling.

---

## Changelog (from earlier README)

- Clarified that ASL are **bootstrap proportions** for Klein and VIF events.  
- Documented the **pairs bootstrap** scheme and robustness to transformed terms (`log()`, `I()`, interactions, factors, `poly()`).  
- Added the **pairwiseKStest** helper with class `c("pairwiseKStest","list")` and a `print()` method that prints the full p-value matrix.
