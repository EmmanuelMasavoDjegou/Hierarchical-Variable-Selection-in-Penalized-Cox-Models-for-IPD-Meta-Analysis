# Response to the Associate Editor

**Manuscript:** Hierarchical Variable Selection for Nonconvex Penalized Cox Models in Individual-Participant-Data Meta-Analysis (previously CSDA-D-26-01273)

We thank the Associate Editor for a careful and constructive assessment. We took every comment as a request for substantive work rather than rewording. The revision is a new version of the paper: the estimator, algorithm, software, simulation study and data analysis have all been redone. In the course of this work we found, and corrected, several problems in the first version that the comments did not mention; we list them first because they change the results.

---

## Problems found in the first version and corrected

1. **The code did not fit the model described in the paper.** The global step used an *unstratified* penalized Cox fit (`ncvreg::ncvsurv` on stacked data, risk sets mixing studies). This contradicts the stratified partial likelihood that the paper describes and that motivates the method. The new solver computes all risk sets within study (§3.1). Its unpenalized fit reproduces `survival::coxph(... + strata(study))` to numerical precision, and this is checked by the package tests.
2. **The zero-sum constraint is incompatible with sparse deviations.** With $\sum_k \varepsilon_{kj}=0$, a single deviating study forces every other study to deviate. The constraint was also not implemented in the code. We replace it with a *sharing rule*: every active covariate has at least one study with no deviation, and among admissible decompositions the one with the smallest penalty is chosen (§2.2, eq. 2.2). Identification in the theory now rests on an explicit assumption (Assumption 1(ii)), and the proofs have been rewritten for this parameterization.
3. **The algorithm described was not the one used.** The paper described local quadratic approximation, while the code alternated `ncvreg` calls. §3.1 and Algorithm 1 now describe exactly what the C code does. Proposition 4.1 is restricted to what the algorithm actually guarantees: feasibility and monotone decrease of the objective. We removed the claim of convergence to the oracle solution.
4. **Tuning was sequential, not joint.** The paper stated a joint two-dimensional grid; the code tuned $\lambda_\alpha$ and $\lambda_\varepsilon$ one after the other, using different criteria for different methods. §3.2 now defines one joint procedure, applied identically to every method.

As a consequence, **the conclusions changed.**
- The first version attributed the lower false discovery rate to the hierarchical decomposition. The new simulation, which adds a stratified common-effect competitor, shows that stratification and the nonconvex penalty drive global selection.
- The hierarchical deviations mainly buy much better estimation of study-specific effects (1.2 to 7.7 times lower squared error when heterogeneity is present, and no loss when it is absent) and identification of the deviating studies. With the one-standard-error rule, they also give the best global selection in most scenarios.
- In the data, *IGF2*, previously reported as the strongest global biomarker, turns out to have an essentially null global effect. Its effect is concentrated in one cohort.

We report these findings as they are.

---

## Point-by-point responses

**1. "The research results reported are too premature … the descriptions on the methodology, simulation study, literature review, and the data example are not written well."**

Every section has been rewritten; see the items below and the list of corrections above.

**2. "A more detailed method for choosing two tuning parameters should be described. The objective function should be defined clearly, and the suggested number of folds should be given."**

- *Objective function.* Eq. (2.6) now states it in full: the negative stratified partial log-likelihood divided by $N$, an SCAD/MCP penalty on the global effects with $\lambda_\alpha$, and a penalty on the deviations with $\lambda_\varepsilon$ and study weights $\omega_k=(n_k/N)^{1/2}$. It is minimized over the set $\mathcal H$ defined by the sharing and hierarchy rules (eqs. 2.2–2.3). The role of each term and of $\omega_k$ is explained.
- *Tuning (§3.2).* The procedure has four steps:
  - Grid: 20 log-spaced values of $\lambda_\alpha$ from $\lambda_{\max}$ to $0.05\lambda_{\max}$, times the ratios $r=\lambda_\varepsilon/\lambda_\alpha\in\{0.5,1,2\}$.
  - Folds: stratified by study and by event status, so every training set contains every study.
  - Criterion: the cross-validated linear-predictor deviance of the stratified partial likelihood (eq. 3.2; Dai & Breheny, 2019).
  - Selection: a minimum rule and a one-standard-error rule. The latter is defined precisely, as the fewest nonzero parameters within one standard error.
- *Number of folds.* We recommend $V=5$ and use it throughout. We explain when $V=10$ is preferable, and §6 reports how the selected genes change with $V=10$. Figure 6 shows the cross-validation curves for the data.

**3. "The proposed method is computationally complex … Without the development of a software package (e.g., an R package), the method is not used by readers."**

We developed the R package **`hiermetacox`** (§3.4), with its solver written in C.
- *Functions:* `hmcox()` fits the path, `cv.hmcox()` performs the two-dimensional cross-validation and returns both rules from one run, `hmcox_inference()` gives hazard ratios, CIs and p-values, `hmcox_boot()` gives selection stability, and `sim_ipd()` generates data from the simulation design.
- *Speed:* a complete five-fold, three-ratio cross-validation takes about 4 seconds for $N=1{,}000$ and $p=100$, and about 18 seconds for the 1,446 patients and 500 genes of the application (Table 5).
- *Tests:* the test suite verifies the solver against exact references, namely `coxph` with strata, separate per-study Cox fits, `glmnet`'s stratified lasso path, and the KKT conditions.

**4. "The code for the simulation and data analysis should be made available … in any unrestricted form."**

Everything is public at https://github.com/EmmanuelMasavoDjegou/Hierarchical-Variable-Selection-in-Penalized-Cox-Models-for-IPD-Meta-Analysis.
- *Contents:* the package, the simulation and application scripts, all per-replicate simulation results, and a consistency-check script.
- *No manual transcription:* every table in the paper is `\input` from a file written by the scripts, and every figure is a script output. The README maps each table and figure to the script that produces it.
- *Data:* the data are public (Bioconductor `curatedOvarianData`).
- A "Code and data availability" statement has been added.

**5. "The base setting of p = 50 is not really high-dimensional. This could be at least 100."**

- *Simulation:* the base setting is now $p=100$. Other scenarios use $p=200$, $p=500$ (with $p$ exceeding every study's sample size and equal to $N/2$), and $n_k=100$ with $p=200$.
- *Application:* it now uses $p=500$ genes, chosen by an outcome-blind filter, instead of 66. This exceeds the sample size of four of the six cohorts.
- *Wording:* we describe the setting precisely rather than as generically "high-dimensional". The theory covers $p=O(N^c)$ with $c<1$, and we state that $p\gg N$ requires screening.

**6. "Section 5.2 is difficult to read. Please avoid acronyms, such as S01–S03."**

- Scenarios are now named descriptively ("Strong correlation", "Heavy censoring", "Small studies, p > n_k") and listed in Table 1 with the single feature changed from the base configuration.
- Methods are named in words, e.g. "MCP, hierarchical (proposed)" and "MCP, stratified common effects".
- The figures display scenario names on the vertical axis.

**7. "The introduction … should focus more on … survival analysis methods for IPD meta-analytic data … Nothing is reviewed on these random-effect meta-analysis methods … a relevant work is [1] Utazirubanda et al. (2021)."**

- The introduction has been rewritten around three topics: one- and two-stage IPD meta-analysis of time-to-event outcomes; handling of baseline and effect heterogeneity by stratification, shared frailties and random slopes; and variable selection in Cox and frailty models. It now cites Tudur Smith et al. 2005, Bowden et al. 2011, Burke et al. 2017, Debray et al. 2015, Vaida & Xu 2000, Ripatti & Palmgren 2000, Therneau et al. 2003, Ha et al. 2001, Fan & Li 2002, Androulakis et al. 2012, Ha et al. 2014, Groll et al. 2017 and Utazirubanda et al. 2021.
- It also reviews fixed-effect multi-study selection methods closest to ours (Ma et al. 2011; He et al. 2016; Cheng et al. 2015; Gross & Tibshirani 2016; Ollier & Viallon 2017).
- The textbook material on ridge, elastic net and LARS has been removed, as have references unrelated to survival analysis.
- A new subsection (§2.3) explains how our model relates to frailty and random-slope models and when each is preferable.
- The data analysis compares stratified, shared gamma frailty and pooled hazard ratios for the selected genes (Table 7).

**8. "Real data: The gene selection results do not have a clear conclusion. Table 4 shows the log-HR without the SE, confidence interval, or P-value … medical researchers may look at the HR."**

- Table 6 now reports hazard ratios per within-cohort SD, 95% confidence intervals, Wald p-values, Holm-adjusted p-values, the number of deviating cohorts, bootstrap selection frequencies, and how many of the eight methods select each gene. Figure 7 shows global and study-specific hazard ratios with CIs.
- We state what these intervals are: post-selection intervals, conditional on the selected model. Their empirical coverage is evaluated in the simulation (Table 3), where they under-cover. This is why we report them alongside bootstrap stability.
- §6.4 states the conclusion explicitly:
  - GMPR and CXCL9 are protective and robust.
  - IGF2's association is confined to one cohort.
  - The signal beyond these genes is modest and tuning-dependent.
  - Dense convex models predict a new cohort slightly better than the sparse models (leave-one-study-out validation, Table 8).

**9. "Below (2.3): … readers are not clear if α … and ε_k … are known or unknown and fixed or random."**

§2.2 now lists every component of the model. The baseline hazards $h_{0k}$ are unknown, unspecified and not estimated. $\boldsymbol\alpha$ and $\boldsymbol\varepsilon_k$ are *unknown fixed* parameters: not random effects, with no distribution assumed. The meaning of $\varepsilon_{kj}=0$ is also stated there, and §2.3 contrasts this with random-effect formulations.

---

We believe the revised manuscript addresses every concern of the Associate Editor, and we would be grateful for the opportunity to submit it.
