# self-test-epidemic
Deterministic SEIR modelling code evaluating self-testing for epidemic-prone pathogens. Repository includes global sensitivity analysis (20,000 LHS), SROC-based diagnostic accuracy sampling, and modular helper functions to generate simulation outputs and quantify mitigation and marginal efficiency under varying compliance and transmission regimes.

**Self-testing SEIR model: Deterministic transmission model with diagnostic uncertainty and global sensitivity analysis**

_Overview_

This repository contains R code implementing a deterministic SEIR transmission model incorporating self-testing and isolation dynamics. 
The model evaluates the epidemiological impact of self-testing under varying behavioural, diagnostic, and transmission conditions. Diagnostic accuracy uncertainty is incorporated using paired sensitivity and specificity draws from a bivariate random-effects SRoC (reitsma) model, truncated at 0.90 specificity to prevent epidemiological impact estimates being artificially inflated under implausibly low specificity. Global Sensitivity Analysis (GSA) is conducted using Latin Hypercube Sampling (n = 20,000) and Partial Rank Correlation Coefficients (PRCC).
The repository provides full transparency of model structure, parameter sampling, and outcome derivation used in the associated manuscript.

_Model structure_

The transmission model extends a standard SEIR framework by including:
	Explicit testing pathways (standard of care and self-testing)
	Behavioural compliance to isolation
	Diagnostic sensitivity and specificity influencing case detection
	Epidemic endpoints including peak infectious prevalence, relative peak reduction cumulative deaths and deaths averted under intervention.
The system is solved using deSolve::ode.

_Model parameters_

The following parameters are sampled within the GSA framework unless otherwise stated:
	R0: Basic reproduction number 
	CFR: Case fatality ratio
	epsilon: latent period 
	D: Duration of pathogen infectiousness
	sigma: self-test flowrate (rate at which infectious individuals access self-testing)
	sens: diagnostic sensitivity
	spec: diagnostic specificity
	c: hazard of adherence to post-test isolation conditional on positive test; influences transmission dynamics in two ways. First, isolated susceptibles may breach isolation, becoming exposed to infection at rate (1-c)λS_i. Second, infectious individuals who are isolating transmit at a reduced effective contact rate proportional to their breach frequency, contributing (1-c)βI_i/N to the overall force of infection (represented as β_i=β(1-c)). 
	r: health system case detection hazard (standard-of-care testing pathway)
Parameter ranges reflect plausible outbreak conditions derived from published literature and continental preparedness frameworks. Ranges are specified in the GSA script.

_Diagnostic accuracy modelling_

Diagnostic uncertainty is incorporated via a bivariate random-effects SRoC model (reitsma approach). 
Study-level 2×2 data are stored in data/diagnostic.csv (derived from systematic reviews reported in the associated manuscripts supplementary material) and must include:
	TP (true positives)
	FP (false positives)
	FN (false negatives)
	TN (true negatives)
No individual-level data are included.
Paired posterior draws are generated and mapped directly onto the LHS parameter matrix, preserving correlation between sensitivity and specificity.
Specificity draws below 0.90 are excluded prior to sampling.

_Global Sensitivity Analysis_

Global uncertainty is explored using Latin Hypercube Sampling (lhs), n = 20,000.
Sensitivity and specificity are not independently sampled. instead, lhs-generated placeholders are replaced by sroc-derived paired draws.
Model outputs are evaluated using Partial Spearman Rank Correlations (prcc).
Random seeds are set to ensure reproducibility.

_Methodological notes_

The model is deterministic and does not incorporate stochastic fade-out.
Parameter ranges reflect plausible values drawn from outbreak literature and continental policy guidance.
Specificity truncation modifies the implied prior distribution of diagnostic performance and is explicitly enforced in sampling.
Simulations are seeded to ensure reproducibility.

_License_

This repository is released under the MIT license.

_Citation_

If you use this code, please cite:
Dunkley Y. Self-testing SEIR model with diagnostic uncertainty and global sensitivity analysis. GitHub repository. 2026. Available at: https://github.com/yasmindunkley/self-test-epidemic 
A versioned DOI will be added upon formal release.


