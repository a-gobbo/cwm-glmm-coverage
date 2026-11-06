COVERAGE OF WALD AND REsim INTERVALS IN A CLUSTER-WEIGHTED MIXED-EFFECTS LOGISTIC MODEL: THE COST OF ESTIMATING sigma_b AND THE CLUSTER ASSIGNMENT

QUESTION
In the ML-CWMd of Caldera et al. (Biometrics, 2025) the p-values of the fixed effects and the intervals of the hospital effects come from glmer fitted on the final cluster assignment, treated as known; the hospital intervals (REsim) also treat the estimated random-effect SD sigma_b as known. How much does this lower the coverage of the intervals?

Four versions are fitted on the same simulated data:
	A	true clusters, sigma_b known
	B	true clusters, sigma_b estimated
	C	CEM clusters, sigma_b known
	D	CEM clusters, sigma_b estimated (method of the paper)

Methods and results: report/report.pdf

CODE
	sim.R	data generator
	fit.R	mixed-effects logistic regression of one cluster, sigma_b estimated (glmer) or fixed; REsim intervals
	cem.R	classification EM (hard assignment, glmer M-step)
	study.R	one replicate: versions A-D, raw estimates
	run.R	8 scenarios x 100 replicates in parallel, then analyse.R
	analyse.R	coverages and Monte Carlo errors from the raw estimates
	report/	report.pdf and its figures

USAGE
	Rscript run.R	# ~10 min on 7 cores; SMOKE=1 for a quick test
	Rscript analyse.R	# metrics only, no refitting

Packages: lme4, merTools, mclust.

RESULTS (results/)
	J.._sb.._sep.._rep....rds	one replicate
	hosp_all.rds	REsim mean and sd and true effect, per hospital
	clus_all.rds	fixed effects, SE, sigma_b and assignment diagnostics, per cluster
	metrics.rds	metrics per replicate
	summary.csv	coverages (cov_b, beta_cov) with Monte Carlo errors (_se), per scenario, cluster and version
