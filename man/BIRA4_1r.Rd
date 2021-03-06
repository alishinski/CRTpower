% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_formulas.R
\name{BIRA4_1r}
\alias{BIRA4_1r}
\title{Model 2.5:  MDES Calculator for 4-Level Random Effects Blocked Individual Random Assignment Design (BIRA4_1r)—Individuals Randomized within Blocks}
\usage{
BIRA4_1r(R2_1, R2_2T, R2_3T, R2_4T, P, g_4, n, J, K, L, rho_4, rho_3, rho_2,
  omega_2, omega_3, omega_4, alpha = 0.05, power = 0.8, tails = 2)
}
\arguments{
\item{R2_1}{Proportion of variance in level 1 outcome explained by level 1 covariates}

\item{R2_2T}{Proportion of between block variance in treatment effect explained by Level 2 covariates}

\item{R2_3T}{Proportion of variance in treatment effect explained by Block and Level 3 covariates}

\item{R2_4T}{Proportion of variance in treatment effect explained by  Block and Level 4 covariates}

\item{P}{Proportion of level 1 units randomized to the treatment}

\item{g_4}{Number of level 4 covariates}

\item{n}{Mean # of level 1 units per level 2 unit}

\item{J}{Mean # of level 2 units per level 3 unit}

\item{K}{Number of level 3 units per level 4 unit}

\item{L}{Number of level 4 units}

\item{rho_4}{proportion of variance in outcome between level 4 units}

\item{rho_3}{proportion of variance in outcome between level 3 units}

\item{rho_2}{proportion of variance in outcome between level 2 units}

\item{omega_2}{Level 2 treatment effect heterogeneity: variance in treatment effect across Level 2 units, standardized by the Level-2 outcome variation}

\item{omega_3}{Level 3 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by the Level-3 outcome variation}

\item{omega_4}{Level 4 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by the Level-4 outcome variation}

\item{alpha}{Probability of a type I error}

\item{power}{Statistical power (1 - probability of type II error)}

\item{tails}{One or two tailed test}
}
\value{
Minimum Detectable Effect Size
}
\description{
MDES for 4-Level Random Effects Blocked Individual Random Assignment Design
}

