% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_formulas.R
\name{BCRA3_2r}
\alias{BCRA3_2r}
\title{Model 4.2: MDES Calculator for 3-Level Random Effects Blocked Cluster Random Assignment Designs (BCRA3_2r)—Treatment at Level 2}
\usage{
BCRA3_2r(R2_1, R2_2, R2_T3, P, g_3, n, J, K, rho_2, rho_3, omega_3,
  alpha = 0.05, power = 0.8, tails = 2)
}
\arguments{
\item{R2_1}{Proportion of variance in level 1 outcome explained by level 1 covariates}

\item{R2_2}{Proportion of variance in level 2 outcome explained by level 2 covariates}

\item{R2_T3}{Proportion of between block variance in treatment effect explained by Level-3 covariates}

\item{P}{Proportion of the level 3 units randomized to the treatment}

\item{g_3}{Number of level 3 Covariates used}

\item{n}{Mean # of level 1 units per level 2 unit}

\item{J}{Mean # of level 2 units per level 3 unit}

\item{K}{Number of level 3 units}

\item{rho_2}{proportion of variance in outcome that is between level 2 units}

\item{rho_3}{proportion of variance in outcome that is between level 3 units}

\item{omega_3}{Level 3 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by Level-3 outcome variation}

\item{alpha}{Probability of a type I error}

\item{power}{Statistical power (1 - probability of type II error)}

\item{tails}{One or two tailed test}
}
\value{
Minimum Detectable Effect Size
}
\description{
MDEs for 3-Level Random Effects Blocked Cluster Random Assignment Designs
}

