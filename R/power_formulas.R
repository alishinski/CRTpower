#' @title Model 1.0: Individual Random assignment
#' @description MDES for Individual Random assignment
#'
#' @param R2 Proportion of variance explained in the outcome by the covariates
#' @param P Proportion of the sample randomized to the treatment
#' @param k Number of Covariates used
#' @param n Total Sample Size
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
IRA <- function(R2, P, k, n, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), n - k - 2), stats::qt(1 - alpha, n - k - 2))
  T2 <- ifelse(power >= 0.50, stats::qt(power, n - k - 2), stats::qt(1 - power, n - k - 2))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes = M * sqrt((1 - R2) / (P * (1 - P) * n))
  mdes
}

#' @title Model 2.1:  MDES Calculator for 2-Level Constant Effects Blocked Individual Random Assignment Designs (BIRA2_1c)—School Intercepts Only
#' @description MDES for 2-Level Constant Effects Blocked Individual Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by block and level 1 covariates
#' @param P Proportion of the level 1 units randomized to the treatment
#' @param g Number of level 1 Covariates used
#' @param n Average Block Size (mean # of level 1 units per level 2 unit)
#' @param J Number of level 2 units
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BIRA2_1c <- function(R2_1, P, g, n, J, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), n*J-J-g-1), stats::qt(1 - alpha, n*J-J-g-1))
  T2 <- ifelse(power>=0.5, stats::qt(power, n*J-J-g-1), stats::qt(1-power,n*J-J-g-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes = M * sqrt((1-R2_1)/(P * (1-P) * n * J))
  mdes
}

#' @title Model 2.2:  MDES Calculator for 2-Level Fixed Effects Blocked Individual Random Assignment Designs (BIRA2_1f)— School Intercepts and Interactions with TREATMENT
#' @description MDES for 2-Level Fixed Effects Blocked Individual Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by block and level 1 covariates
#' @param P Proportion of level 1 units randomized to the treatment
#' @param g Number of level 1 covariates
#' @param n Average Block Size (mean # of level 1 units per level 2 unit)
#' @param J Number of level 2 units
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BIRA2_1f <- function(R2_1, P, g, n, J, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), n*J-2*J-g-1), stats::qt(1 - alpha, n*J-2*J-g-1))
  T2 <- ifelse(power>=0.5, stats::qt(power, n*J-2*J-g-1), stats::qt(1-power,n*J-2*J-g-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  mdes = M * sqrt((1-R2_1)/(P * (1-P) * n * J))
  mdes
}



#' @title Model 2.3:  MDES Calculator for 2-Level Random Effects Blocked Individual Random Assignment (BIRA2_1r) Designs— Individuals Randomized within Blocks
#' @description MDES for 2-Level Random Effects Blocked Individual Random Assignment
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2T Proportion of between block variance in treatment effect explained by Level 2 covariates
#' @param P Proportion of level 1 units randomized to the treatment
#' @param g Number of level 2 covariates
#' @param n Average Block Size (mean # of level 1 units per level 2 unit)
#' @param J Number of level 2 units
#' @param rho proportion of variance in outcome between clusters
#' @param omega Treatment effect heterogeneity:  variability in treatment effects across Level 2 units, standardized by the variability in the Level-2 outcome
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BIRA2_1r <- function(R2_1, R2_2T, P, g, n, J, rho, omega, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), J - g - 1), stats::qt(1 - alpha, J - g - 1))
  T2 <- ifelse(power>=0.5, stats::qt(power, J-g-1), stats::qt(1-power,J-g-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes = M * sqrt(rho * omega * (1-R2_2T)/ J + (1 - rho) * (1 - R2_1)/(P*(1-P)*J*n))
  mdes
}


#' @title Model 2.4:  MDES Calculator for 3-Level Random Effects Blocked Individual Random Assignment Design (BIRA3_1r)—Individuals Randomized within Blocks
#' @description MDES for 3-Level Random Effects Blocked Individual Random Assignment Design
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2T Proportion of between block variance in treatment effect explained by Level 2 covariates
#' @param R2_3T Proportion of between block variance in treatment effect explained by Level 3 covariates
#' @param P Proportion of level 1 units randomized to the treatment
#' @param g_3 Number of level 3 covariates
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units
#' @param rho_3 proportion of variance in outcome between level 3 units
#' @param rho_2 proportion of variance in outcome between level 2 units
#' @param omega_3 Level 3 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by the Level-3 outcome variation
#' @param omega_2 Level 2 treatment effect heterogeneity: variance in treatment effect across Level 2 units, standardized by the Level-2 outcome variation
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BIRA3_1r <- function(R2_1, R2_2T, R2_3T, P, g_3, n, J, K, rho_3, rho_2, omega_2, omega_3, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), K - g_3 - 1), stats::qt(1 - alpha, K - g_3 - 1))
  T2 <- ifelse(power>=0.5, stats::qt(power, K-g_3-1), stats::qt(1-power,K-g_3-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_3*omega_3*(1-R2_3T)/K+rho_2*omega_2*(1-R2_2T)/(J*K)+(1-rho_2-rho_3)*(1-R2_1)/(P*(1-P)*J*n*K))
  mdes
}

#' @title Model 2.5:  MDES Calculator for 4-Level Random Effects Blocked Individual Random Assignment Design (BIRA4_1r)—Individuals Randomized within Blocks
#' @description MDES for 4-Level Random Effects Blocked Individual Random Assignment Design
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2T Proportion of between block variance in treatment effect explained by Level 2 covariates
#' @param R2_3T Proportion of variance in treatment effect explained by Block and Level 3 covariates
#' @param R2_4T Proportion of variance in treatment effect explained by  Block and Level 4 covariates
#' @param P Proportion of level 1 units randomized to the treatment
#' @param g_4 Number of level 4 covariates
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units per level 4 unit
#' @param L Number of level 4 units
#' @param rho_4 proportion of variance in outcome between level 4 units
#' @param rho_3 proportion of variance in outcome between level 3 units
#' @param rho_2 proportion of variance in outcome between level 2 units
#' @param omega_4 Level 4 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by the Level-4 outcome variation
#' @param omega_3 Level 3 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by the Level-3 outcome variation
#' @param omega_2 Level 2 treatment effect heterogeneity: variance in treatment effect across Level 2 units, standardized by the Level-2 outcome variation
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BIRA4_1r <- function(R2_1, R2_2T, R2_3T, R2_4T, P, g_4, n, J, K, L, rho_4, rho_3, rho_2, omega_2, omega_3, omega_4, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), L - g_4 - 1), stats::qt(1 - alpha, L - g_4 - 1))
  T2 <- ifelse(power>=0.5, stats::qt(power, L-g_4-1), stats::qt(1-power,L-g_4-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_4*omega_4*(1-R2_4T)/L+rho_3*omega_3*(1-R2_3T)/(L*K)+rho_2*omega_2*(1-R2_2T)/(J*K*L)+(1-rho_2-rho_3-rho_4)*(1-R2_1)/(P*(1-P)*J*n*K*L))
  mdes
}

#' @title Model 3.1:  MDES Calculator for Two-Level Cluster Random Assignment Design (CRA2_2)— Treatment at Level 2
#' @description MDEs for Two-Level Cluster Random Assignment Design
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2T Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param P Proportion of the level 2 units randomized to the treatment
#' @param g Number of level 2 Covariates used
#' @param n Mean # of level 1 units per level 2 cluster
#' @param J Number of level 2 units
#' @param rho proportion of variance in outcome that is between clusters
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
CRA2_2r <- function(R2_1, R2_2T, P, g, n, J, rho, alpha = 0.05, power = 0.80, tails = 2){
    T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), J - g - 2), stats::qt(1 - alpha, J - g - 2))
    T2 <- ifelse(power>=0.5, stats::qt(power, J-g-2), stats::qt(1-power,J-g-2))
    M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

    mdes <- M*sqrt(rho*(1-R2_2T)/(P*(1-P)*J)+(1-rho)*(1-R2_1)/(P*(1-P)*J*n))
    mdes
}


#' @title Model 3.2: MDES Calculator for 3 Level Cluster Random Assignment Designs (CRA3_3r) Treatment at Level 3
#' @description MDES for 3 Level Cluster Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param R2_3 Proportion of variance in level 3 outcome explained by level 3 covariates
#' @param P Proportion of the level 3 units randomized to the treatment
#' @param g_3 Number of level 3 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param rho_3 proportion of variance in outcome that is between level 3 units
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
CRA3_3r <- function(R2_1, R2_2, R2_3, P, g_3, n, J, K, rho_2, rho_3, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), K-g_3-2), stats::qt(1 - alpha, K-g_3-2))
  T2 <- ifelse(power>=0.5, stats::qt(power, K-g_3-2), stats::qt(1-power,K-g_3-2))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_3*(1-R2_3)/(P*(1-P)*K)+rho_2*(1-R2_2)/(P*(1-P)*J*K)+(1-rho_2-rho_3)*(1-R2_1)/(P*(1-P)*J*n*K))
  mdes
}

#' @title Model 3.3:  MDES Calculator for 4-Level Cluster Random Assignment Designs (CRA4_4r)—Treatment at Level 4
#' @description MDES for 4-Level Cluster Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param R2_3 Proportion of variance in level 3 outcome explained by level 3 covariates
#' @param R2_4 Proportion of variance in level 4 outcome explained by level 4 covariates
#' @param P Proportion of the level 4 units randomized to the treatment
#' @param g_4 Number of level 4 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units per level 4 unit
#' @param L Number of level 4 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param rho_3 proportion of variance in outcome that is between level 3 units
#' @param rho_4 proportion of variance in outcome that is between level 4 units
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
CRA4_4r <- function(R2_1, R2_2, R2_3, R2_4, P, g_4, n, J, K, L, rho_2, rho_3, rho_4, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), L-g_4-2), stats::qt(1 - alpha, L-g_4-2))
  T2 <- ifelse(power>=0.5, stats::qt(power, L-g_4-2), stats::qt(1-power,L-g_4-2))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_4*(1-R2_4)/(P*(1-P)*L)+rho_3*(1-R2_3)/(P*(1-P)*L*K)+rho_2*(1-R2_2)/(P*(1-P)*J*K*L)+(1-rho_2-rho_3-rho_4)*(1-R2_1)/(P*(1-P)*J*n*K*L))
  mdes
}

#'  @title Model 4.1:  MDES Calculator for 3-Level Fixed Effects Blocked Cluster Random Assignment Design (BCRA3_2f)—Treatment at Level 2
#'  @description MDES for 3-Level Fixed Effects Blocked Cluster Random Assignment Design
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param P Proportion of the level 2 units randomized to the treatment
#' @param g_2 Number of level 2 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BCRA3_2f <- function(R2_1, R2_2, P, g_2, n, J, K, rho_2, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), K*(J-2)-g_2), stats::qt(1 - alpha, K*(J-2)-g_2))
  T2 <- ifelse(power>=0.5, stats::qt(power, K*(J-2)-g_2), stats::qt(1-power,K*(J-2)-g_2))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_2*(1-R2_2)/(P*(1-P)*J*K)+(1-rho_2)*(1-R2_1)/(P*(1-P)*J*n*K))
  mdes
}

#' @title Model 4.2: MDES Calculator for 3-Level Random Effects Blocked Cluster Random Assignment Designs (BCRA3_2r)—Treatment at Level 2
#' @description MDEs for 3-Level Random Effects Blocked Cluster Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param R2_T3 Proportion of between block variance in treatment effect explained by Level-3 covariates
#' @param P Proportion of the level 3 units randomized to the treatment
#' @param g_3 Number of level 3 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param rho_3 proportion of variance in outcome that is between level 3 units
#' @param omega_3 Level 3 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by Level-3 outcome variation
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BCRA3_2r <- function(R2_1, R2_2, R2_T3, P, g_3, n, J, K, rho_2, rho_3, omega_3, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), K-g_3-1), stats::qt(1 - alpha, K-g_3-1))
  T2 <- ifelse(power>=0.5, stats::qt(power, K-g_3-1), stats::qt(1-power,K-g_3-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_3*omega_3*(1-R2_T3)/K+rho_2*(1-R2_2)/(P*(1-P)*J*K)+(1-rho_2-rho_3)*(1-R2_1)/(P*(1-P)*J*n*K))
  mdes
}

#' @title Model 4.3:  MDES Calculator for 4-Level Random Effects Block Random Assignment Designs (BCRA4_2r)—Treatment at Level 2
#' @description MDES for 4-Level Random Effects Block Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param R2_3T Proportion of between block variance in impacts on outcome between Level-3 blocks explained by Level-3 covariates
#' @param R2_4T Proportion of variance in impacts on outcomes between Level-4 blocks explained by Level-4 covariates
#' @param P Proportion of the level 4 units randomized to the treatment
#' @param g_4 Number of level 4 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units per level 4 unit
#' @param L Number of level 4 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param rho_3 proportion of variance in outcome that is between level 3 units
#' @param rho_4 proportion of variance in outcome that is between level 4 units
#' @param omega_3 Level 3 treatment effect heterogeneity: variance in treatment effect across Level 3 units, standardized by the Level 3 outcome variation
#' @param omega_4 Level 4 treatment effect heterogeneity: variance in treatment effect across Level 4 units, standardized by the Level 4 outcome variation
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BCRA4_2r <- function(R2_1, R2_2, R2_3T, R2_4T, P, g_4, n, J, K, L, rho_4, rho_3, rho_2, omega_3, omega_4, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), L - g_4 - 1), stats::qt(1 - alpha, L - g_4 - 1))
  T2 <- ifelse(power>=0.5, stats::qt(power, L-g_4-1), stats::qt(1-power,L-g_4-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_4*omega_4*(1-R2_4T)/L+rho_3*omega_3*(1-R2_3T)/(L*K)+rho_2*(1-R2_2)/(P*(1-P)*J*K*L)+(1-rho_2-rho_3-rho_4)*(1-R2_1)/(P*(1-P)*J*n*K*L))
  mdes
}

#' @title Model 4.4:  MDES Calculator for 4-Level Fixed Effects Blocked Cluster Random Assignment Designs (BCRA4_3f)—Treatment at Level 3
#' @description MDES for 4-Level Fixed Effects Blocked Cluster Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param R2_3 Proportion of variance in the Level 3 mean outcome explained by Block and Level 3 covariates
#' @param P Proportion of the level 3 units randomized to the treatment
#' @param g_3 Number of level 3 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units per level 4 unit
#' @param L Number of level 4 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param rho_3 proportion of variance in outcome that is between level 3 units
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BCRA4_3f <- function(R2_1, R2_2, R2_3, P, g_3, n, J, K, L, rho_2, rho_3, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), L*(K-2)-g_3), stats::qt(1 - alpha, L*(K-2)-g_3))
  T2 <- ifelse(power>=0.5, stats::qt(power, L*(K-2)-g_3), stats::qt(1-power, L*(K-2)-g_3))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_3*(1-R2_3)/(P*(1-P)*L*K)+rho_2*(1-R2_2)/(P*(1-P)*J*K*L)+(1-rho_2-rho_3)*(1-R2_1)/(P*(1-P)*J*n*K*L))
  mdes
}


#' @title Model 4.5:  MDES Calculator for 4-Level Random Effects Blocked Cluster Random Assignment Designs (BCRA4_3r)—Treatment at Level 3
#' @description MDES for 4-Level Random Effects Blocked Cluster Random Assignment Designs
#'
#' @param R2_1 Proportion of variance in level 1 outcome explained by level 1 covariates
#' @param R2_2 Proportion of variance in level 2 outcome explained by level 2 covariates
#' @param R2_3 Proportion of variance in the Level 3 outcome explained by the Level 3 covariates
#' @param R2_4T Proportion of variance in treatment effect between Level-4 blocks explained by Level-4 covariates
#' @param P Proportion of the level 3 units randomized to the treatment
#' @param g_4 Number of level 4 Covariates used
#' @param n Mean # of level 1 units per level 2 unit
#' @param J Mean # of level 2 units per level 3 unit
#' @param K Number of level 3 units per level 4 unit
#' @param L Number of level 4 units
#' @param rho_2 proportion of variance in outcome that is between level 2 units
#' @param rho_3 proportion of variance in outcome that is between level 3 units
#' @param rho_4 proportion of variance in outcome that is between level 4 units
#' @param omega_4 Level 4 treatment effect heterogeneity: variance in treatment effect across Level 4 units, standardized by the Level 4 outcome variation
#' @param alpha Probability of a type I error
#' @param power Statistical power (1 - probability of type II error)
#' @param tails One or two tailed test
#' @return Minimum Detectable Effect Size
#' @importFrom stats qt
BCRA4_3r <- function(R2_1, R2_2, R2_3, R2_4T, P, g_4, n, J, K, L, rho_4, rho_3, rho_2, omega_4, alpha = 0.05, power = 0.80, tails = 2){
  T1 <- ifelse(tails == 2, stats::qt(1 - (alpha/2), L - g_4 - 1), stats::qt(1 - alpha, L - g_4 - 1))
  T2 <- ifelse(power>=0.5, stats::qt(power, L-g_4-1), stats::qt(1-power,L-g_4-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)

  mdes <- M*sqrt(rho_4*omega_4*(1-R2_4T)/L+rho_3*(1-R2_3)/(P*(1-P)*L*K)+rho_2*(1-R2_2)/(P*(1-P)*J*K*L)+(1-rho_2-rho_3-rho_4)*(1-R2_1)/(P*(1-P)*J*n*K*L))
  mdes
}

