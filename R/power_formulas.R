
# 1.0 Individual Random assignment


IRA <- function(R2, P, k, n, alpha = 0.05, power = 0.80, tails = 1){
  T1 <- ifelse(tails == 2, qt(1 - (alpha/2), n - k - 2), qt(1 - alpha, n - k - 2))
  T2 <- ifelse(power >= 0.50, qt(power, n - k - 2), qt(1 - power, n - k - 2))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  mdes = M * sqrt((1 - R2) / (P * (1 - P) * n))
  mdes
}


# Model 2.1:  MDES Calculator for 2-Level Constant Effects Blocked Individual Random Assignment Designs (BIRA2_1c)—School Intercepts Only

BIRA2_1c <- function(R2_1, P, g, n, J, alpha = 0.05, power = 0.80, tails = 1){
  T1 <- ifelse(tails == 2, qt(1 - (alpha/2), n*J-J-g-1), qt(1 - alpha, n*J-J-g-1))
  T2 <- ifelse(power>=0.5, qt(power, n*J-J-g-1), qt(1-power,n*J-J-g-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  mdes = M * sqrt((1-R2_1)/(P * (1-P) * n * J))
  mdes
}



# Model 2.2:  MDES Calculator for 2-Level Fixed Effects Blocked Individual Random Assignment Designs (BIRA2_1f)— School Intercepts and Interactions with TREATMENT

BIRA2_1f <- function(R2_1, P, g, n, J, alpha = 0.05, power = 0.80, tails = 1){
  T1 <- ifelse(tails == 2, qt(1 - (alpha/2), n*J-2*J-g-1), qt(1 - alpha, n*J-2*J-g-1))
  T2 <- ifelse(power>=0.5, qt(power, n*J-2*J-g-1), qt(1-power,n*J-2*J-g-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  mdes = M * sqrt((1-R2_1)/(P * (1-P) * n * J))
  mdes
}

BIRA2_1f(R2_1 = 0.50, P = 0.50, g = 0, n = 55, J = 3, tails = 2)


# Model 2.3:  MDES Calculator for 2-Level Random Effects Blocked Individual Random Assignment (BIRA2_1r) Designs— Individuals Randomized within Blocks

BIRA2_1r <- function(R2_1, R2_2T, P, g, n, J, rho, omega, alpha = 0.05, power = 0.80, tails = 1){
  T1 <- ifelse(tails == 2, qt(1 - (alpha/2), J - g - 1), qt(1 - alpha, J - g - 1))
  T2 <- ifelse(power>=0.5, qt(power, J-g-1), qt(1-power,J-g-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  mdes = M * sqrt(rho * omega * (1-R2_2T)/ J + (1 - rho) * (1 - R2_1)/(P*(1-P)*J*n))
  mdes

}

BIRA2_1r(R2_1 = 0.0, R2_2T = 0, P = 0.50, g = 0, n = 80, J = 480, tails = 2, rho = 0.35, omega = 0.10)

#Model 2.4:  MDES Calculator for 3-Level Random Effects Blocked Individual Random Assignment Design (BIRA3_1r)—Individuals Randomized within Blocks

BIRA3_1r <- function(R2_1, R2_2T, R2_3T, P, g_3, n, J, K, rho_3, rho_2, omega_2, omega_3, alpha = 0.05, power = 0.80, tails = 1){
  T1 <- ifelse(tails == 2, qt(1 - (alpha/2), K - g_3 - 1), qt(1 - alpha, K - g_3 - 1))
  T2 <- ifelse(power>=0.5, qt(power, K-g_3-1), qt(1-power,K-g_3-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  mdes <- M*sqrt(rho_3*omega_3*(1-R2_3T)/K+rho_2*omega_2*(1-R2_2T)/(J*K)+(1-rho_2-rho_3)*(1-R2_1)/(P*(1-P)*J*n*K))
  mdes
}


BIRA3_1r(R2_1 = 0.0, R2_2T = 0, R2_3T = 0, P = 0.50, g_3 = 0, n = 80, J = 10, K = 100, rho_2 = 0.15, rho_3 = 0.20, omega_2 = 0.10, omega_3 = 0.10, tails = 2)

# Model 2.5:  MDES Calculator for 4-Level Random Effects Blocked Individual Random Assignment Design (BIRA4_1r)—Individuals Randomized within Blocks

BIRA4_1r <- function(R2_1, R2_2T, R2_3T, R2_4T, P, g_4, n, J, K, L, rho_4, rho_3, rho_2, omega_2, omega_3, omega_4, alpha = 0.05, power = 0.80, tails = 1){
  T1 <- ifelse(tails == 2, qt(1 - (alpha/2), L - g_4 - 1), qt(1 - alpha, L - g_4 - 1))
  T2 <- ifelse(power>=0.5, qt(power, L-g_4-1), qt(1-power,L-g_4-1))
  M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
  print(T1)
  print(T2)
  print(M)
  mdes <- M*sqrt(rho_4*omega_4*(1-R2_4T)/L+rho_3*omega_3*(1-R2_3T)/(L*K)+rho_2*omega_2*(1-R2_2T)/(J*K*L)+(1-rho_2-rho_3-rho_4)*(1-R2_1)/(P*(1-P)*J*n*K*L))
  mdes
}


BIRA4_1r(R2_1 = 0.50, R2_2T = 0.50, R2_3T = 0.50, R2_4T = 0.50, P = 0.50, g_4 = 1, n = 10, J = 4, K = 4, L = 20, rho_2 = 0.15, rho_3 = 0.15, rho_4 = 0.05, omega_2 = 0.50, omega_3 = 0.50, omega_4 = 0.50 ,tails = 2)


# Model 3.1:  MDES Calculator for Two-Level Cluster Random Assignment Design (CRA2_2)— Treatment at Level 2

CRA2_2r <- function(R2_1, R2_2T, P, g, n, J, rho, alpha = 0.05, power = 0.80, tails = 1){
    T1 <- ifelse(tails == 2, qt(1 - (alpha/2), J - g - 2), qt(1 - alpha, J - g - 2))
    T2 <- ifelse(power>=0.5, qt(power, J-g-2), qt(1-power,J-g-2))
    M <- ifelse(power >= 0.50, T1 + T2, T1 - T2)
    mdes <- M*sqrt(rho*(1-R2_2T)/(P*(1-P)*J)+(1-rho)*(1-R2_1)/(P*(1-P)*J*n))
    mdes

}

CRA2_2r(R2_1 = 0.01, R2_2T = 0.13, P = 0.50, g = 4, n = 90, J = 9, rho = 0.02, tails = 2)




