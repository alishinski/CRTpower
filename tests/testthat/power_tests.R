# Individual Random Assignment
expect_equal(round(IRA(R2 = 0, P = 0.50,  n = 400, k = 0), digits = 3), 0.281)

# 2-Level Constant Effects Blocked Individual Random Assignment Designs (BIRA2_1c)—School Intercepts Only
expect_equal(round(BIRA2_1c(R2_1 = 0.00, P = 0.50, g = 0, n = 55, J = 14), digits = 3), 0.202)

# 2-Level Fixed Effects Blocked Individual Random Assignment Designs (BIRA2_1f)— School Intercepts and Interactions with TREATMENT
expect_equal(round(BIRA2_1f(R2_1 = 0.50, P = 0.50, g = 0, n = 55, J = 3), digits = 3), 0.310)

# 2-Level Random Effects Blocked Individual Random Assignment (BIRA2_1r) Designs— Individuals Randomized within Blocks
expect_equal(round(BIRA2_1r(R2_1 = 0.0, R2_2T = 0, P = 0.50, g = 0, n = 80, J = 480, tails = 2, rho = 0.35, omega = 0.10), digits = 3), 0.033)

# 3-Level Random Effects Blocked Individual Random Assignment Design (BIRA3_1r)—Individuals Randomized within Blocks
expect_equal(round(BIRA3_1r(R2_1 = 0.0, R2_2T = 0, R2_3T = 0, P = 0.50, g_3 = 0, n = 80, J = 10, K = 100, rho_2 = 0.15, rho_3 = 0.20, omega_2 = 0.10, omega_3 = 0.10, tails = 2), digits = 3), 0.045)

# 4-Level Random Effects Blocked Individual Random Assignment Design (BIRA4_1r)—Individuals Randomized within Blocks
expect_equal(round(BIRA4_1r(R2_1 = 0.50, R2_2T = 0.50, R2_3T = 0.50, R2_4T = 0.50, P = 0.50, g_4 = 1, n = 10, J = 4, K = 4, L = 20, rho_2 = 0.15, rho_3 = 0.15, rho_4 = 0.05, omega_2 = 0.50, omega_3 = 0.50, omega_4 = 0.50 ,tails = 2), digits = 3), 0.119)

# Two-Level Cluster Random Assignment Design (CRA2_2)— Treatment at Level 2
expect_equal(round(CRA2_2r(R2_1 = 0.01, R2_2T = 0.13, P = 0.50, g = 4, n = 90, J = 9, rho = 0.02, tails = 2), digits = 3), 0.466)

# 3 Level Cluster  Random Assignment Designs (CRA3_3r) Treatment at Level 3
expect_equal(round(CRA3_3r(R2_1 = 0.37, R2_2 = 0.53, R2_3 = 0.87, P = 0.50, g_3 = 1, n = 20, J = 2, K = 66, rho_3 = 0.38, rho_2 = 0.10, tails = 2), digits = 3), 0.199)

# 4-Level Cluster  Random Assignment Designs (CRA4_4r)—Treatment at Level 4
expect_equal(round(CRA4_4r(R2_1 = 0.50, R2_2 = 0.50, R2_3 = 0.50, R2_4 = 0.50, P = 0.50, g_4 = 1, n = 10, J = 2, K = 3, L = 20, rho_4 = 0.05, rho_3 = 0.05, rho_2 = 0.10, tails = 2) , digits = 3), 0.292)

# 3-Level Random Effects Blocked Cluster Random Assignment Designs (BCRA3_2r)—Treatment at Level 2
expect_equal(round(BCRA3_2f(R2_1 = 0.50, R2_2 = 0.50, P = 0.50, g_2 = 1, n = 20, J = 44, K = 5, rho_2 = 0.10, tails = 2), digits = 3), 0.102)

# 3-Level Random Effects Blocked Cluster Random Assignment Designs (BCRA3_2r)—Treatment at Level 2
expect_equal(round(BCRA3_2r(R2_1 = 0.37, R2_2 = 0.53, R2_T3 = 0.00, P = 0.50, g_3 = 0, n = 20, J = 2, K = 64, rho_2 = 0.10, rho_3 = 0.38, omega_3 = 0.50, tails = 2), digits = 3), 0.200)

# 4-Level Random Effects Block Random Assignment Designs (BCRA4_2r)—Treatment at Level 2
expect_equal(round(BIRA4_2r(R2_1 = 0.50, R2_2 = 0.50, R2_3T = 0.50, R2_4T = 0.50, P = 0.50, g_4 = 1, n = 10, J = 4, K = 4, L = 20, rho_4 = 0.05, rho_3 = 0.15, rho_2 = 0.15, omega_4 = 0.50, omega_3 = 0.50, tails = 2), digits = 3), 0.146)

# 4-Level Fixed Effects Blocked Cluster Random Assignment Designs (BCRA4_3f)—Treatment at Level 3
expect_equal(round(BCRA4_2r(R2_1 = 0.50, R2_2 = 0.50, R2_3 = 0.50, P = 0.50, g_3 = 2, n = 10, J = 4, K = 4, L = 15, rho_2 = 0.15, rho_3 = 0.15, tails = 2) , digits = 3), 0.240)

# 4-Level Random Effects Blocked Cluster Random Assignment Designs (BCRA4_3r)—Treatment at Level 3
expect_equal(round(BCRA4_3r(R2_1 = 0.50, R2_2 = 0.50, R2_3 = 0.50, R2_4T = 0.50, P = 0.50, g_4 = 3, n = 10, J = 4, K = 20, L = 20, rho_4 = 0.05, rho_3 = 0.15, rho_2 = 0.15, omega_4 = 0.50, tails = 2), digits = 3), 0.121)
