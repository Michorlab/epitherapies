

#################

graphics.off()  # Close all open graphics devices
rm(list = ls())     # clear all
cat("\014")     # clc


##################  


# Load necessary libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(doParallel)
library(tidyverse)
library(readxl)
library(viridis)
library(DiceKriging)
library(DiceOptim)
library(GA)  # Genetic Algorithm package
library(future.apply)
plan(multisession)   # portable parallelization


set.seed(1234) # Ensure reproducibility

# Output location

output_location <- "/Users/..." # to be modified
output_location_figures <- "/Users/..." # to be modified
if(!file.exists(output_location)) dir.create(output_location)
if(!file.exists(output_location_figures)) dir.create(output_location_figures)


# Load ODE PD model parameters from saved Stan fitting
posterior_samples <- readRDS(paste0(output_location, "posterior_samples.rds"))




patient_to_plot <- 1 #Select 1,2,3,4


# Load patient-specific posterior samples - Select patietn 1,2,3, or 4

theta_patient <- readRDS(
  file.path(output_location, "posterior_samples_patient_4.Rds")
)


# ------------------------------------------------------------
# Posterior mean of the four patient-specific parameters
# (old k1, k11, k14, k17)
# ------------------------------------------------------------

theta_patient_mean <- colMeans(theta_patient)

# ------------------------------------------------------------
# Population means for the remaining parameters
# ------------------------------------------------------------

theta_population_mean <- colMeans(posterior_samples)

# Final individualized parameter vector
k_patient <- theta_population_mean

k_patient[c(1,11,14,17)] <- theta_patient_mean





# Draw parameter samples and run the ODE for each sample
n_patients <- 1


# Initialize a data frame to store parameters for each in-silico patient
params_samples <- data.frame(matrix(ncol = 30, nrow = n_patients))
colnames(params_samples) <- paste0("k", 1:30)



# Loop to draw random parameter samples for each in-silico patient
for (i in 1:n_patients) {
  # Sample parameter values from the posterior
  
  # Store the sampled parameters in the data frame
  params_samples[i, ] <- k_patient
}


################################################################################
# ============================================================
# Generate patient-specific PK parameters for EZH2i
# ============================================================

# Population values
Vd_pop   <- 681     # L (FDA Vss/F)

CL0_pop  <- 126      # L/h (Day 1)
CLss_pop <- 274      # L/h (Day 15 steady state)

ka_pop   <- 2.42     # h^-1

# Autoinduction rate constant
# ~95% induction complete by Day 15
k_ind_pop <- 3/(15*24)

# ============================================================
# Interpatient variability
# ============================================================

# Uncomment to disable PK variability --> variability no needed for patient specific simulations
CV_Vd <- 0
CV_CL <- 0

# Log-normal standard deviations
omega_Vd <- sqrt(log(1 + CV_Vd^2))
omega_CL <- sqrt(log(1 + CV_CL^2))

# ============================================================
# Storage
# ============================================================

pk_params_samples <- data.frame(
  Vd    = numeric(n_patients),
  CL0   = numeric(n_patients),
  CLss  = numeric(n_patients),
  ka    = numeric(n_patients),
  k_ind = numeric(n_patients)
)

# ============================================================
# Generate one PK parameter set per virtual patient
# ============================================================

for(i in 1:n_patients){
  
  # Volume
  Vd_i <- Vd_pop * exp(
    rnorm(1, mean = 0, sd = omega_Vd)
  )
  
  # Day 1 clearance
  CL0_i <- CL0_pop * exp(
    rnorm(1, mean = 0, sd = omega_CL)
  )
  
  # Day 15 clearance
  CLss_i <- CLss_pop * exp(
    rnorm(1, mean = 0, sd = omega_CL)
  )
  
  # Ensure autoinduction always increases clearance
  CLss_i <- max(CLss_i, CL0_i)
  
  # Fixed values
  ka_i    <- ka_pop
  k_ind_i <- k_ind_pop
  
  # Store
  pk_params_samples[i, ] <- c(
    Vd_i,
    CL0_i,
    CLss_i,
    ka_i,
    k_ind_i
  )

  
}


# ============================================================
# Derived quantities
# ============================================================

pk_params_samples$ke0 <-
  pk_params_samples$CL0 /
  pk_params_samples$Vd

pk_params_samples$ke_ss <-
  pk_params_samples$CLss /
  pk_params_samples$Vd

pk_params_samples$Tmax_ss <-
  log(pk_params_samples$ka /
        pk_params_samples$ke_ss) /
  (pk_params_samples$ka -
     pk_params_samples$ke_ss)

pk_params_samples$t_half_ss <-
  log(2) / pk_params_samples$ke_ss


################################################################################
# ============================================================
# Generate patient-specific PK parameters for AKTi
# ============================================================

# Population values
CL_pop <- 162      # L/h
ka_pop <- 1.84     # 1/h
t0_pop <- 0.348    # h

# Between-subject variances (reported in the paper)
# omega2_CL <- 0.0729
# omega2_ka <- 1.74
# omega2_t0 <- 2.23
# omega2_FI <- 0.139


# Uncomment to disable PK variability --> variability no needed for patient specific simulations
omega2_CL <- 0
omega2_ka <- 0
omega2_t0 <- 0
omega2_FI <- 0

# Corresponding standard deviations
omega_CL <- sqrt(omega2_CL)
omega_ka <- sqrt(omega2_ka)
omega_t0 <- sqrt(omega2_t0)
omega_FI <- sqrt(omega2_FI)

# ============================================================
# Storage
# ============================================================

pk_AKTi_params_samples <- data.frame(
  CL    = numeric(n_patients),
  ka    = numeric(n_patients),
  t0    = numeric(n_patients),
  FI_eta = numeric(n_patients)
)

# ============================================================
# Generate one PK parameter set per virtual patient
# ============================================================

for(i in 1:n_patients){
  
  # Sample patient-specific PK parameters
  CL_i <- CL_pop * exp(rnorm(1, mean = 0, sd = omega_CL))
  ka_i <- ka_pop * exp(rnorm(1, mean = 0, sd = omega_ka))
  t0_i <- pmin(4, pmax(0.05, t0_pop * exp(rnorm(1, mean = 0, sd = omega_t0))))
  
  # Patient-specific relative bioavailability term
  FI_eta_i <- exp(rnorm(1, mean = 0, sd = omega_FI))
  
  # Store
  pk_AKTi_params_samples[i, ] <- c(
    CL_i,
    ka_i,
    t0_i,
    FI_eta_i
  )
  
}


###############################################################################

# Define the pharmacokinetic model function for EZH2i

define_PK_EZHi_fun <- function(D0_e, dose_times, total_time, times_between_doses, pk_params) {


  pk_EZHi_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      
      
      # Time-dependent clearance
      CL_t <- CL0 + (CLss - CL0) * (1 - exp(-k_ind * time))
      
      # Time-dependent elimination rate
      ke_t <- CL_t / Vd
      
      
      # Differential equations for the one-compartment model with first-order absorption
      dA_GI <- -ka * A_GI  # Amount in gut decreases due to absorption
      dC_plasma <- (ka * A_GI) / Vd - ke_t * C_plasma  # Change in concentration in plasma

      list(c(dA_GI, dC_plasma))
    })
  }

  # Parameters
  parameters <- c(
    ka    = as.numeric(pk_params$ka), # Absorption rate constant (1/hour)
    CL0   = as.numeric(pk_params$CL0)* 0.33,
    CLss  = as.numeric(pk_params$CLss)* 0.33,
    k_ind = as.numeric(pk_params$k_ind),
    Vd    = as.numeric(pk_params$Vd) * 0.33  # Volume of distribution (L)
  )
  
  # Initial conditions
  initial_conditions <- c(
    A_GI = D0_e[1]*0.33,   # Initial amount of drug in GI tract* Bioavailability
    C_plasma = 0  # Initial concentration in plasma 
  )

  # Initialize output storage
  output_df_EZHi <- data.frame()

  # Current state of the system
  current_state <- initial_conditions

  # Current time to start the loop
  current_time <- 0


  #index for dose
  iijjii <- 2

  # Loop over each dosing time
  for (dose_time in dose_times) {
    # Define time sequence for current loop
    time_segment <- seq(current_time, dose_time, by = times_between_doses)

    # Simulate the model for the current segment
    output_segment <- ode(
      func = pk_EZHi_model,
      times = time_segment,
      y = current_state,
      parms = parameters,
      method = "bdf",
      maxsteps = 5000  # Increased maximum number of steps
    )

    # Convert output to a data frame and store it
    output_df_EZHi <- rbind(output_df_EZHi, as.data.frame(output_segment))

    # Update current state to the last row of the output segment
    current_state_first <- tail(output_segment, 1)[-1]  # Remove time column

    # Add dose to the gut compartment
    current_state_first[1] <- current_state_first[1] + D0_e[iijjii]*0.33 #bioavailability

    current_state <- c(
      A_GI = current_state_first[1],   # Initial amount of drug in GI tract
      C_plasma = current_state_first[2]  # Initial concentration in plasma
    )

    # Update current time
    current_time <- dose_time


    #Update index for dose
    iijjii <- iijjii + 1
  }

  # Simulate the final segment after the last dose until the end of the total time
  time_segment <- seq(current_time, total_time, by = times_between_doses)
  output_segment <- ode(
    func = pk_EZHi_model,
    times = time_segment,
    y = current_state,
    parms = parameters,
    method = "bdf",
    maxsteps = 5000  # Increased maximum number of steps
  )


  output_df_EZHi <- rbind(output_df_EZHi, as.data.frame(output_segment))

  
  # Concentrations
  
  Me <- 572.75
  output_df_EZHi <- output_df_EZHi %>%
    mutate(C_plasma = (C_plasma * 10000) / Me)


  # Return the data frame
  return(output_df_EZHi)

}



################################################################################

# Define the pharmacokinetic model function for AKTi
define_PK_AKTi_fun <- function(D0_a, dose_times, total_time, times_between_doses, pk_params) {


  # ============================================================
  # PARAMETERS
  # ============================================================

  parameters <- list(
    CLI = as.numeric(pk_params$CL),    # Clearance from Central compartment (L/h)
    V2I = 1230,   # Volume of distribution of Central compartment (L)
    V3I = 2590,   # Volume of distribution Peripheral 1 (L)
    V4I = 4340,   # Volume of distribution Peripheral 2 (L)
    Q3I = 76.6,   # Intercompartmental clearance (Central <-> P1)
    Q4I = 3.26,   # Intercompartmental clearance (Central <-> P2)
    kaI = as.numeric(pk_params$ka),   # First-order absorption rate constant (1/h)

    t0 = as.numeric(pk_params$t0),          # Zero-order absorption duration (h)
    FI_eta = as.numeric(pk_params$FI_eta),
    theta_FI_MD = 0.201,
    
    dose_times = dose_times,
    D0_a = D0_a
  )
  

  # ============================================================
  # INITIAL CONDITIONS
  # ============================================================

  state <- c(
    A_depot_2 = 0,
    A_central = 0,
    A_peripheral1 = 0,
    A_peripheral2 = 0
  )

  # ============================================================
  # PK MODEL
  # ============================================================

  pk_model <- function(time, state, parameters) {
    with(
      c(as.list(state), parameters),
      {

      # -----------------------------
      # Zero-order input rate R_zero
      # -----------------------------
      R_zero <- 0
      
      MD <- ifelse(time < 48, 0, 1)
      FI_t <- (1 + theta_FI_MD * MD) * FI_eta

      for (i in seq_along(dose_times)) {
        t_i <- dose_times[i]
        D_i <- D0_a[i]

        if (time >= t_i && time < (t_i + t0)) {
          
          R_zero <- R_zero + (FI_t * D_i / t0)
        }
      }

      # -----------------------------
      # ODEs
      # -----------------------------
      
      dA_depot_2   <- R_zero - kaI * A_depot_2

      dA_central <- kaI * A_depot_2 -
        (CLI / V2I) * A_central -
        (Q3I / V2I) * A_central + (Q3I / V3I) * A_peripheral1 -
        (Q4I / V2I) * A_central + (Q4I / V4I) * A_peripheral2

      dA_peripheral1 <- (Q3I / V2I) * A_central - (Q3I / V3I) * A_peripheral1
      dA_peripheral2 <- (Q4I / V2I) * A_central - (Q4I / V4I) * A_peripheral2

      list(c(dA_depot_2, dA_central, dA_peripheral1, dA_peripheral2))
    })
  }

  # ============================================================
  # SOLVE THE SYSTEM
  # ============================================================

  times <- seq(0, total_time, by = times_between_doses)

  output_df_AKTi <- as.data.frame(ode(
    y = state,
    times = times,
    func = pk_model,
    parms = parameters,
    method = "bdf",
    maxsteps = 5000
  ))

  # ============================================================
  # CONCENTRATIONS
  # ============================================================

  output_df_AKTi$C_central <- output_df_AKTi$A_central / parameters$V2I
  
  output_df_AKTi$C_total <- (output_df_AKTi$A_central + output_df_AKTi$A_peripheral1 + output_df_AKTi$A_peripheral2) /
    (parameters$V2I + parameters$V3I + parameters$V4I)
  
  MW <- 458
  output_df_AKTi <- output_df_AKTi %>% mutate(C_central = (C_central * 10000) / MW)
  output_df_AKTi <- output_df_AKTi %>% mutate(C_total = (C_total * 10000) / MW)


  # ============================================================
  # RETURN
  # ============================================================

  return(output_df_AKTi)
}


################################################################################
################################################################################

################################################################################
# AKTi schedule:
# - no dosing during first 5 pretreatment days
# - then 3 weeks
#   days_on = c(days in week1, days in week2, days in week3)
#   dose_a  = c(dose_week1, dose_week2, dose_week3)
################################################################################

generate_dosing_schedule_D0a <- function(dose_a, days_on,
                                         pretreat_days = 5, treat_days = 21) {
  
  total_days <- pretreat_days + treat_days
  schedule <- rep(0, total_days)
  
  # ---- Week 1 ----
  if (days_on[1] > 0) {
    schedule[(pretreat_days + 1):(pretreat_days + days_on[1])] <- dose_a[1]
  }
  
  # ---- Week 2 ----
  if (days_on[2] > 0) {
    week2_start <- pretreat_days + 8
    schedule[week2_start:(week2_start + days_on[2] - 1)] <- dose_a[2]
  }
  
  # ---- Week 3 ----
  if (days_on[3] > 0) {
    week3_start <- pretreat_days + 15
    schedule[week3_start:(week3_start + days_on[3] - 1)] <- dose_a[3]
  }
  
  return(schedule)
}


################################################################################
# EZH2i schedule:
# dose_e = c(dose_pretreat, dose_week1, dose_week2, dose_week3)
################################################################################

generate_dosing_schedule_D0e <- function(dose_e,
                                         pretreat_days = 5, treat_days = 21) {
  
  total_days <- pretreat_days + treat_days
  schedule <- numeric(total_days)
  
  # Pretreatment: days 1–5
  schedule[1:pretreat_days] <- dose_e[1]
  
  # Week 1: days 6–12
  schedule[(pretreat_days + 1):(pretreat_days + 7)] <- dose_e[2]
  
  # Week 2: days 13–19
  schedule[(pretreat_days + 7 + 1):(pretreat_days + 14)] <- dose_e[3]
  
  # Week 3: days 20–26
  schedule[(pretreat_days + 14 + 1):(pretreat_days + 21)] <- dose_e[4]
  
  return(schedule)
}



################################################################################
# BID administration
################################################################################

expand_daily_doses <- function(d) {
  rep(d, each = 2)   # turns 26 days → 52 BID doses
}




################################################################################

# Global storage to keep track of y1 values
results_storage <- list(y1_values = numeric(), cost_values = numeric(), D0_e_pt = numeric(), D0_e_w1 = numeric(), D0_e_w2 = numeric(), D0_e_w3 = numeric(), D0_a_w1 = numeric(), D0_a_w2 = numeric(), D0_a_w3 = numeric(), AKTi_Days_Week1 = numeric(), AKTi_Days_Week2 = numeric(), AKTi_Days_Week3 = numeric())



evaluation_counter <- 0  # Global counter

# GA Objective Function
ga_objective <- function(x, limit_toxicity_e, limit_toxicity_a) {
  
  D0_e_day <- generate_dosing_schedule_D0e(c(x[1], x[2], x[3], x[4]))
  D0_a_day <- generate_dosing_schedule_D0a(c(x[5], x[6], x[7]), c(x[8], x[9], x[10]))
  
  # Expand to two doses/day
  D0_e <- expand_daily_doses(D0_e_day)
  D0_a <- expand_daily_doses(D0_a_day)
  
  
  dose_times <- sort(c(
    1,      1 + 12,
    1*24,   1*24 + 12,
    2*24,   2*24 + 12,
    3*24,   3*24 + 12,
    4*24,   4*24 + 12,
    5*24,   5*24 + 12,
    6*24,   6*24 + 12,
    7*24,   7*24 + 12,
    8*24,   8*24 + 12,
    9*24,   9*24 + 12,
    10*24,  10*24 + 12,
    11*24,  11*24 + 12,
    12*24,  12*24 + 12,
    13*24,  13*24 + 12,
    14*24,  14*24 + 12,
    15*24,  15*24 + 12,
    16*24,  16*24 + 12,
    17*24,  17*24 + 12,
    18*24,  18*24 + 12,
    19*24,  19*24 + 12,
    20*24,  20*24 + 12,
    21*24,  21*24 + 12,
    22*24,  22*24 + 12,
    23*24,  23*24 + 12,
    24*24,  24*24 + 12,
    25*24,  25*24 + 12
  ))

  total_time <- 26 * 24

results_list <- future_lapply(1:n_patients, function(i) {

  params <- as.list(params_samples[i, ])
  
  pki <- pk_AKTi_params_samples[i, ]
  
  dt_i <- max(
    0.05,
    min(as.numeric(pki$t0), 0.1)
  )
  # dt_i <- 0.1

  # Solve the PK models for AKTi and EZHi
  
  output_df_AKTi <- define_PK_AKTi_fun(D0_a, dose_times, total_time, dt_i, pk_params = pk_AKTi_params_samples[i, ])

  output_df_EZHi <- define_PK_EZHi_fun(D0_e, dose_times[-1], total_time, dt_i, pk_params = pk_params_samples[i, ])
  # 
  
  
  
  

  # Ensure unique time points
  pk_AKTi_unique <- output_df_AKTi[!duplicated(output_df_AKTi[, "time"]), ]
  pk_EZHi_unique <- output_df_EZHi[!duplicated(output_df_EZHi[, "time"]), ]

  # Interpolating functions for concentrations
  AKTi_interp <- approxfun(pk_AKTi_unique[, "time"], pk_AKTi_unique[, "C_central"], rule = 2)
  EZHi_interp <- approxfun(pk_EZHi_unique[, "time"], pk_EZHi_unique[, "C_plasma"], rule = 2)
  
  
  
  # Define the time sequence for the PD model simulation
  times <- seq(0, 624, by = dt_i)
  
  PD_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      k31 <- AKTi_interp(time)
      k32 <- EZHi_interp(time)
      dy1_dt <- (k1 * (y6 / (k2 + y6)) * exp(-k16 * time) * y8) - k15 * y1 * (1 - (y1 / 100))
      dy2_dt <- k3 * (k4 / (k4 + y4)) * y3 - k5 * y2
      dy3_dt <- (k6 +k25*y3) * (50 - y3) - (k29+k7 * (y5 / (k8 + y5)) *(1+k26*(50-y3)))* y3
      dy4_dt <- k9 - k10 * (k31 / (k11 + k31)) * y4 - k5 * y4
      dy5_dt <- k12 - k13 * (k32 / (k14 + k32)) * y5 - k5 * y5
      dy6_dt <- k19 * (k20 / (k20 + y4)) * y7 - k21 * y6
      dy7_dt <- (k22 + k27 * y7) * (50 - y7) - (k30 + k23 * (y5 / (k24 + y5)) * (1 + k28 * (50 - y7))) * y7
      dy8_dt <- k17 * (y2 / (k18 + y2)) * (100 - y1 - y8) - (k1 * (y6 / (k2 + y6)) * exp(-k16 * time) - k15 * (y1 / 100)) * y8
      list(c(dy1_dt, dy2_dt, dy3_dt, dy4_dt, dy5_dt, dy6_dt, dy7_dt, dy8_dt))
    })
  }
  
  
  
  y_initial <- c(y1 = 10.6765, y2 = 0, y3 = 0, y4 = 1, y5 = 1, y6 = 0, y7 = 0, y8 = 0)

  pd_output <- ode(
    y     = y_initial,
    times = times,
    func  = PD_model,
    parms = params
  )

  data.frame(
    patient_id = i,
    time       = tail(pd_output[, "time"], 1),
    y1         = tail(pd_output[, "y1"], 1)
  )
})

  # Combine all rows
  results_df <- do.call(rbind, results_list)
  
  # Compute original y1 metric
  y1_value <- mean(results_df$y1)
  
  penalty <- 10000
  
  # Compute total cost function
  
  # Compute cumulative drug dosage (use daily schedules, not BID-expanded)
  total_EZ <- sum(D0_e)
  total_AK <- sum(D0_a)
  
  # Base objective (maximize % dying cells = minimize surviving y1)
  cost_base <- -(y1_value / 100)
  
  # Apply toxicity penalty if either drug exceeds its limit
  if (total_EZ > limit_toxicity_e || total_AK > limit_toxicity_a) {
    cost_function <- cost_base + penalty
  } else {
    cost_function <- cost_base
  }
  
  
  
  # Store y1 for later analysis
  results_storage$y1_values <<- c(results_storage$y1_values, y1_value)
  results_storage$cost_values <<- c(results_storage$cost_values, cost_function)
  results_storage$D0_e_pt <<- c(results_storage$D0_e_pt, x[1])
  results_storage$D0_e_w1 <<- c(results_storage$D0_e_w1, x[2])
  results_storage$D0_e_w2 <<- c(results_storage$D0_e_w2, x[3])
  results_storage$D0_e_w3 <<- c(results_storage$D0_e_w3, x[4])
  results_storage$D0_a_w1 <<- c(results_storage$D0_a_w1, x[5])
  results_storage$D0_a_w2 <<- c(results_storage$D0_a_w2, x[6])
  results_storage$D0_a_w3 <<- c(results_storage$D0_a_w3, x[7])
  results_storage$AKTi_Days_Week1 <<- c(results_storage$AKTi_Days_Week1, x[8])
  results_storage$AKTi_Days_Week2 <<- c(results_storage$AKTi_Days_Week2, x[9])
  results_storage$AKTi_Days_Week3 <<- c(results_storage$AKTi_Days_Week3, x[10])
  
  
  return(cost_function)
  
}



################### GA Optimization


# Set optimization bounds
lower_bounds <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
upper_bounds <- c(800, 800, 800, 800, 400, 400, 400, 7, 7, 7)


# Set weights
limit_toxicity_e_L3 <- 41600
limit_toxicity_a_L3 <- 9600

limit_toxicity_ee <- limit_toxicity_e_L3
limit_toxicity_aa <- limit_toxicity_a_L3


# Run GA Optimization
ga_result <- ga(
  type = "real-valued",  # Specifies that the variables are continuous (not binary or integer)

  # fitness = function(x) -ga_objective(x, lambda_days = lambdad, lambda_concentration = lambdaconc),
  fitness = function(x) {
    # Ensure last three variables are integers and within bounds
    x[8] <- pmin(pmax(round(x[8]), lower_bounds[8]), upper_bounds[8])
    x[9] <- pmin(pmax(round(x[9]), lower_bounds[9]), upper_bounds[9])
    x[10] <- pmin(pmax(round(x[10]), lower_bounds[10]), upper_bounds[10])

    -ga_objective(x, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)
  },
  # Fitness function: The objective function to minimize.
  # Since GA maximizes by default, we negate it (-ga_objective) to perform minimization.

  lower = lower_bounds,  # Lower bounds for each parameter (minimum values they can take)
  upper = upper_bounds,  # Upper bounds for each parameter (maximum values they can take)

  popSize = 70,  
  # Population size: Number of candidate solutions (individuals) per generation.
  
  maxiter = 70,  #100
  # Maximum number of generations (iterations of evolution).
  
  pmutation = 0.2,  
  # Probability of mutation: Controls randomness in solution evolution.
  
  pcrossover = 0.8,  
  # Probability of crossover: Controls how solutions are recombined.
  
  run = 30,  
  # Stopping criteria: GA stops if there is no improvement for this many consecutive generations.

  monitor = TRUE  # This enables the default monitoring output
)




# Extract best solution
best_solution <- ga_result@solution
best_solution[8:10] <- round(best_solution[8:10])  # Round last three variables


best_cost <- ga_result@fitnessValue

# Store results
tested_treatments <- data.frame(
  Dose_EZHi_pt = results_storage$D0_e_pt,
  Dose_EZHi_w1 = results_storage$D0_e_w1,
  Dose_EZHi_w2 = results_storage$D0_e_w2,
  Dose_EZHi_w3 = results_storage$D0_e_w3,
  Dose_AKTi_w1 = results_storage$D0_a_w1,
  Dose_AKTi_w2 = results_storage$D0_a_w2,
  Dose_AKTi_w3 = results_storage$D0_a_w3,
  AKTi_Days_Week1 = results_storage$AKTi_Days_Week1,
  AKTi_Days_Week2 = results_storage$AKTi_Days_Week2,
  AKTi_Days_Week3 = results_storage$AKTi_Days_Week3,
  Y1_T624 = results_storage$y1_values,  # Store y1 values separately
  Cost_Function = results_storage$cost_values
)

best_treatment <- data.frame(
  Dose_EZHi_pt = best_solution[1],
  Dose_EZHi_w1 = best_solution[2],
  Dose_EZHi_w2 = best_solution[3],
  Dose_EZHi_w3 = best_solution[4],
  Dose_AKTi_w1 = best_solution[5],
  Dose_AKTi_w2 = best_solution[6],
  Dose_AKTi_w3 = best_solution[7],
  AKTi_Days_Week1 = best_solution[8],
  AKTi_Days_Week2 = best_solution[9],
  AKTi_Days_Week3 = best_solution[10],
  Cost_Function = best_cost,
  Y1_T624 = results_storage$y1_values[which.min(results_storage$cost_values)]
)



###############################################
# MANUAL TESTS: EZH2-ONLY AND AKTi-ONLY CASES
###############################################

cat("Running manual EZH2-only and AKTi-only tests...\n")

n_before <- length(results_storage$y1_values)

# 5 doses for each drug
ezh2_doses <- seq(lower_bounds[1], upper_bounds[1], length.out = 3)
akti_doses <- seq(lower_bounds[5], upper_bounds[5], length.out = 3)

# =====================================================
# C) STANDARD - LEVEL 3
# =====================================================

dos_e <- 800
dos_a <- 400

x_manual_std <- c(
  dos_e, dos_e, dos_e, dos_e,     # EZH2 dose
  dos_a, dos_a, dos_a,   # AKTi dose
  4,4,4   # AKTi schedule (4-4-4)
)
-ga_objective(x_manual_std, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)

# =====================================================
# Append manual results into tested_treatments
# =====================================================

n_after <- length(results_storage$y1_values)
manual_idx <- (n_before + 1):n_after

manual_treatments <- data.frame(
  Dose_EZHi_pt = results_storage$D0_e_pt[manual_idx],
  Dose_EZHi_w1 = results_storage$D0_e_w1[manual_idx],
  Dose_EZHi_w2 = results_storage$D0_e_w2[manual_idx],
  Dose_EZHi_w3 = results_storage$D0_e_w3[manual_idx],
  Dose_AKTi_w1 = results_storage$D0_a_w1[manual_idx],
  Dose_AKTi_w2 = results_storage$D0_a_w2[manual_idx],
  Dose_AKTi_w3 = results_storage$D0_a_w3[manual_idx],
  AKTi_Days_Week1 = results_storage$AKTi_Days_Week1[manual_idx],
  AKTi_Days_Week2 = results_storage$AKTi_Days_Week2[manual_idx],
  AKTi_Days_Week3 = results_storage$AKTi_Days_Week3[manual_idx],
  Y1_T624   = results_storage$y1_values[manual_idx],
  Cost_Function = results_storage$cost_values[manual_idx]
)

tested_treatments <- rbind(tested_treatments, manual_treatments)

####################################################################
# Removed repeated rows

tested_treatments <- tested_treatments %>%
  distinct(
    Dose_EZHi_pt,
    Dose_EZHi_w1,
    Dose_EZHi_w2,
    Dose_EZHi_w3,
    Dose_AKTi_w1,
    Dose_AKTi_w2,
    Dose_AKTi_w3,
    AKTi_Days_Week1,
    AKTi_Days_Week2,
    AKTi_Days_Week3,
    .keep_all = TRUE
  )

#####################################################################

# SAVE FINALs

# Define file paths
tested_treatments_file <- "/Users/....rds"
best_treatment_file <- "/Users/....rds"


## Save the objects
# saveRDS(tested_treatments, file = tested_treatments_file)
# saveRDS(best_treatment, file = best_treatment_file)

## To reload them later
tested_treatments <- readRDS(tested_treatments_file)
best_treatment <- readRDS(best_treatment_file)
#
#
set.seed(1234) # Ensure reproducibility


# ###################### SUMMARY PLOTs



tested_treatments <- tested_treatments %>%
  mutate(
    EZHi_total = round(Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3, 0),
    AKTi_total = round(Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3, 0),
    AKTi_days_total = AKTi_Days_Week1 + AKTi_Days_Week2 + AKTi_Days_Week3
  )
#
############################################################
# 2) Remove NEGATIVE cost-function solutions
############################################################


tested_treatments <- tested_treatments %>%
  filter(Cost_Function <= 0)

# # Uncomment for patient 1
# tested_treatments <- tested_treatments[1:4039, ] #_1_2

# # Uncomment for patient 2
# tested_treatments <- tested_treatments[1:4062, ] #_2_2

# # Uncomment for patient 3
# tested_treatments <- tested_treatments[1:4040, ] #_3_2

# # Uncomment for patient 4
# tested_treatments <- tested_treatments[1:4039, ] #_4_2


############################################################
# 3) Prepare BEST treatment with numeric axes
############################################################

best_treatment <- best_treatment %>%
  mutate(
    EZHi_total = round(Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3, 0),
    AKTi_total = round(Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3, 0),
    AKTi_days_total = AKTi_Days_Week1 + AKTi_Days_Week2 + AKTi_Days_Week3
  )

############################################################
# 4) FIX: ensure best_treatment is inside plotting range
############################################################

# Update the total

tested_treatments <- tested_treatments %>%
  mutate(
    EZHi_total = 2*Dose_EZHi_pt*5 + 2*Dose_EZHi_w1*7 + 2*Dose_EZHi_w2*7 + 2*Dose_EZHi_w3*7,
    AKTi_total = (2*Dose_AKTi_w1 * AKTi_Days_Week1) +
      (2*Dose_AKTi_w2 * AKTi_Days_Week2) +
      (2*Dose_AKTi_w3 * AKTi_Days_Week3)
  )


best_treatment <- best_treatment %>%
  mutate(
    EZHi_total = 2*Dose_EZHi_pt*5 + 2*Dose_EZHi_w1*7 + 2*Dose_EZHi_w2*7 + 2*Dose_EZHi_w3*7,
    AKTi_total = (2*Dose_AKTi_w1 * AKTi_Days_Week1) +
      (2*Dose_AKTi_w2 * AKTi_Days_Week2) +
      (2*Dose_AKTi_w3 * AKTi_Days_Week3)
  )



# Extend axis limits if needed
x_min <- min(tested_treatments$EZHi_total, best_treatment$EZHi_total)
x_max <- max(tested_treatments$EZHi_total, best_treatment$EZHi_total)

y_min <- min(tested_treatments$AKTi_total, best_treatment$AKTi_total)
y_max <- max(tested_treatments$AKTi_total, best_treatment$AKTi_total)

############################################################
# 5) FINAL PLOT (numeric axes)
############################################################


ggplot(tested_treatments, aes(x = EZHi_total, y = AKTi_total, color = -Cost_Function)) +
  geom_point(size = 5) +
  scale_color_gradientn(
    # colours = c("#d7191c","#ffe08a","#2c7bb6"),
    colours = c("#d7191c","gold","#2c7bb6"),
    # values  = c(0, 0.9, 1),   # <-- white at 80% of the scale #_1
    values  = c(0, 0.8, 1),   # <-- white at 80% of the scale
    trans = scales::pseudo_log_trans(sigma = 0.02),
    guide = guide_colorbar(
      barheight = unit(4, "cm"),
      barwidth  = unit(0.4, "cm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  )+
labs(
  title = "Genetic Algorithm Optimization Results",
  subtitle = "Cost Function (negative = better) vs Total Drug Exposure",
  x = "Total EZH2i Exposure (mg over full schedule)",
  y = "Total AKTi Exposure (mg over full schedule)",
  color = "- Cost Function"
) +
  scale_y_continuous(limits = c(0, 10500)) +
  theme_minimal() +
  geom_point(
    data = best_treatment,
    aes(x = EZHi_total, y = AKTi_total),
    color = "darkorange",
    size = 8,
    shape = 22
  ) +
  annotate(
    "text",
    x = best_treatment$EZHi_total,
    y = best_treatment$AKTi_total,
    label = "Optimal schedule",
    color = "darkorange",
    vjust = -1.5,
    size = 5
  )+ theme(
    text = element_text(family = "Arial"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# ################################################################################
# ################################################################################
# ################################################################################
# # NEW SECTION — RE-EVALUATE SELECTED SCHEDULES FOR WATERFALL PLOT
# ################################################################################
#

ordered_idx <- order(tested_treatments$Cost_Function, decreasing = FALSE)
M <- length(ordered_idx)

# Optimal schedule - patient level
 selected_schedule_indices <- c(ordered_idx[1]) #_2

 selected_schedules <- tested_treatments[selected_schedule_indices, ]


# Optimal schedule - population level
 selected_schedules <- bind_rows(
   selected_schedules,
   tibble(
     Dose_EZHi_pt     = 608.0943,
     Dose_EZHi_w1     = 671.5265,
     Dose_EZHi_w2     = 543.2221,
     Dose_EZHi_w3     = 468.53446,
     Dose_AKTi_w1     = 397.2609,
     Dose_AKTi_w2     = 265.1715,
     Dose_AKTi_w3     = 198.4646,
     AKTi_Days_Week1  = 6,
     AKTi_Days_Week2  = 6,
     AKTi_Days_Week3  = 4,
     Y1_T624          = 8.632480,
     Cost_Function    = -0.08632480,
     EZHi_total       =  29646.91,
     AKTi_total       = 9536.906,
     AKTi_days_total  = 16
   )
 )

 selected_schedule_indices <- c(1,2)


###############################################################################
# BOX PLOTS FOR SELECTED SCHEDULES
###############################################################################


set.seed(1234) # Ensure reproducibility


waterfall_data <- selected_schedules %>%
  mutate(schedule_id = row_number())



schedule_info <- selected_schedules %>%
  mutate(schedule_id = row_number()) %>%
  mutate(
    EZHi_total = Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3,
    AKTi_total = Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3
  ) %>%
  group_by(schedule_id) %>%
  summarise(
    avg_cost = mean(Cost_Function),
    EZHi_total = mean(EZHi_total),
    AKTi_total = mean(AKTi_total),
    .groups = "drop"
  ) %>%
  arrange(avg_cost) %>%                     # order schedules by cost (best → worst)
  mutate(
    legend_label = paste0(
      "ID ", schedule_id,
      " | cost=", round(avg_cost, 3),
      " | EZHi=", round(EZHi_total),
      " | AKTi=", round(AKTi_total)
    ),
    legend_label = factor(legend_label,     # enforce legend order by avg_cost
                          levels = legend_label)
  )



cat("Generating point plot of y1(T=624) for selected schedules...\n")

# Merge existing schedule_info into boxplot_df
boxplot_df <- waterfall_data %>%
  dplyr::select(schedule_id, Y1_T624) %>%
  left_join(schedule_info, by = "schedule_id")   # <-- USE SAME schedule_info

# Reorder schedules by average cost (ascending = best first)
boxplot_df <- boxplot_df %>%
  arrange(avg_cost)

boxplot_df$legend_label <- factor(
  boxplot_df$legend_label,
  levels = unique(boxplot_df$legend_label)
)

# Plot
ggplot(
  boxplot_df,
  aes(
    x = legend_label,
    y = Y1_T624,
    color = legend_label
  )
) +

  # Plot one point per schedule
  geom_point(
    size = 6.5
  ) +


  theme_minimal() +

  labs(
    title = "y1(T=624) for Selected Schedules",
    subtitle = paste("Digital twin prediction for one patient"),
    # x = "Schedule (ID | AvgCost | EZHi | AKTi)",
    # y = "% Dying Cells (y1 at T=624)",
    color = "Schedule (AvgCost | EZHi | AKTi)"
  ) +

  scale_color_viridis_d(
    option = "plasma",
    begin = 0.05,
    end = 0.9
  ) +

  scale_y_continuous(limits = c(8.5, 10.5)) +
  theme(
    axis.text.x = element_blank(),    # <-- hides the long schedule labels
    axis.ticks.x = element_blank(),   # <-- optional: hides ticks too
    text = element_text(family = "Arial"),
    # axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 1),
    axis.title.y = element_text(size = 1)
  )




###############################################
# Detailed study of weekly doses
###############################################
###############################################
# DEFINE SCHEDULES OF INTEREST
###############################################


set.seed(1234) # Ensure reproducibility


df_sub <- selected_schedules %>%
  mutate(schedule_id = row_number())


###############################################
# ORDER BY COST FUNCTION (best = lowest cost)
###############################################
df_sub <- df_sub %>%
  arrange(Cost_Function) %>%
  mutate(order_rank = row_number(),
         legend_label = paste0("ID ", schedule_id,
                               " | cost=", round(Cost_Function, 3)))

# enforce ordered legend
df_sub$legend_label <- factor(df_sub$legend_label,
                              levels = df_sub$legend_label)

###############################################
# BUILD LONG DATA FOR AKTi — WEEKLY TOTAL DOSE
###############################################

akti_long <- df_sub %>%
  dplyr::select(schedule_id, legend_label,
         AKTi_w1 = Dose_AKTi_w1, AKTi_d1 = AKTi_Days_Week1,
         AKTi_w2 = Dose_AKTi_w2, AKTi_d2 = AKTi_Days_Week2,
         AKTi_w3 = Dose_AKTi_w3, AKTi_d3 = AKTi_Days_Week3) %>%

  # compute weekly totals only (no pretreat for AKTi)
  mutate(
    AKTi_week1_total = 2*AKTi_w1 * AKTi_d1,
    AKTi_week2_total = 2*AKTi_w2 * AKTi_d2,
    AKTi_week3_total = 2*AKTi_w3 * AKTi_d3
  ) %>%

  # select only week 1–3 totals
  dplyr::select(schedule_id, legend_label,
         `Week 1` = AKTi_week1_total,
         `Week 2` = AKTi_week2_total,
         `Week 3` = AKTi_week3_total) %>%

  pivot_longer(cols = c(`Week 1`, `Week 2`, `Week 3`),
               names_to = "timepoint",
               values_to = "dose") %>%

  mutate(timepoint = factor(timepoint,
                            levels = c("Week 1", "Week 2", "Week 3")))


###############################################
# BUILD LONG DATA FOR EZH2i — WEEKLY TOTAL DOSE
###############################################

ezh2_long <- df_sub %>%
  dplyr::select(schedule_id, legend_label,
         EZ_pt = Dose_EZHi_pt,
         EZ_w1 = Dose_EZHi_w1,
         EZ_w2 = Dose_EZHi_w2,
         EZ_w3 = Dose_EZHi_w3) %>%

  mutate(
    Pretreat_total = 2*EZ_pt * 5,   # 5 pretreat days
    EZ_week1_total = 2*EZ_w1 * 7,
    EZ_week2_total = 2*EZ_w2 * 7,
    EZ_week3_total = 2*EZ_w3 * 7
  ) %>%

  dplyr::select(schedule_id, legend_label,
         Pretreat = Pretreat_total,
         `Week 1` = EZ_week1_total,
         `Week 2` = EZ_week2_total,
         `Week 3` = EZ_week3_total) %>%

  pivot_longer(cols = c(Pretreat, `Week 1`, `Week 2`, `Week 3`),
               names_to = "timepoint",
               values_to = "dose") %>%

  mutate(timepoint = factor(timepoint,
                            levels = c("Pretreat", "Week 1", "Week 2", "Week 3")))


###############################################
# PLOT 1 — AKTi schedule over time
###############################################
ggplot(akti_long,
       aes(x = timepoint, y = dose,
           group = legend_label, color = legend_label)) +
  geom_line(size = 1.3) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  labs(
    title = "AKTi Weekly Total Dose (Week 1 → Week 3)",
    x = "Treatment Week",
    y = "Weekly AKTi Dose (mg)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
  text = element_text(family = "Arial"),
  axis.text.x  = element_text(size = 10, angle = 30, hjust = 1),
  axis.text.y  = element_text(size = 10),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12)
)



###############################################
# PLOT 2 — EZH2i schedule over time
###############################################
ggplot(ezh2_long,
       aes(x = timepoint, y = dose, group = legend_label, color = legend_label)) +
  geom_line(size = 1.3) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  labs(
    title = "EZH2i Dosing Schedule (Pretreat → Week 3)",
    x = "Timepoint",
    y = "EZH2i Dose (mg)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x  = element_text(size = 10, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

####################################################################################

###############################################
# COMPUTE TOTAL DOSE PER DRUG FOR EACH SCHEDULE
###############################################

df_totals <- df_sub %>%
  mutate(
    EZ_total = 2*Dose_EZHi_pt*5 + 2*Dose_EZHi_w1*7 + 2*Dose_EZHi_w2*7 + 2*Dose_EZHi_w3*7,
    AKT_total = (2*Dose_AKTi_w1 * AKTi_Days_Week1) +
      (2*Dose_AKTi_w2 * AKTi_Days_Week2) +
      (2*Dose_AKTi_w3 * AKTi_Days_Week3)
  ) %>%
  # Long format: one row per drug per schedule
  dplyr::select(schedule_id, legend_label, EZ_total, AKT_total) %>%
  pivot_longer(cols = c(EZ_total, AKT_total),
               names_to = "drug",
               values_to = "total_dose") %>%
  mutate(
    drug = factor(drug,
                  levels = c("AKT_total", "EZ_total"),
                  labels = c("AKTi", "EZH2i"))
  )


ggplot(df_totals,
       aes(x = drug, y = total_dose, color = legend_label, fill = legend_label)) +
  geom_point(size = 4, shape = 22) +   # filled square in schedule color
  scale_color_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  labs(
    title = "Total Drug Dose per Schedule",
    x = "Drug",
    y = "Total Dose Over Pretreatment + 3 Weeks (mg)",
    color = "Schedule",
    fill  = "Schedule"
  ) +
  theme_minimal(base_size = 14)


df_totals %>%
  filter(drug == "AKTi") %>%
  ggplot(aes(x = drug, y = total_dose,
             color = legend_label, fill = legend_label)) +
  # geom_point(size = 4, shape = 22) +   # filled square
  geom_point(size = 3, shape = 0, stroke = 1.5) +   # empty square
  scale_color_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  labs(
    title = "Total AKTi Dose per Schedule",
    x = "Drug",
    y = "Total Dose (mg)"
  ) +
  theme_minimal(base_size = 14)+
  coord_cartesian(clip = "off") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x  = element_text(size = 10, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )




df_totals %>%
  filter(drug == "EZH2i") %>%
  ggplot(aes(x = drug, y = total_dose,
             color = legend_label, fill = legend_label)) +
  # geom_point(size = 4, shape = 22) +   # filled square
  geom_point(size = 3, shape = 0, stroke = 1.5) +   # empty square
  scale_color_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  labs(
    title = "Total EZH2i Dose per Schedule",
    x = "Drug",
    y = "Total Dose (mg)"
  ) +
  theme_minimal(base_size = 14)+
  coord_cartesian(clip = "off") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x  = element_text(size = 10, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )