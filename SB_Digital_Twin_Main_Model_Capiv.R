

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


# Draw parameter samples and run the ODE for each sample
n_patients <- 500 #500

# Extract posterior samples as a data frame (assuming posterior_samples is a list or matrix)
posterior_samples_df <- as.data.frame(posterior_samples)

# Initialize a data frame to store parameters for each in-silico patient
params_samples <- data.frame(matrix(ncol = 30, nrow = n_patients))
colnames(params_samples) <- paste0("k", 1:30)

# Loop to draw random parameter samples for each in-silico patient
for (i in 1:n_patients) {
  # Sample parameter values from the posterior
  sample_idx <- sample(1:nrow(posterior_samples_df), 1)
  theta_sample <- posterior_samples_df[sample_idx, ]
  
  # Store the sampled parameters in the data frame
  params_samples[i, ] <- theta_sample[1:30]
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

CV_Vd <- 0.46
CV_CL <- 0.49

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
# Generate patient-specific PK parameters for AKTi / Capivasertib
# ============================================================

# Population values
CL0_pop  <- 62.2      # L/h
V2_pop   <- 47.9      # L
Imax_pop <- -1.54     # dimensionless
ka_pop   <- 0.417     # 1/h
t0_pop   <- 45.1      # h

# Coefficients of variation (expressed as fractions)
CV_CL0  <- 0.393
CV_V2   <- 1.14
CV_Imax <- 0.706
CV_ka   <- 0.15
CV_t0   <- 0.15

# Log-normal standard deviations
omega_CL0  <- sqrt(log(1 + CV_CL0^2))
omega_V2   <- sqrt(log(1 + CV_V2^2))
omega_Imax <- sqrt(log(1 + CV_Imax^2))
omega_ka   <- sqrt(log(1 + CV_ka^2))
omega_t0   <- sqrt(log(1 + CV_t0^2))

# ============================================================
# Storage
# ============================================================

pk_Capivasertib_params_samples <- data.frame(
  CL0  = numeric(n_patients),
  V2   = numeric(n_patients),
  Imax = numeric(n_patients),
  ka   = numeric(n_patients),
  t0   = numeric(n_patients)
)

# ============================================================
# Generate one PK parameter set per virtual patient
# ============================================================

for(i in 1:n_patients){
  
  CL0_i  <- CL0_pop  * exp(rnorm(1, mean = 0, sd = omega_CL0))
  V2_i   <- V2_pop   * exp(rnorm(1, mean = 0, sd = omega_V2))
  Imax_i <- Imax_pop * exp(rnorm(1, mean = 0, sd = omega_Imax))
  ka_i   <- ka_pop   * exp(rnorm(1, mean = 0, sd = omega_ka))
  t0_i   <- t0_pop   * exp(rnorm(1, mean = 0, sd = omega_t0))
  
  pk_Capivasertib_params_samples[i, ] <- c(
    CL0_i,
    V2_i,
    Imax_i,
    ka_i,
    t0_i
  )
}


################################################################################

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
    A_GI = D0_e[1]*0.33,   # Initial amount of drug in GI tract * Bioavailability
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
      A_GI = current_state_first[1],   # Initial amount of drug in GI tract (mg)
      C_plasma = current_state_first[2]  # Initial concentration in plasma (mg/L)
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
    CLI  = as.numeric(pk_params$CL0),
    V2I  = as.numeric(pk_params$V2),
    V3I = 113,
    V4I = 94.7,
    Q3I = 2.66,
    Q4I = 21.8,
    kaI  = as.numeric(pk_params$ka),
    Imax = as.numeric(pk_params$Imax),
    
    t0   = as.numeric(pk_params$t0),             # zero-order duration
    T50 = 67.4,
    Imax_dose = -0.00183,
    dose_times = dose_times,
    D0_a = D0_a
  )

  
  # ============================================================
  # INITIAL CONDITIONS
  # ============================================================
  
  state <- c(
    A_depot_1    = 0,   # zero-order depot
    A_depot_2    = 0,   # first-order depot
    A_central    = 0,
    A_peripheral1 = 0,
    A_peripheral2 = 0
  )
  
  
  # ============================================================
  # ODE MODEL
  # ============================================================
  
  pk_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      
      # ---------------------------------------
      # ZERO-ORDER INPUT (20% of each dose)
      # ---------------------------------------
      R_zero <- 0
      current_dose <- 0
      
      for (i in seq_along(dose_times)) {
        
        t_i <- dose_times[i]
        D_i <- D0_a[i]
        
        if (time >= t_i && D_i > 0) {
          current_dose <- D_i
        }
        
        if (time >= t_i && time < t_i + t0) {
          R_zero <- R_zero + (0.20 * D_i / t0)
        }
        
      }
      
      Imax_D <- Imax * (1 + (current_dose - 480) * Imax_dose)
      
      CL_t <- CLI * ( 1 - (exp(Imax_D) * time^5) / (T50^5 + time^5))
      
      # ---------------------------------------
      # ODEs
      # ---------------------------------------
      dA_depot_1 <- -R_zero                  # depot1 only releases
      
      dA_depot_2 <- -kaI * A_depot_2        # absorbs first-order
      
      dA_central <-
        kaI * A_depot_2 +                   # first-order absorption
        R_zero -                            # zero-order absorption
        (CL_t / V2I) * A_central -
        (Q3I / V2I) * A_central + (Q3I / V3I) * A_peripheral1 -
        (Q4I / V2I) * A_central + (Q4I / V4I) * A_peripheral2
      
      dA_peripheral1 <-(Q3I / V2I) * A_central - (Q3I / V3I) * A_peripheral1
      
      dA_peripheral2 <- (Q4I / V2I) * A_central - (Q4I / V4I) * A_peripheral2
      
      list(c(dA_depot_1, dA_depot_2, dA_central, dA_peripheral1, dA_peripheral2))
    })
  }
  
  
  # ============================================================
  # EVENT TABLE (80% into depot2, 20% into depot1)
  # ============================================================
  
  event_table <- data.frame(
    var  = rep(c("A_depot_1", "A_depot_2"), each = length(dose_times)),
    time = rep(dose_times, times = 2),
    value = c(0.20 * D0_a, 0.80 * D0_a),
    method = "add"
  )
  
  event_table <- event_table[order(event_table$time), ]
  
  
  # ============================================================
  # SOLVE SYSTEM
  # ============================================================
  
  times <- seq(0, total_time, by = times_between_doses)
  
  output_df_AKTi <- as.data.frame(ode(
    y = state,
    times = times,
    func = pk_model,
    parms = parameters,
    events = list(data = event_table),
    method = "bdf",
    maxsteps = 5000
  ))
  
  
  # ============================================================
  # CONCENTRATIONS
  # ============================================================
  
  output_df_AKTi$C_central <- output_df_AKTi$A_central / parameters$V2I
  
  MW <- 428.9
  output_df_AKTi <- output_df_AKTi %>% mutate(C_central = (C_central * 10000) / MW)

  
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

# Storage for results
results_storage <- list(y1_values = numeric(), cost_values = numeric(), D0_e_pt = numeric(), D0_e_w1 = numeric(), D0_e_w2 = numeric(), D0_e_w3 = numeric(), D0_a_w1 = numeric(), D0_a_w2 = numeric(), D0_a_w3 = numeric(), AKTi_Days_Week1 = numeric(), AKTi_Days_Week2 = numeric(), AKTi_Days_Week3 = numeric())



evaluation_counter <- 0  # Global counter

# GA Objective Function
ga_objective <- function(x, limit_toxicity_e, limit_toxicity_a) {
  
  D0_e_day <- generate_dosing_schedule_D0e(c(x[1], x[2], x[3], x[4]))
  D0_a_day <- generate_dosing_schedule_D0a(c(x[5], x[6], x[7]), c(x[8], x[9], x[10]))
  
  # New: expand to two doses/day
  D0_e <- expand_daily_doses(D0_e_day)
  D0_a <- expand_daily_doses(D0_a_day)
  
  
  dose_times <- sort(c(
    1,      1    + 12,
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
  times_between_doses <- 0.1
  
  # Solve the PK models for AKTi and EZHi
  output_df_AKTi <- define_PK_AKTi_fun(D0_a, dose_times, total_time, times_between_doses, pk_params = pk_Capivasertib_params_samples[i, ])
  
  output_df_EZHi <- define_PK_EZHi_fun(D0_e, dose_times[-1], total_time, times_between_doses, pk_params = pk_params_samples[i, ])
  
  
  # Ensure unique time points
  pk_AKTi_unique <- output_df_AKTi[!duplicated(output_df_AKTi[, "time"]), ]
  pk_EZHi_unique <- output_df_EZHi[!duplicated(output_df_EZHi[, "time"]), ]
  
  # Interpolating functions for concentrations
  AKTi_interp <- approxfun(pk_AKTi_unique[, "time"], pk_AKTi_unique[, "C_central"], rule = 2)
  EZHi_interp <- approxfun(pk_EZHi_unique[, "time"], pk_EZHi_unique[, "C_plasma"], rule = 2)
  
  # Define the time sequence for the PD model simulation
  times <- seq(0, 624, by = 0.1)
  
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
  
  
  results_list <- future_lapply(1:n_patients, function(i) {
    
    params <- as.list(params_samples[i, ])
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
  
  maxiter = 70,
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
# A) EZH2 ONLY (AKTi = 0)
# =====================================================
for (dose in ezh2_doses) {

  x_manual <- c(
    dose, dose, dose, dose,   # EZH2 dose
    0, 0, 0,      # AKTi dose
    0, 0, 0 # no AKTi schedule
  )

  -ga_objective(x_manual, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)
}

# =====================================================
# B) AKTi ONLY (EZH2 = 0)
# =====================================================

for (dose in akti_doses) {

  # ---- Daily AKTi : 7-7-7 ----
  x_manual_daily <- c(
    0, 0, 0, 0,      # EZH2 dose
    dose, dose, dose,   # AKTi dose
    7,7,7   # AKTi schedule (daily)
  )
  -ga_objective(x_manual_daily, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)

  # ---- AKTi 4-4-4 ----
  x_manual_444 <- c(
    0, 0, 0, 0,      # EZH2 dose
    dose, dose, dose,   # AKTi dose
    4,4,4   # AKTi schedule (4-4-4)
  )
  -ga_objective(x_manual_444, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)
}



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
# D) effect of EZH2i and AKTi
# =====================================================

# 5 doses for each drug
ezh2_doses <- seq(lower_bounds[1], upper_bounds[1], length.out = 5)
akti_doses <- seq(lower_bounds[5], upper_bounds[5], length.out = 5)


akti_doses <- c(5, 10, 20, 40, 80)
ezh2_doses <- c(5, 10, 20, 40, 80)

dos_e <- 800
dos_a <- 400

for (dose in ezh2_doses) {

  x_manual_ezh2_effect <- c(
    dose, dose, dose, dose,   # EZH2 dose
    dos_a, dos_a, dos_a,   # AKTi dose
    4,4,4   # AKTi schedule (4-4-4)
  )

  -ga_objective(x_manual_ezh2_effect, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)
}

for (dose in akti_doses) {

  x_manual_akt_effect <- c(
    dos_e, dos_e, dos_e, dos_e,   # EZH2 dose
    dose, dose, dose,   # AKTi dose
    4,4,4   # AKTi schedule (4-4-4)
  )

  -ga_objective(x_manual_akt_effect, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)
}




# =====================================================
# D) effect of EZH2i and AKTi
# =====================================================

# 5 doses for each drug

akti_doses <- c(1, 100, 200, 300)

dos_e <- 800


for (dose in akti_doses) {

  x_manual_akt_effect <- c(
    dos_e, dos_e, dos_e, dos_e,   # EZH2 dose
    dose, dose, dose,   # AKTi dose
    4,4,4   # AKTi schedule (4-4-4)
  )

  -ga_objective(x_manual_akt_effect, limit_toxicity_e = limit_toxicity_ee, limit_toxicity_a = limit_toxicity_aa)
}


# =====================================================
# Append manual results into tested_treatments
# =====================================================

n_after <- length(results_storage$y1_values)
manual_idx <- (n_before + 1):n_after
# manual_idx <- 1:4

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

# # # Save the objects
# saveRDS(tested_treatments, file = tested_treatments_file)
# saveRDS(best_treatment, file = best_treatment_file)

# # To reload them later
tested_treatments <- readRDS(tested_treatments_file)
best_treatment <- readRDS(best_treatment_file)
# 

set.seed(1234) # Ensure reproducibility


# ###################### SUMMARY PLOTS


tested_treatments <- tested_treatments %>%
  mutate(
    EZHi_total = round(Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3, 0),
    AKTi_total = round(Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3, 0),
    AKTi_days_total = AKTi_Days_Week1 + AKTi_Days_Week2 + AKTi_Days_Week3
  )

############################################################
# 2) Remove NEGATIVE cost-function solutions
############################################################


tested_treatments <- tested_treatments %>%
  filter(Cost_Function <= 0)

tested_treatments <- tested_treatments[1:4034, ] #2

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
  scale_x_continuous(limits = c(0,40000)) +
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




################################################################################
################################################################################
################################################################################
# RE-EVALUATE SELECTED SCHEDULES FOR WATERFALL PLOT
################################################################################

cat("Re-running selected treatment schedules on all patients...\n")

##################################################################
# STEP 1 — SELECT N SCHEDULES (equally spaced + special cases)
##################################################################

set.seed(1234) # Ensure reproducibility

N <- 5

n_patients_waterfall <- 500 # selected for final box plots


# ## General study

selected_schedule_indices <- c(4039, 4040, 4037, 4043, 3871, 2893, 334) #_2
selected_schedules <- tested_treatments[selected_schedule_indices, ]



# # # AKTi study
# 
# selected_schedule_indices_AKTi_study  <- c(163, 169, 175, 176, 177, 178, 179, 180, 181, 182, 183) 
# selected_schedule_indices <- selected_schedule_indices_AKTi_study
# 
# selected_schedules <- tested_treatments[selected_schedule_indices, ]




# ##################################################################
# # Uncomment to make doses discreet
# 
# # ------------------------------------------------------------
# # Round doses to the nearest clinically allowed dose
# # ------------------------------------------------------------
# 
# taz_levels <- c(0, 200, 400, 600, 800)
# ipa_levels <- c(0, 100, 200, 300, 400)
# 
# 
# # CLOEST LEVEL
# round_to_levels <- function(x, levels) {
#   levels[which.min(abs(levels - x))]
# }
# 
# 
# 
# selected_schedules_discrete <- selected_schedules %>%
#   mutate(
#     Dose_EZHi_pt = sapply(Dose_EZHi_pt, round_to_levels, levels = taz_levels),
#     Dose_EZHi_w1 = sapply(Dose_EZHi_w1, round_to_levels, levels = taz_levels),
#     Dose_EZHi_w2 = sapply(Dose_EZHi_w2, round_to_levels, levels = taz_levels),
#     Dose_EZHi_w3 = sapply(Dose_EZHi_w3, round_to_levels, levels = taz_levels),
# 
#     Dose_AKTi_w1 = sapply(Dose_AKTi_w1, round_to_levels, levels = ipa_levels),
#     Dose_AKTi_w2 = sapply(Dose_AKTi_w2, round_to_levels, levels = ipa_levels),
#     Dose_AKTi_w3 = sapply(Dose_AKTi_w3, round_to_levels, levels = ipa_levels)
#   )
# 
# selected_schedules <- selected_schedules_discrete
# ##################################################################


#################################################################
STEP 2 — SIMULATION FUNCTION
#################################################################


set.seed(1234) # Ensure reproducibility

simulate_schedule_for_patient <- function(D0_e_single, D0_a_single, days_on_vec, patient_params, pk_akti_params, pk_ezh2i_params, dt) {

  D0_e <- expand_daily_doses(generate_dosing_schedule_D0e(D0_e_single))
  D0_a <- expand_daily_doses(generate_dosing_schedule_D0a(D0_a_single, days_on_vec))

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



  pkA <- define_PK_AKTi_fun(D0_a, dose_times, total_time, dt, pk_params = pk_akti_params)

  pkE <- define_PK_EZHi_fun(D0_e, dose_times[-1], total_time, dt, pk_params = pk_ezh2i_params)

  pkA <- pkA[!duplicated(pkA$time), ]
  pkE <- pkE[!duplicated(pkE$time), ]

  AKTi_interp <- approxfun(pkA$time, pkA$C_central, rule = 2)
  EZHi_interp <- approxfun(pkE$time, pkE$C_plasma, rule = 2)

  times <- seq(0, 624, by = dt)
  y_initial <- c(y1=10.6765,y2=0,y3=0,y4=1,y5=1,y6=0,y7=0,y8=0)







  PD_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      k31 <- AKTi_interp(time)
      k32 <- EZHi_interp(time)

      dy1_dt <- (k1*(y6/(k2+y6))*exp(-k16*time)*y8) - k15*y1*(1-(y1/100))
      dy2_dt <- k3*(k4/(k4+y4))*y3 - k5*y2
      dy3_dt <- (k6+k25*y3)*(50-y3) - (k29+k7*(y5/(k8+y5))*(1+k26*(50-y3)))*y3
      dy4_dt <- k9 - k10*(k31/(k11+k31))*y4 - k5*y4
      dy5_dt <- k12 - k13*(k32/(k14+k32))*y5 - k5*y5
      dy6_dt <- k19*(k20/(k20+y4))*y7 - k21*y6
      dy7_dt <- (k22+k27*y7)*(50-y7) -
        (k30+k23*(y5/(k24+y5))*(1+k28*(50-y7)))*y7
      dy8_dt <- k17*(y2/(k18+y2))*(100-y1-y8) -
        (k1*(y6/(k2+y6))*exp(-k16*time) - k15*(y1/100))*y8

      list(c(dy1_dt,dy2_dt,dy3_dt,dy4_dt,dy5_dt,dy6_dt,dy7_dt,dy8_dt))
    })
  }

  pd_output <- ode(y = y_initial, times = times, func = PD_model, parms = patient_params)

  return(tail(pd_output[, "y1"], 1))
}


##################################################################
# STEP 3 — RUN ALL SELECTED SCHEDULES (PARALLEL)
##################################################################

patient_indices_waterfall <- sample(1:n_patients, n_patients_waterfall)

waterfall_data <- data.frame()
NN <- length(selected_schedule_indices)

for (s in seq_len(NN)) {

  row <- selected_schedules[s, ]

  D0_e_single <- c(row$Dose_EZHi_pt, row$Dose_EZHi_w1, row$Dose_EZHi_w2, row$Dose_EZHi_w3)
  D0_a_single <- c(row$Dose_AKTi_w1, row$Dose_AKTi_w2, row$Dose_AKTi_w3)
  days_on_vec <- c(row$AKTi_Days_Week1, row$AKTi_Days_Week2, row$AKTi_Days_Week3)



  cat("\n========================\n")
  cat("Schedule", s, "\n")
  print(D0_e_single)
  print(D0_a_single)
  cat("========================\n")


  ###PARALLELIZED
  results_list <- future_lapply(patient_indices_waterfall, function(i) {

    ###NOT PARALLELIZED
    # results_list <- lapply(patient_indices_waterfall,function(i) {

    patient_par <- as.list(params_samples[i, ])

    pk_akti_par <- pk_Capivasertib_params_samples[i, ]

    pk_ezh2i_par <- pk_params_samples[i, ]

    dt_i <- 0.1



    y1_val <- simulate_schedule_for_patient(D0_e_single, D0_a_single, days_on_vec, patient_par, pk_akti_par, pk_ezh2i_par, dt_i)

    cost_val <- -(y1_val / 100)

    data.frame(
      schedule_rank = s,
      schedule_id = selected_schedule_indices[s],
      patient_id = i,
      y1_T624 = y1_val,
      cost = cost_val
    )
  })

  waterfall_data <- rbind(waterfall_data, do.call(rbind, results_list))
}


##################################################################
# STEP 4 — BUILD waterfall_cost (NOW EXISTS)
##################################################################

waterfall_cost <- waterfall_data %>%
  arrange(cost) %>%
  mutate(sim_index = row_number())


##################################################################
# STEP 5 — BUILD schedule_info FOR LEGEND
##################################################################

schedule_info <- tested_treatments %>%
  mutate(schedule_id = row_number()) %>%
  filter(schedule_id %in% unique(waterfall_cost$schedule_id)) %>%
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

waterfall_cost <- waterfall_cost %>%
  left_join(schedule_info, by = "schedule_id")


############################################################
# SAVE WATERFALL / SCHEDULE SUMMARY OBJECTS
############################################################

# # Main study
# schedule_info_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/06_2026_CAPIV_main_study_schedule_info.rds"
#
# waterfall_data_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/06_2026_CAPIV_main_study_waterfall_data.rds"
#
# # AKTi study
# schedule_info_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/06_2026_CAPIV_AKTi_study_schedule_info_1.rds"
#
# waterfall_data_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/06_2026_CAPIV_AKTi_study_waterfall_data_1.rds"
#
# # Main study - DISCREET - CLOSEST
schedule_info_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/06_2026_CAPIV_main_study_DISCREET_CLOSEST_schedule_info_1.rds"

waterfall_data_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/06_2026_CAPIV_main_study_DISCREET_CLOSEST_waterfall_data_1.rds"



# # Main study
# schedule_info_file <- "/Users/.../Capiv_main_study_schedule_info.rds"
# 
# waterfall_data_file <- "/Users/.../Capiv_main_study_waterfall_data.rds"

# # Main study - DISCREET
# schedule_info_file <- "/Users/.../Capiv_main_study_discreet_schedule_info.rds"
# 
# waterfall_data_file <- "/Users/.../Capiv_main_study_discreet_waterfall_data.rds"

# # # AKTi study
# schedule_info_file <- "/Users/.../Capiv_AKTi_study_schedule_info.rds"
# 
# waterfall_data_file <- "/Users/.../Capiv_AKTi_study_waterfall_data.rds"

# Save objects
# saveRDS(schedule_info,  file = schedule_info_file)
# saveRDS(waterfall_data, file = waterfall_data_file)


# 
# # # To reload them later
schedule_info <- readRDS(schedule_info_file)
waterfall_data <- readRDS(waterfall_data_file)
# #



###############################################################################
#  BOX PLOTS FOR SELECTED SCHEDULES (using waterfall_data)
###############################################################################


set.seed(1234) # Ensure reproducibility

cat("Generating boxplots of y1(T=624) for selected schedules...\n")

# Merge existing schedule_info into boxplot_df
boxplot_df <- waterfall_data %>%
  dplyr::select(schedule_id, y1_T624) %>%
  left_join(schedule_info, by = "schedule_id")   # <-- USE SAME schedule_info


##########################################################
##### Uncomment for AKTi study

# # Schedule IDs to keep
# selected_schedule_ids <- c(163, 180,181, 182, 183, 169)
#
# boxplot_df <- waterfall_data %>%
#   dplyr::select(schedule_id, y1_T624) %>%
#   left_join(schedule_info, by = "schedule_id") %>%
#   filter(schedule_id %in% selected_schedule_ids)
##########################################################


# Reorder schedules by average cost (ascending = best first)
boxplot_df <- boxplot_df %>%
  arrange(avg_cost)

boxplot_df$legend_label <- factor(boxplot_df$legend_label,
                                  levels = unique(boxplot_df$legend_label))

# Plot
ggplot(boxplot_df,
       aes(x = legend_label, y = y1_T624, fill = legend_label)) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Box Plot of y1(T=624) for Selected Schedules",
    subtitle = paste("Based on", n_patients_waterfall, "patients per schedule"),
    x = "Schedule (ID | AvgCost | EZHi | AKTi)",
    y = "% Dying Cells (y1 at T=624)",
    fill = "Schedule (AvgCost | EZHi | AKTi)"
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  # scale_fill_brewer(
  #   palette = "Blues",   # light → dark BLUE (AKTi study)
  #   direction = -1         # higher factor levels = darker
  # ) +
  theme(
    axis.text.x = element_blank(),   # <-- hides the long schedule labels
    axis.ticks.x = element_blank(),   # <-- optional: hides ticks too
    text = element_text(family = "Arial"),
    # axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


# ######################################################################
# ######################################################################

###############################################
# Detailed study of weekly doses
###############################################
###############################################
# DEFINE SCHEDULES OF INTEREST
###############################################
set.seed(1234) # Ensure reproducibility

sched_ids  <- c(4039, 4040, 4037, 4043, 3871, 2893, 334) #_2

df_sub <- tested_treatments %>%
  mutate(schedule_id = row_number()) %>%
  filter(schedule_id %in% sched_ids)

# # ##############################################################
# ######### FOR DISCREET DOSES
# 
# # ------------------------------------------------------------
# # Round doses to the nearest clinically allowed dose
# # ------------------------------------------------------------
# 
# taz_levels <- c(0, 200, 400, 600, 800)
# ipa_levels <- c(0, 100, 200, 300, 400)
# 
# round_to_levels <- function(x, levels) {
#   levels[which.min(abs(levels - x))]
# }
# 
# 
# df_sub <- df_sub %>%
#   mutate(
#     Dose_EZHi_pt = sapply(Dose_EZHi_pt, round_to_levels, levels = taz_levels),
#     Dose_EZHi_w1 = sapply(Dose_EZHi_w1, round_to_levels, levels = taz_levels),
#     Dose_EZHi_w2 = sapply(Dose_EZHi_w2, round_to_levels, levels = taz_levels),
#     Dose_EZHi_w3 = sapply(Dose_EZHi_w3, round_to_levels, levels = taz_levels),
# 
#     Dose_AKTi_w1 = sapply(Dose_AKTi_w1, round_to_levels, levels = ipa_levels),
#     Dose_AKTi_w2 = sapply(Dose_AKTi_w2, round_to_levels, levels = ipa_levels),
#     Dose_AKTi_w3 = sapply(Dose_AKTi_w3, round_to_levels, levels = ipa_levels)
#   )
# 
# # ##############################################################


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


# # ##########################################################################################################################################
# # ##########################################################################################################################################
# # ##########################################################################################################################################
# # ##########################################################################################################################################

# MUTIPLE CYCLE SIMULATIONS


set.seed(1234)

# Number of 21-day treatment cycles
n_cycles <- 2

# Number of virtual patients
n_patients_waterfall <- 500

# Fixed schedule structure
pretreat_days <- 5
cycle_days <- 21

# Evaluation day (measured from simulation start, including pretreatment)
evaluation_day <- 26


################################################################################
# EXTEND ONE-CYCLE DAILY SCHEDULE TO N TREATMENT CYCLES
################################################################################

extend_schedule_to_n_cycles <- function(
    one_cycle_schedule,
    n_cycles,
    pretreat_days = 5
) {
  
  if (!is.numeric(n_cycles) ||
      length(n_cycles) != 1 ||
      n_cycles < 1 ||
      n_cycles != floor(n_cycles)) {
    
    stop("n_cycles must be a positive integer.")
  }
  
  total_days <- length(one_cycle_schedule)
  
  if (total_days <= pretreat_days) {
    stop("Treatment portion not found.")
  }
  
  treatment_days <- total_days - pretreat_days
  
  pretreatment_schedule <- one_cycle_schedule[
    seq_len(pretreat_days)
  ]
  
  treatment_schedule <- one_cycle_schedule[
    (pretreat_days + 1):total_days
  ]
  
  extended_schedule <- c(
    pretreatment_schedule,
    rep(treatment_schedule, n_cycles)
  )
  
  return(extended_schedule)
}



################################################################################
# FUNCTION TO GENERATE BID DOSE TIMES
################################################################################

generate_bid_dose_times <- function(total_days) {
  
  first_daily_dose_times <- 1 + 24 * (0:(total_days - 1))
  second_daily_dose_times <- 13 + 24 * (0:(total_days - 1))
  
  sort(
    c(
      first_daily_dose_times,
      second_daily_dose_times
    )
  )
}



# ============================================================================
# GENERALIZED simulate_schedule_for_patient()
# ============================================================================

simulate_schedule_for_patient <- function(
    D0_e_single,
    D0_a_single,
    days_on_vec,
    patient_params,
    pk_akti_params,
    pk_ezh2i_params,
    dt,
    n_cycles = 1
) {
  
  # --------------------------------------------------------------------------
  # Build one-cycle schedules
  # --------------------------------------------------------------------------
  
  D0_e_daily <- generate_dosing_schedule_D0e(
    D0_e_single,
    pretreat_days = pretreat_days,
    treat_days = cycle_days
  )
  
  D0_a_daily <- generate_dosing_schedule_D0a(
    D0_a_single,
    days_on_vec,
    pretreat_days = pretreat_days,
    treat_days = cycle_days
  )
  
  # --------------------------------------------------------------------------
  # Extend schedules
  # --------------------------------------------------------------------------
  
  D0_e_daily <- extend_schedule_to_n_cycles(
    D0_e_daily,
    n_cycles,
    pretreat_days
  )
  
  D0_a_daily <- extend_schedule_to_n_cycles(
    D0_a_daily,
    n_cycles,
    pretreat_days
  )
  
  # --------------------------------------------------------------------------
  # BID expansion
  # --------------------------------------------------------------------------
  
  D0_e <- expand_daily_doses(D0_e_daily)
  D0_a <- expand_daily_doses(D0_a_daily)
  
  total_days <- length(D0_e_daily)
  total_time <- total_days * 24
  
  dose_times <- generate_bid_dose_times(total_days)
  
  # Sanity check
  
  length(D0_e)
  length(dose_times[-1])
  
  # --------------------------------------------------------------------------
  # PK
  # --------------------------------------------------------------------------
  
  pkA <- define_PK_AKTi_fun(
    D0_a,
    dose_times,
    total_time,
    dt,
    pk_params = pk_akti_params
  )
  
  pkE <- define_PK_EZHi_fun(
    D0_e,
    dose_times[-1],
    total_time,
    dt,
    pk_params = pk_ezh2i_params
  )
  
  pkA <- pkA[!duplicated(pkA$time), ]
  pkE <- pkE[!duplicated(pkE$time), ]
  
  AKTi_interp <- approxfun(pkA$time, pkA$C_central, rule = 2)
  EZHi_interp <- approxfun(pkE$time, pkE$C_plasma, rule = 2)
  
  times <- seq(0, total_time, by = dt)
  
  y_initial <- c(
    y1 = 10.6765,
    y2 = 0,
    y3 = 0,
    y4 = 1,
    y5 = 1,
    y6 = 0,
    y7 = 0,
    y8 = 0
  )
  
  # --------------------------------------------------------------------------
  # COPY YOUR EXISTING PD_model HERE UNCHANGED
  # --------------------------------------------------------------------------
  
  PD_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      
      k31 <- AKTi_interp(time)
      k32 <- EZHi_interp(time)
      
      dy1_dt <- (k1*(y6/(k2+y6))*exp(-k16*time)*y8) - k15*y1*(1-(y1/100))
      dy2_dt <- k3*(k4/(k4+y4))*y3 - k5*y2
      dy3_dt <- (k6+k25*y3)*(50-y3) - (k29+k7*(y5/(k8+y5))*(1+k26*(50-y3)))*y3
      dy4_dt <- k9 - k10*(k31/(k11+k31))*y4 - k5*y4
      dy5_dt <- k12 - k13*(k32/(k14+k32))*y5 - k5*y5
      dy6_dt <- k19*(k20/(k20+y4))*y7 - k21*y6
      dy7_dt <- (k22+k27*y7)*(50-y7) -
        (k30+k23*(y5/(k24+y5))*(1+k28*(50-y7)))*y7
      dy8_dt <- k17*(y2/(k18+y2))*(100-y1-y8) -
        (k1*(y6/(k2+y6))*exp(-k16*time) - k15*(y1/100))*y8
      
      list(c(dy1_dt,dy2_dt,dy3_dt,dy4_dt,dy5_dt,dy6_dt,dy7_dt,dy8_dt))
      
    })
  }
  
  pd_output <- ode(
    y = y_initial,
    times = times,
    func = PD_model,
    parms = patient_params
  )
  
  # --------------------------------------------------------------------------
  # Compute additional metrics
  # --------------------------------------------------------------------------
  
  pd_output <- as.data.frame(pd_output)
  
  death_flux <- with(
    pd_output,
    patient_params$k1 *
      (y6 / (patient_params$k2 + y6)) *
      exp(-patient_params$k16 * time) *
      y8
  )
  
  pd_output$cumulative_death <- c(
    0,
    cumsum(
      diff(pd_output$time) *
        (head(death_flux, -1) + tail(death_flux, -1)) / 2
    )
  )
  
  pd_output$surviving_cells <- 100 - pd_output$y1
  
  pd_output$survival_auc <- c(
    0,
    cumsum(
      diff(pd_output$time) *
        (
          head(pd_output$surviving_cells, -1) +
            tail(pd_output$surviving_cells, -1)
        ) / 2
    )
  )
  
  # --------------------------------------------------------------------------
  # Evaluate metrics
  # --------------------------------------------------------------------------
  
  evaluation_time <- evaluation_day * 24
  
  interp_metric <- function(variable) {
    approx(
      pd_output$time,
      pd_output[[variable]],
      xout = evaluation_time,
      rule = 2
    )$y
  }
  
  list(
    
    metrics = list(
      
      total_days = total_days,
      total_time = total_time,
      
      y1_evaluation = interp_metric("y1"),
      y1_final = tail(pd_output$y1, 1),
      
      cumulative_death_evaluation = interp_metric("cumulative_death"),
      cumulative_death_final = tail(pd_output$cumulative_death, 1),
      
      surviving_cells_evaluation = interp_metric("surviving_cells"),
      surviving_cells_final = tail(pd_output$surviving_cells, 1),
      
      survival_auc_evaluation = interp_metric("survival_auc"),
      survival_auc_final = tail(pd_output$survival_auc, 1)
      
    )
    
  )
}


##################################################################
# STEP 3 — RUN ALL SELECTED SCHEDULES
##################################################################

patient_indices_waterfall <- sample(1:n_patients, n_patients_waterfall)

waterfall_data <- data.frame()

NN <- length(selected_schedule_indices)

for (s in seq_len(NN)) {
  
  row <- selected_schedules[s, ]
  
  D0_e_single <- c(
    row$Dose_EZHi_pt,
    row$Dose_EZHi_w1,
    row$Dose_EZHi_w2,
    row$Dose_EZHi_w3
  )
  
  D0_a_single <- c(
    row$Dose_AKTi_w1,
    row$Dose_AKTi_w2,
    row$Dose_AKTi_w3
  )
  
  days_on_vec <- c(
    row$AKTi_Days_Week1,
    row$AKTi_Days_Week2,
    row$AKTi_Days_Week3
  )
  
  cat("\n========================\n")
  cat("Schedule", s, "\n")
  print(D0_e_single)
  print(D0_a_single)
  cat("========================\n")
  
  results_list <- future_lapply(patient_indices_waterfall, function(i) {
    
    patient_par <- as.list(params_samples[i, ])
    
    pk_akti_par <- pk_Capivasertib_params_samples[i, ]
    
    pk_ezh2i_par <- pk_params_samples[i, ]
    
    dt_i <- 0.1
    
    
    sim <- simulate_schedule_for_patient(
      D0_e_single      = D0_e_single,
      D0_a_single      = D0_a_single,
      days_on_vec      = days_on_vec,
      patient_params   = patient_par,
      pk_akti_params   = pk_akti_par,
      pk_ezh2i_params  = pk_ezh2i_par,
      dt               = dt_i,
      n_cycles         = n_cycles
    )
    
    y1_final <- sim$metrics$y1_final
    y1_eval  <- sim$metrics$y1_evaluation
    
    cumdeath_final <- sim$metrics$cumulative_death_final
    cumdeath_eval  <- sim$metrics$cumulative_death_evaluation
    
    surviving_final <- sim$metrics$surviving_cells_final
    surviving_eval  <- sim$metrics$surviving_cells_evaluation
    
    auc_final <- sim$metrics$survival_auc_final
    auc_eval  <- sim$metrics$survival_auc_evaluation
    
    cost_val <- -(y1_eval / 100)
    
    data.frame(
      schedule_rank = s,
      schedule_id = selected_schedule_indices[s],
      patient_id = i,
      
      y1_final = y1_final,
      y1_evaluation = y1_eval,
      
      cumulative_death_final = cumdeath_final,
      cumulative_death_evaluation = cumdeath_eval,
      
      surviving_cells_final = surviving_final,
      surviving_cells_evaluation = surviving_eval,
      
      survival_auc_final = auc_final,
      survival_auc_evaluation = auc_eval,
      
      cost = cost_val
    )
    
  })
  
  waterfall_data <- rbind(
    waterfall_data,
    do.call(rbind, results_list)
  )
  
}


##################################################################
# STEP 4 — BUILD waterfall_cost
##################################################################

waterfall_cost <- waterfall_data %>%
  arrange(cost) %>%
  mutate(
    sim_index = row_number()
  )


##################################################################
# STEP 5 — BUILD schedule_info
##################################################################

schedule_info <- tested_treatments %>%
  mutate(schedule_id = row_number()) %>%
  filter(schedule_id %in% unique(waterfall_cost$schedule_id)) %>%
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
  arrange(avg_cost) %>%
  mutate(
    legend_label = paste0(
      "ID ", schedule_id,
      " | cost=", round(avg_cost,3),
      " | EZHi=", round(EZHi_total),
      " | AKTi=", round(AKTi_total)
    ),
    legend_label = factor(
      legend_label,
      levels = legend_label
    )
  )

waterfall_cost <- waterfall_cost %>%
  left_join(schedule_info, by = "schedule_id")


###############################################################################
# STEP 6 — BOX PLOTS FOR SELECTED SCHEDULES (using waterfall_data)
###############################################################################

set.seed(1234) # Ensure reproducibility

cat("Generating boxplots for selected schedules...\n")

plot_y1_boxplot <- function(y_variable, title_suffix, output_suffix) {

  # Merge existing schedule_info into boxplot_df
  boxplot_df <- waterfall_data %>%
    dplyr::select(schedule_id, !!sym(y_variable)) %>%
    left_join(schedule_info, by = "schedule_id")

  # Reorder schedules by average cost
  boxplot_df <- boxplot_df %>%
    arrange(avg_cost)

  boxplot_df$legend_label <- factor(
    boxplot_df$legend_label,
    levels = unique(boxplot_df$legend_label)
  )

  p <- ggplot(
    boxplot_df,
    aes(
      x = legend_label,
      y = .data[[y_variable]],
      fill = legend_label
    )
  ) +
    geom_boxplot(outlier.shape = 16, outlier.alpha = 0.6) +
    theme_minimal() +
    labs(
      title = paste("Box Plot of", title_suffix, "for Selected Schedules"),
      subtitle = paste("Based on", n_patients_waterfall, "patients per schedule"),
      x = "Schedule (ID | AvgCost | EZHi | AKTi)",
      y = title_suffix,
      fill = "Schedule (AvgCost | EZHi | AKTi)"
    ) +
    scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(family = "Arial"),
      axis.text.y  = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )

  print(p)

}

# Evaluation day
plot_y1_boxplot(
  y_variable = "y1_evaluation",
  title_suffix = paste0("% Dying Cells (Evaluation Day = ", evaluation_day, ")"),
  output_suffix = "evaluation"
)

# Final simulation time
plot_y1_boxplot(
  y_variable = "y1_final",
  title_suffix = paste0("% Dying Cells (Final Day = ", unique(waterfall_data$total_days), ")"),
  output_suffix = "final"
)

###############################################################################
# Cumulative death (Evaluation day)
###############################################################################

plot_y1_boxplot(
  y_variable = "cumulative_death_evaluation",
  title_suffix = paste0("Cumulative Death (Evaluation Day = ", evaluation_day, ")"),
  output_suffix = "cumdeath_evaluation"
)


###############################################################################
# Cumulative death (Final time)
###############################################################################

plot_y1_boxplot(
  y_variable = "cumulative_death_final",
  title_suffix = paste0("Cumulative Death (Final Day = ", unique(waterfall_data$total_days), ")"),
  output_suffix = "cumdeath_final"
)


###############################################################################
# Survival AUC (Evaluation day)
###############################################################################

plot_y1_boxplot(
  y_variable = "survival_auc_evaluation",
  title_suffix = paste0("Survival AUC (Evaluation Day = ", evaluation_day, ")"),
  output_suffix = "survivalAUC_evaluation"
)


###############################################################################
# Survival AUC (Final time)
###############################################################################

plot_y1_boxplot(
  y_variable = "survival_auc_final",
  title_suffix = paste0("Survival AUC (Final Day = ", unique(waterfall_data$total_days), ")"),
  output_suffix = "survivalAUC_final"
)
