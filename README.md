A mechanistic digital twin model for epigenetic therapy optimization in triple-negative breast cancer 

This repository contains the codes used to develop and simulate digital twin model for epigenetic therapy optimization in triple-negative breast cancer , as presented in the associated manuscript.
The digital twin model integrates mechanistic modeling of chromatin dynamics and tumor progression with pharmacological inputs to study treatment effects, resistance mechanisms, and alternative therapeutic strategies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Repository Structure
•	SB_Digital_Twin_Main_Model_Population.R
Main digital twin model for combined epigenetic and signaling therapy (tazemetostat + ipatasertib) for population optimization

•	SB_Digital_Twin_Main_Model_Patient_Specific.R
Main digital twin model for combined epigenetic and signaling therapy (tazemetostat + ipatasertib) for patient-specific optimization

•	SB_Digital_Twin_Resist_Model.R
Extension of the model to study resistance mechanisms 

•	SB_Digital_Twin_Capiv.R
Alternative treatment scenario using capivasertib

•	Posterior_samples_population.Rds
Pre-estimated parameters for the population study

•	Posterior_samples_patient_i.Rds
Pre-estimated parameters for the patient-specific study, with i = 1,2,3,4 for patient 1,2,3,4, respectively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup
Before running any script, you must update the output directories in each file:
output_location <- "/Users/..."           # to be modified
output_location_figures <- "/Users/..."  # to be modified
Set these paths to your desired local directories where results and figures will be saved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reproducibility
•	All simulations use estimated parameter values stored in:
Posterior_samples_population.Rds
Posterior_samples_patient_i.Rds
with i = 1,2,3,4 for patient 1,2,3,4, respectively.

