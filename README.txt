To run this model you must install R and R studio. Our analyses were run in R version 4.3.0

Before running any models you should add a folder called "Results" to this directory with three subfolders called "Table 2 - Main Incremental Cost-effectiveness Analysis", "Table 3 - Sub-group analysis", "Table 4 - Sensitivity Analysis". 
This folder is ingnored by git version control and contains the actual model results. 

To simply run the analyses detailed in the attached report, open the GLoW Economics Analysis.Rproj in R studio. Load and then run the .R script you want to run. 

The main file to run is the Main Model Run file.R file. This replicates the analyses undergone in the GLoW trial economic analysis (Reference to be added).

The other files contained in the directory are: 
all_model_files.R, this file calls in all R functions used in the model (all of the scripts containing these functions are located in the folder called R)
Global Options.R, this file generates a matrix containing all Global Options in the model (e.g. discount rate) and a column explaining what each Global Option is. This is pulled through the simulation to control what functions return

The other folders and files contained within are:
data - this folder contains the data used to populate the simulations
LifeTables.csv, is a reformatting of ONS lifetables for use in the simulation
parameter.rda, is a file containing all of the model parameters, with thier mean values in the first row and sampled uncertainty in the subsequent rows
PopulationVariables.csv, is a file containing a list of all vairables that each patient in the simulation will have and a description of what each varaible is

Populations - this contains the populations that are run through the model, prior to cleaning. We store these data in .csv files.
csv files contained here contained simulated populations, if modifying with data that is directly bootstrapped from individual patient level data, make sure that you do not load your population csv file into this folder and make it
available under seperate licence terms. 


R - this folder contains all of the user written functions in the simulation model
build_population.R - this contains the functions that convert the population csv files into the populations run throught the simulation
Cancer Risks.R - this contains the functions that estimate the probability of developing cancer in each year of the simulation
Costs.R - this contains the function that calculates the cost per year in the simulation
Depression.R - this contains the function that calculates the probability of developing depression in a given year
Generate Results Template.R - this contains the function that generates the correct matrix to store the model results top be outputted by the simulation
Generate Results.R - this contains functions that populate specific resulsts templates (althought  the defaul method is contained in Run model.R as this is necessary for most results collection methods)
generate_random.R - this contains the function that generates the common random numbers. This is different for each patient for each event for up to 100 years. These are the same across PSA runs and different analysis arms
intervention.R - this contains the functions that estimate the matrices given the intervention effects on A1c, BMI, SBP, HDL and LDL cholesterol for each patient in each year.  
LifeTableMortality.R - this contains the function that matches each patient's gender and age to all cause mortality from ONS life tables
Oestoarthritis functions.R - this contains the function that estimates each patient's probability of developing oestoarthritis
QALYs.R - this contains the function that estimates each patients QALY in accrued in each of the simulation. For patients who are alive this is equivalent to their utility in each year, but patients who die are assumed to die halfway
through the year. 
Run model.R - this contains the function that will conduct a full model run for a single row of the parameters matrix
Run simulation.R - this contains the function that runs the model for either the first row of the parameters matrix (deterministic) or across multiple rows to run the PSA analysis
UKPDS 82 risk functions.R - this contains the functions that predict the probability of all events in the UKPDS 82 risk functions
UKPDS 90 risk functions.R - this contains the functions that implement all of the equations in UKPDS 90 
Update Events.R - this contains the fuction that happens at the end of each simulation year. It ensures that all events that year are recorded as historical events and events are set to not happening for the next simulation year
Update Pat Chars.R - this contains the function that happens at the end of each simulation year. This ensures that all patient characteristics are set to what they should be for the next year in the simulation. 
Cleaning Code.R - this file contains the code that summarizes the outputs from each arm for each model run.
Combining Populations and Cleaning Results.R - this file contains the code to combine the population subsets and the results. It needs to be ran after all models are run through for all 10 population subsets. 

Changes made from original model:
Populations are sampled using trial data means and covariances.
Parameters have been updated to include updated and additional costs, the new intervention effect, GP utilisation equation, new relative risks, updated utility decrement.
R files:
build_population.R – Update to new available sampled population data. Set baseline LDL and HDL to ensure consistent trajectories.
Cancer Risks.R –U.
Costs.R – Updated to include costs of intervention and diabetes remission and to estimate the GP utilisation and related costs.
Depression.R – Unchanged.
Generate Results Template.R – Updated to add new results into output.
Generate Results.R - Updated to add new results into output.
generate_random.R – Added random numbers for GP utilisation function.
intervention.R – This has been updated to add intervention effect from the trial analysis or using the beta regression.
LifeTableMortality.R – Unchanged.
Oestoarthritis functions.R – Unchanged.
QALYs.R – Unchanged.
Run model.R – Updated to record different outputs.
Run simulation.R – Updated to include new intervention effect and decay.
UKPDS 82 risk functions.R – Unchanged.
UKPDS 90 risk functions.R – Unchanged.
Update Events.R – Added relative risks of MI and stroke and calibration.
Update Pat Chars.R – Updated to estimate individuals going into diabetes remission. 
Cleaning Code.R – Added to summarise outputs.
Combining Populations and Cleaning Results.R – Added to combine populations ran separately.

