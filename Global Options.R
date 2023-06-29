#List of the global variables
glo_vars <- c("n","run_psa","psa_count", "disc_rate_costs", "disc_rate_QALYs", "disc_rate_LYs",
              "Population", "treatment", "Treatment_effect_duration_HBA","Treatment_effect_duration_BMI","Results_output", "KMmethod", "Number of cores")

##get the length of glo_vars, for checking below
length(glo_vars)

#create the global variables matrix to store the global varaibles
GlobalVars <- matrix(data=NA, nrow = length(glo_vars), ncol = 2)
#update the row names to be the list of global variables above
rownames(GlobalVars) <- glo_vars
colnames(GlobalVars) <- c("Value", "Description")

#drop the glo_vars varaible
rm(glo_vars)

GlobalVars["n", "Value"] <- 10000
GlobalVars["n", "Description"] <- "The number of patients to run through the model"

GlobalVars["run_psa", "Value"] <- F
GlobalVars["run_psa", "Description"] <- "T = run psa, F = run deterministic"

GlobalVars["psa_count", "Value"] <- 2
#Stop the model here if there is a discrepancy in the model mode and declared number 
#of PSA runs
if(is.numeric(GlobalVars["psa_count", "Value"])==T&GlobalVars["run_psa", "Value"]==F|
   is.na(GlobalVars["psa_count", "Value"])==T&GlobalVars["run_psa", "Value"]==T){
  stop("discrepancy between the psa mode and the number of psa runs")  
}
GlobalVars["psa_count", "Description"] <- "number of PSA runs. set to NA for deterministic runs"

GlobalVars["disc_rate_costs", "Value"] <- 0.035
GlobalVars["disc_rate_costs", "Description"] <- "discount rate for costs, writen as a numeric number"

GlobalVars["disc_rate_QALYs", "Value"] <- 0.035
GlobalVars["disc_rate_QALYs", "Description"] <- "discount rate for QALYs, writen as a numeric number"

GlobalVars["disc_rate_LYs", "Value"] <- 0.035
GlobalVars["disc_rate_LYs", "Description"] <- "discount rate for QALYs, writen as a numeric number"

GlobalVars["Population", "Value"] <- "NICE First"
GlobalVars["Population", "Description"] <- "text telling the model which population file to read in. Options are:
NICE First, NICE Second, NICE Third. Triple check this before running the model."

GlobalVars["treatment", "Value"] <- "GLOW"
GlobalVars["treatment", "Description"] <- "text term for the treatment. Cross reference that this is called correctly
in your code for creating an interventions effect matrix (should be in R/Interventions.R)"

GlobalVars["Treatment_effect_duration_HBA", "Value"] <-5
GlobalVars["Treatment_effect_duration_HBA", "Description"] <- "number of years of the treatment effect duration. Cross reference that this is called correctly
in your code for creating an interventions effect matrix (should be in R/Interventions.R)"

GlobalVars["Treatment_effect_duration_BMI", "Value"] <-10
GlobalVars["Treatment_effect_duration_BMI", "Description"] <- "number of years of the treatment effect duration. Cross reference that this is called correctly
in your code for creating an interventions effect matrix (should be in R/Interventions.R)"


GlobalVars["Results_output", "Value"] <- "Patient Level"
GlobalVars["Results_output", "Description"] <- "text term to produce different output types. Default is summary for PSA,
set to Patient Level to get the patient matrix. If not specfied, the determinisitc model will return a detailed year by year summary"

#Checks, Time_to_Event and Detailed modes are only run in a deterministic model
if(GlobalVars["run_psa", "Value"]==T&GlobalVars["Results_output", "Value"]=="Time_to_Event"|
   GlobalVars["run_psa", "Value"]==T&GlobalVars["Results_output", "Value"]=="Detailed"){
  stop("model outputs are not compatible with outputs from a PSA")  
}


GlobalVars["KMmethod", "Value"] <- "Events"
GlobalVars["KMmethod", "Description"] <- "Text term to indicate the method for calculating a Kaplan-Meier curve from
the simulation. Options are Events or Probabilities. Events uses the events in the model, probabilities uses the estimated probabilities.
Default is Events."

GlobalVars["Number of cores", "Value"] <- 100
GlobalVars["KMmethod", "Description"] <- "Number to declare the maximum number of cores to use when processing PSA results. Use this to decrease memory usage."
