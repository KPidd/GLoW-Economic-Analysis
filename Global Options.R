#    Sheffield Type 2 Diabetes Treatment Model version 4.0: with the GLoW Trial health economic analysis implemented.
#    Copyright (C) 2023 Pidd, Pollard, Breeze, Bates, Thomas, Mueller, Ahern, Griffin, Brennan

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#    Contact person: Katharine Pidd, Email: k.pidd@sheffield.ac.uk, 
#    Address: Regent Court, 30 Regent Court, Sheffield, United Kingdom, S1 4DA
#List of the global variables
glo_vars <- c("n","Population Subset","run_psa","psa_count", "disc_rate_costs", "disc_rate_QALYs", "disc_rate_LYs",
              "Population", "Treatment_effect_duration_HBA","Treatment_effect_duration_BMI","Results_output", "KMmethod", "Number of cores")

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

GlobalVars["Population Subset", "Value"] <- "A"
GlobalVars["Population Subset", "Description"] <- "Which population subset to run through the model (A=1:20000, B=20001:40000, C=40001:60000, D=60001:80000, E=80001:100000)"

GlobalVars["run_psa", "Value"] <- F
GlobalVars["run_psa", "Description"] <- "T = run psa, F = run deterministic"

GlobalVars["psa_count", "Value"] <- 2000
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

GlobalVars["Treatment_effect_duration_HBA", "Value"] <-5
GlobalVars["Treatment_effect_duration_HBA", "Description"] <- "number of years of the treatment effect duration. Cross reference that this is called correctly
in your code for creating an interventions effect matrix (should be in R/Interventions.R)"

GlobalVars["Treatment_effect_duration_BMI", "Value"] <-10
GlobalVars["Treatment_effect_duration_BMI", "Description"] <- "number of years of the treatment effect duration. Cross reference that this is called correctly
in your code for creating an interventions effect matrix (should be in R/Interventions.R)"


GlobalVars["Results_output", "Value"] <- "Summary"
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

GlobalVars["Number of cores", "Value"] <- 4
GlobalVars["KMmethod", "Description"] <- "Number to declare the maximum number of cores to use when processing PSA results. Use this to decrease memory usage."
