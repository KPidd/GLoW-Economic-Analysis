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
####'This function runs the model for a given set of patients and parameters
####'@param population_, is the population matrix
####'@param parameters_, is the full parameters matrix
####'@param endtime_, is the number of years to run the simulation for
####'@param GlobalVars_, is the matrix giving the global variables
####'@param random_numbs_, is an array of common random numbers giving a random number
####'@param LifeTables_, is a dataframe containing life table information in formate
####'that can easily be matched to the population matrix
####'@param SOUR_, is the current second order uncertainty run
####'draw for each patient in each year for every event in which a random number 
####'is required 
####'@return results, is the results matrix
####'@return psaresults, is a summary set of results to produce in the PSA

run_simulation <- function(population_, parameters_, endtime_, treatment_, GlobalVars_, random_numbs_,LifeTables_, SOUR_){
  
  ##reduce the parameters matrix down to the correct row
  if(GlobalVars_["run_psa","Value"]==F){
    parameters_ <- parameters_[1,]
  }else{
    parameters_ <- parameters_[SOUR_+1,]
  }
  ##Create the underlying UKPDS trajectory matrices
  HBA1c_underlying  <- UKPDS_90_contrisk_A1c(population_,parameters_,endtime_)
  BMI_underlying    <- UKPDS_90_contrisk_BMI(population_,parameters_,endtime_)
  SBP_underlying    <- UKPDS_90_contrisk_SBP(population_,parameters_,endtime_)
  HDL_underlying    <- UKPDS_90_contrisk_HDL(population_,parameters_,endtime_)
  LDL_underlying    <- UKPDS_90_contrisk_LDL(population_,parameters_,endtime_)
  HEARTR_underlying <- UKPDS_90_HEARTR(population_,parameters_,endtime_)
  WBC_underlying    <- UKPDS_90_WBC(population_,parameters_,endtime_)
  HAEM_underlying   <- UKPDS_90_HAEM(population_,parameters_,endtime_)
  
  #Place to add Intervention effects
  HBA1c_INTV        <- initialise_intervention_dt_HbA1c(length(population_[,"ID"]),treatment_,parameters_,population_,endtime_, HBA1c_underlying, random_numbs_)
  BMI_INTV          <- initialise_intervention_dt_BMI(length(population_[,"ID"]), population_,treatment_,parameters_,endtime_)
  SBP_INTV          <- initialise_intervention_dt_SBP(length(population_[,"ID"]),treatment_,parameters_,endtime_)
  HDL_INTV          <- initialise_intervention_dt_HDL(length(population_[,"ID"]),treatment_,parameters_,endtime_)
  LDL_INTV          <- initialise_intervention_dt_LDL(length(population_[,"ID"]),treatment_,parameters_,endtime_)
  
  #start year at 0
  year <- 0
  
  #update population to have 12 month outcomes.
  population_[,"HBA"]<-population_[,"HBA"]+HBA1c_INTV[,2]
  population_[,"BMI"]<-population_[,"BMI"]+BMI_INTV[,2]
  
  #initialise the results matrix
  results <- GenerateResultsMatrix(GlobalVars_, endtime_)
  #initialise a PSA results matrix, if it is a PSA
  if(GlobalVars_["run_psa","Value"]==T){
  psaresults <- matrix(data=NA, nrow = as.numeric(GlobalVars_["psa_count", "Value"]), ncol = 5)  
  colnames(psaresults) <- c("Life Years per Patient", "Undiscounted QALYs per Patient", 
                            "QALYs per Patient", "Undiscounted Costs per Patient",
                            "Discounted Costs per Patient")
  }
  
  #run the model up to the specified end time or so long as at least two people are alive
  while (year < endtime_ & 
         (sum(is.na(population_[,"F_ALLCAUSE"]))) >= 1){
  #Get a logical vector indicating if people are alive, use the fact that FOTH
  #in the population matrix is NA if alive
  alive <- is.na(population_[,"F_ALLCAUSE"])
  #create matrix of ture = dead too for unit tests
  dead <- is.na(population_[,"F_ALLCAUSE"])==F
  
  #Estimate Diabetes Related complications and all cause deaths for this year
  population_ <- update_events_UKPDS82(population_,parameters_, treatment_, year, alive, random_numbs_, LifeTables_)
  #Estimate PVD,ATFIB,MMALB
  population_ <- update_events_UKPDS90(population_,parameters_, year, alive,random_numbs_)
  #Estimate Depression
  population_ <- update_events_SPHR_depression(population_,parameters_,year,alive,random_numbs_)
  #Estimate Osteoarthritis
  population_ <- update_events_SPHR_osteoarthritis(population_,parameters_,year,alive,random_numbs_)
  #Estimate Cancer incidence
  population_ <- update_events_SPHR_cancer(population_, parameters_,year, alive,random_numbs_)
  
  #Stop the model if dead people have events
  if(sum(population_[,"MI_E"][dead])  !=0|
     sum(population_[,"MI2_E"][dead]) != 0|
     sum(population_[,"AMP_E"][dead]) != 0| 
     sum(population_[,"AMP2_E"][dead]) != 0|
     sum(population_[,"BLIND_E"][dead]) != 0|
     sum(population_[,"ULCER_E"][dead]) != 0|
     sum(population_[,"RENAL_E"][dead]) != 0|
     sum(population_[,"CHF_E"][dead]) != 0|
     sum(population_[,"IHD_E"][dead]) != 0|
     sum(population_[,"STRO_E"][dead]) != 0|
     sum(population_[,"STRO2_E"][dead]) != 0|
     sum(population_[,"ATFIB_E"][dead]) != 0|
     sum(population_[,"PVD_E"][dead]) != 0|
     sum(population_[,"MMALB_E"][dead]) != 0|
     sum(population_[,"CANB_E"][dead]) != 0|
     sum(population_[,"CANC_E"][dead]) != 0|
     sum(population_[,"DEP_E"][dead]) != 0|
     sum(population_[,"OST_E"][dead]) != 0){
    #stop the model if the number of people with dead people are recorded as
    #having an event
    #push everything to the global enviroment for bug checking
    for (variable in ls()) {
      assign(variable, get(variable), envir = .GlobalEnv)
    }
    
    stop("dead people are getting events")
  }
  
  ##QALYs
  population_ <- calculate_QALYs(population_, parameters_,  year, alive, GlobalVars_)
  ##Costs
  population_ <- calculate_costs(population_, parameters_, year, alive, GlobalVars_, random_numbs_, treatment_)
  
  #Record results
  results <- GenerateDetailedresults(results,population_, year, alive, GlobalVars_)
  
  
  #update histories
  population_ <- update_history(population_,
                                HBA1c_underlying,
                                BMI_underlying,
                                SBP_underlying,    
                                HDL_underlying,
                                LDL_underlying, 
                                HEARTR_underlying,
                                WBC_underlying,
                                HAEM_underlying,
                                HBA1c_INTV,
                                BMI_INTV,
                                SBP_INTV,
                                HDL_INTV,
                                LDL_INTV,
                                year)
  
  population_ <- update_patchars(population_, parameters_, alive)
  HBA1c_INTV  <- INTE_HBA1c_decay(GlobalVars_, treatment_,HBA1c_INTV,year, endtime_)
  BMI_INTV  <- INTE_BMI_decay(GlobalVars_, treatment_,BMI_INTV,year, endtime_)
  #Unit test
  #Stop the model if there are events still in the population matrix
 if(sum(population_[,"MI_E"])  !=0|
     sum(population_[,"MI2_E"]) != 0|
     sum(population_[,"AMP_E"]) != 0| 
     sum(population_[,"AMP2_E"]) != 0|
     sum(population_[,"BLIND_E"]) != 0|
     sum(population_[,"ULCER_E"]) != 0|
     sum(population_[,"RENAL_E"]) != 0|
     sum(population_[,"CHF_E"]) != 0|
     sum(population_[,"IHD_E"]) != 0|
     sum(population_[,"STRO_E"]) != 0|
     sum(population_[,"STRO2_E"]) != 0){
   
     #push everything to the global enviroment for bug checking
     for (variable in ls()) {
       assign(variable, get(variable), envir = .GlobalEnv)
     }
   #stop the model if the events are not reset to 0
    stop("events are not set to zero before reset")
  }
  
  year <- year + 1
  
  }

  
  #For now return the Detailed results table or the population matrix if the run is deterministic
  if(GlobalVars_["Results_output", "Value"] == "Summary"&
     GlobalVars_["run_psa", "Value"]==T){
    psaresults <- matrix(data=NA,nrow=1,ncol=426)
    #Life Years
    psaresults[,1] <- sum((results["Undiscounted life years accrued",])/length(population_[,"ID"]),na.rm=TRUE)
    psaresults[,2] <- sum((results["Discounted life years accrued",])/length(population_[,"ID"]),na.rm=TRUE)
    psaresults[,3] <- sum((results["Undiscounted QALYs",])/length(population_[,"ID"]),na.rm=TRUE)
    psaresults[,4] <- sum((results["Discounted QALYs",])/length(population_[,"ID"]),na.rm=TRUE)
    psaresults[,5] <- sum((results["Undiscounted Costs",])/length(population_[,"ID"]),na.rm=TRUE)
    psaresults[,6] <- sum((results["Discounted Costs",])/length(population_[,"ID"]),na.rm=TRUE)
    #event rates
    psaresults[,7] <- results["Death", 1]
    psaresults[,8] <- results["Alive_year_end", 1]
    psaresults[,9] <- results["1st MI Hist", 1]
    psaresults[,10] <- results["2nd MI Hist",1]
    psaresults[,11] <- results["1st Stroke Hist",1]
    psaresults[,12] <- results["2nd Stroke Hist",1]
    psaresults[,13] <- results["CHF Hist",1]
    psaresults[,14] <- results["IHD Hist",1]
    psaresults[,15] <- results["Blindness Hist",1]
    psaresults[,16] <- results["Ulcer Hist",1]
    psaresults[,17] <- results["1st Amputation Hist",1]
    psaresults[,18] <- results["2nd Amputation Hist",1]
    psaresults[,19] <- results["Renal Failure Hist",1]
    psaresults[,20] <- results["PVD Hist",1]
    psaresults[,21] <- results["MMALB Hist",1]
    psaresults[,22] <- results["ATFIB Hist",1]
    psaresults[,23] <- results["Breast Cancer Hist",1]
    psaresults[,24] <- results["Colorectal Cancer Hist",1]
    psaresults[,25] <- results["Depression Hist",1]
    psaresults[,26] <- results["Osteoarthritis Hist",1]
    psaresults[,27] <- results["Death", 2]
    psaresults[,28] <- results["Alive_year_end", 2]
    psaresults[,29] <- results["1st MI Hist", 2]
    psaresults[,30] <- results["2nd MI Hist",2]
    psaresults[,31] <- results["1st Stroke Hist",2]
    psaresults[,32] <- results["2nd Stroke Hist",2]
    psaresults[,33] <- results["CHF Hist",2]
    psaresults[,34] <- results["IHD Hist",2]
    psaresults[,35] <- results["Blindness Hist",2]
    psaresults[,36] <- results["Ulcer Hist",2]
    psaresults[,37] <- results["1st Amputation Hist",2]
    psaresults[,38] <- results["2nd Amputation Hist",2]
    psaresults[,39] <- results["Renal Failure Hist",2]
    psaresults[,40] <- results["PVD Hist",2]
    psaresults[,41] <- results["MMALB Hist",2]
    psaresults[,42] <- results["ATFIB Hist",2]
    psaresults[,43] <- results["Breast Cancer Hist",2]
    psaresults[,44] <- results["Colorectal Cancer Hist",2]
    psaresults[,45] <- results["Depression Hist",2]
    psaresults[,46] <- results["Osteoarthritis Hist",2]
    psaresults[,47] <- results["Death", 3]
    psaresults[,48] <- results["Alive_year_end", 3]
    psaresults[,49] <- results["1st MI Hist", 3]
    psaresults[,50] <- results["2nd MI Hist",3]
    psaresults[,51] <- results["1st Stroke Hist",3]
    psaresults[,52] <- results["2nd Stroke Hist",3]
    psaresults[,53] <- results["CHF Hist",3]
    psaresults[,54] <- results["IHD Hist",3]
    psaresults[,55] <- results["Blindness Hist",3]
    psaresults[,56] <- results["Ulcer Hist",3]
    psaresults[,57] <- results["1st Amputation Hist",3]
    psaresults[,58] <- results["2nd Amputation Hist",3]
    psaresults[,59] <- results["Renal Failure Hist",3]
    psaresults[,60] <- results["PVD Hist",3]
    psaresults[,61] <- results["MMALB Hist",3]
    psaresults[,62] <- results["ATFIB Hist",3]
    psaresults[,63] <- results["Breast Cancer Hist",3]
    psaresults[,64] <- results["Colorectal Cancer Hist",3]
    psaresults[,65] <- results["Depression Hist",3]
    psaresults[,66] <- results["Osteoarthritis Hist",3]
    psaresults[,67] <- results["Death", 4]
    psaresults[,68] <- results["Alive_year_end", 4]
    psaresults[,69] <- results["1st MI Hist", 4]
    psaresults[,70] <- results["2nd MI Hist",4]
    psaresults[,71] <- results["1st Stroke Hist",4]
    psaresults[,72] <- results["2nd Stroke Hist",4]
    psaresults[,73] <- results["CHF Hist",4]
    psaresults[,74] <- results["IHD Hist",4]
    psaresults[,75] <- results["Blindness Hist",4]
    psaresults[,76] <- results["Ulcer Hist",4]
    psaresults[,77] <- results["1st Amputation Hist",4]
    psaresults[,78] <- results["2nd Amputation Hist",4]
    psaresults[,79] <- results["Renal Failure Hist",4]
    psaresults[,80] <- results["PVD Hist",4]
    psaresults[,81] <- results["MMALB Hist",4]
    psaresults[,82] <- results["ATFIB Hist",4]
    psaresults[,83] <- results["Breast Cancer Hist",4]
    psaresults[,84] <- results["Colorectal Cancer Hist",4]
    psaresults[,85] <- results["Depression Hist",4]
    psaresults[,86] <- results["Osteoarthritis Hist",4]
    psaresults[,87] <- results["Death", 5]
    psaresults[,88] <- results["Alive_year_end", 5]
    psaresults[,89] <- results["1st MI Hist", 5]
    psaresults[,90] <- results["2nd MI Hist",5]
    psaresults[,91] <- results["1st Stroke Hist",5]
    psaresults[,92] <- results["2nd Stroke Hist",5]
    psaresults[,93] <- results["CHF Hist",5]
    psaresults[,94] <- results["IHD Hist",5]
    psaresults[,95] <- results["Blindness Hist",5]
    psaresults[,96] <- results["Ulcer Hist",5]
    psaresults[,97] <- results["1st Amputation Hist",5]
    psaresults[,98] <- results["2nd Amputation Hist",5]
    psaresults[,99] <- results["Renal Failure Hist",5]
    psaresults[,100] <- results["PVD Hist",5]
    psaresults[,101] <- results["MMALB Hist",5]
    psaresults[,102] <- results["ATFIB Hist",5]
    psaresults[,103] <- results["Breast Cancer Hist",5]
    psaresults[,104] <- results["Colorectal Cancer Hist",5]
    psaresults[,105] <- results["Depression Hist",5]
    psaresults[,106] <- results["Osteoarthritis Hist",5]
    psaresults[,107] <- results["Death", 6]
    psaresults[,108] <- results["Alive_year_end", 6]
    psaresults[,109] <- results["1st MI Hist", 6]
    psaresults[,110] <- results["2nd MI Hist",6]
    psaresults[,111] <- results["1st Stroke Hist",6]
    psaresults[,112] <- results["2nd Stroke Hist",6]
    psaresults[,113] <- results["CHF Hist",6]
    psaresults[,114] <- results["IHD Hist",6]
    psaresults[,115] <- results["Blindness Hist",6]
    psaresults[,116] <- results["Ulcer Hist",6]
    psaresults[,117] <- results["1st Amputation Hist",6]
    psaresults[,118] <- results["2nd Amputation Hist",6]
    psaresults[,119] <- results["Renal Failure Hist",6]
    psaresults[,120] <- results["PVD Hist",6]
    psaresults[,121] <- results["MMALB Hist",6]
    psaresults[,122] <- results["ATFIB Hist",6]
    psaresults[,123] <- results["Breast Cancer Hist",6]
    psaresults[,124] <- results["Colorectal Cancer Hist",6]
    psaresults[,125] <- results["Depression Hist",6]
    psaresults[,126] <- results["Osteoarthritis Hist",6]
    psaresults[,127] <- results["Death", 7]
    psaresults[,128] <- results["Alive_year_end", 7]
    psaresults[,129] <- results["1st MI Hist", 7]
    psaresults[,130] <- results["2nd MI Hist",7]
    psaresults[,131] <- results["1st Stroke Hist",7]
    psaresults[,132] <- results["2nd Stroke Hist",7]
    psaresults[,133] <- results["CHF Hist",7]
    psaresults[,134] <- results["IHD Hist",7]
    psaresults[,135] <- results["Blindness Hist",7]
    psaresults[,136] <- results["Ulcer Hist",7]
    psaresults[,137] <- results["1st Amputation Hist",7]
    psaresults[,138] <- results["2nd Amputation Hist",7]
    psaresults[,139] <- results["Renal Failure Hist",7]
    psaresults[,140] <- results["PVD Hist",7]
    psaresults[,141] <- results["MMALB Hist",7]
    psaresults[,142] <- results["ATFIB Hist",7]
    psaresults[,143] <- results["Breast Cancer Hist",7]
    psaresults[,144] <- results["Colorectal Cancer Hist",7]
    psaresults[,145] <- results["Depression Hist",7]
    psaresults[,146] <- results["Osteoarthritis Hist",7]
    psaresults[,147] <- results["Death", 8]
    psaresults[,148] <- results["Alive_year_end", 8]
    psaresults[,149] <- results["1st MI Hist", 8]
    psaresults[,150] <- results["2nd MI Hist",8]
    psaresults[,151] <- results["1st Stroke Hist",8]
    psaresults[,152] <- results["2nd Stroke Hist",8]
    psaresults[,153] <- results["CHF Hist",8]
    psaresults[,154] <- results["IHD Hist",8]
    psaresults[,155] <- results["Blindness Hist",8]
    psaresults[,156] <- results["Ulcer Hist",8]
    psaresults[,157] <- results["1st Amputation Hist",8]
    psaresults[,158] <- results["2nd Amputation Hist",8]
    psaresults[,159] <- results["Renal Failure Hist",8]
    psaresults[,160] <- results["PVD Hist",8]
    psaresults[,161] <- results["MMALB Hist",8]
    psaresults[,162] <- results["ATFIB Hist",8]
    psaresults[,163] <- results["Breast Cancer Hist",8]
    psaresults[,164] <- results["Colorectal Cancer Hist",8]
    psaresults[,165] <- results["Depression Hist",8]
    psaresults[,166] <- results["Osteoarthritis Hist",8]
    psaresults[,167] <- results["Death", 9]
    psaresults[,168] <- results["Alive_year_end", 9]
    psaresults[,169] <- results["1st MI Hist", 9]
    psaresults[,170] <- results["2nd MI Hist",9]
    psaresults[,171] <- results["1st Stroke Hist",9]
    psaresults[,172] <- results["2nd Stroke Hist",9]
    psaresults[,173] <- results["CHF Hist",9]
    psaresults[,174] <- results["IHD Hist",9]
    psaresults[,175] <- results["Blindness Hist",9]
    psaresults[,176] <- results["Ulcer Hist",9]
    psaresults[,177] <- results["1st Amputation Hist",9]
    psaresults[,178] <- results["2nd Amputation Hist",9]
    psaresults[,179] <- results["Renal Failure Hist",9]
    psaresults[,180] <- results["PVD Hist",9]
    psaresults[,181] <- results["MMALB Hist",9]
    psaresults[,182] <- results["ATFIB Hist",9]
    psaresults[,183] <- results["Breast Cancer Hist",9]
    psaresults[,184] <- results["Colorectal Cancer Hist",9]
    psaresults[,185] <- results["Depression Hist",9]
    psaresults[,186] <- results["Osteoarthritis Hist",9]
    psaresults[,187] <- results["Death", 10]
    psaresults[,188] <- results["Alive_year_end", 10]
    psaresults[,189] <- results["1st MI Hist", 10]
    psaresults[,190] <- results["2nd MI Hist",10]
    psaresults[,191] <- results["1st Stroke Hist",10]
    psaresults[,192] <- results["2nd Stroke Hist",10]
    psaresults[,193] <- results["CHF Hist",10]
    psaresults[,194] <- results["IHD Hist",10]
    psaresults[,195] <- results["Blindness Hist",10]
    psaresults[,196] <- results["Ulcer Hist",10]
    psaresults[,197] <- results["1st Amputation Hist",10]
    psaresults[,198] <- results["2nd Amputation Hist",10]
    psaresults[,199] <- results["Renal Failure Hist",10]
    psaresults[,200] <- results["PVD Hist",10]
    psaresults[,201] <- results["MMALB Hist",10]
    psaresults[,202] <- results["ATFIB Hist",10]
    psaresults[,203] <- results["Breast Cancer Hist",10]
    psaresults[,204] <- results["Colorectal Cancer Hist",10]
    psaresults[,205] <- results["Depression Hist",10]
    psaresults[,206] <- results["Osteoarthritis Hist",10]
    psaresults[,207] <- results["Death", 11]
    psaresults[,208] <- results["Alive_year_end", 11]
    psaresults[,209] <- results["1st MI Hist", 11]
    psaresults[,210] <- results["2nd MI Hist",11]
    psaresults[,211] <- results["1st Stroke Hist",11]
    psaresults[,212] <- results["2nd Stroke Hist",11]
    psaresults[,213] <- results["CHF Hist",11]
    psaresults[,214] <- results["IHD Hist",11]
    psaresults[,215] <- results["Blindness Hist",11]
    psaresults[,216] <- results["Ulcer Hist",11]
    psaresults[,217] <- results["1st Amputation Hist",11]
    psaresults[,218] <- results["2nd Amputation Hist",11]
    psaresults[,219] <- results["Renal Failure Hist",11]
    psaresults[,220] <- results["PVD Hist",11]
    psaresults[,221] <- results["MMALB Hist",11]
    psaresults[,222] <- results["ATFIB Hist",11]
    psaresults[,223] <- results["Breast Cancer Hist",11]
    psaresults[,224] <- results["Colorectal Cancer Hist",11]
    psaresults[,225] <- results["Depression Hist",11]
    psaresults[,226] <- results["Osteoarthritis Hist",11]
    psaresults[,227] <- results["Death", 12]
    psaresults[,228] <- results["Alive_year_end", 12]
    psaresults[,229] <- results["1st MI Hist", 12]
    psaresults[,230] <- results["2nd MI Hist",12]
    psaresults[,231] <- results["1st Stroke Hist",12]
    psaresults[,232] <- results["2nd Stroke Hist",12]
    psaresults[,233] <- results["CHF Hist",12]
    psaresults[,234] <- results["IHD Hist",12]
    psaresults[,235] <- results["Blindness Hist",12]
    psaresults[,236] <- results["Ulcer Hist",12]
    psaresults[,237] <- results["1st Amputation Hist",12]
    psaresults[,238] <- results["2nd Amputation Hist",12]
    psaresults[,239] <- results["Renal Failure Hist",12]
    psaresults[,240] <- results["PVD Hist",12]
    psaresults[,241] <- results["MMALB Hist",12]
    psaresults[,242] <- results["ATFIB Hist",12]
    psaresults[,243] <- results["Breast Cancer Hist",12]
    psaresults[,244] <- results["Colorectal Cancer Hist",12]
    psaresults[,245] <- results["Depression Hist",12]
    psaresults[,246] <- results["Osteoarthritis Hist",12]
    psaresults[,247] <- results["Death", 13]
    psaresults[,248] <- results["Alive_year_end", 13]
    psaresults[,249] <- results["1st MI Hist", 13]
    psaresults[,250] <- results["2nd MI Hist",13]
    psaresults[,251] <- results["1st Stroke Hist",13]
    psaresults[,252] <- results["2nd Stroke Hist",13]
    psaresults[,253] <- results["CHF Hist",13]
    psaresults[,254] <- results["IHD Hist",13]
    psaresults[,255] <- results["Blindness Hist",13]
    psaresults[,256] <- results["Ulcer Hist",13]
    psaresults[,257] <- results["1st Amputation Hist",13]
    psaresults[,258] <- results["2nd Amputation Hist",13]
    psaresults[,259] <- results["Renal Failure Hist",13]
    psaresults[,260] <- results["PVD Hist",13]
    psaresults[,261] <- results["MMALB Hist",13]
    psaresults[,262] <- results["ATFIB Hist",13]
    psaresults[,263] <- results["Breast Cancer Hist",13]
    psaresults[,264] <- results["Colorectal Cancer Hist",13]
    psaresults[,265] <- results["Depression Hist",13]
    psaresults[,266] <- results["Osteoarthritis Hist",13]
    psaresults[,267] <- results["Death", 14]
    psaresults[,268] <- results["Alive_year_end", 14]
    psaresults[,269] <- results["1st MI Hist", 14]
    psaresults[,270] <- results["2nd MI Hist",14]
    psaresults[,271] <- results["1st Stroke Hist",14]
    psaresults[,272] <- results["2nd Stroke Hist",14]
    psaresults[,273] <- results["CHF Hist",14]
    psaresults[,274] <- results["IHD Hist",14]
    psaresults[,275] <- results["Blindness Hist",14]
    psaresults[,276] <- results["Ulcer Hist",14]
    psaresults[,277] <- results["1st Amputation Hist",14]
    psaresults[,278] <- results["2nd Amputation Hist",14]
    psaresults[,279] <- results["Renal Failure Hist",14]
    psaresults[,280] <- results["PVD Hist",14]
    psaresults[,281] <- results["MMALB Hist",14]
    psaresults[,282] <- results["ATFIB Hist",14]
    psaresults[,283] <- results["Breast Cancer Hist",14]
    psaresults[,284] <- results["Colorectal Cancer Hist",14]
    psaresults[,285] <- results["Depression Hist",14]
    psaresults[,286] <- results["Osteoarthritis Hist",14]
    psaresults[,287] <- results["Death", 15]
    psaresults[,288] <- results["Alive_year_end", 15]
    psaresults[,289] <- results["1st MI Hist", 15]
    psaresults[,290] <- results["2nd MI Hist",15]
    psaresults[,291] <- results["1st Stroke Hist",15]
    psaresults[,292] <- results["2nd Stroke Hist",15]
    psaresults[,293] <- results["CHF Hist",15]
    psaresults[,294] <- results["IHD Hist",15]
    psaresults[,295] <- results["Blindness Hist",15]
    psaresults[,296] <- results["Ulcer Hist",15]
    psaresults[,297] <- results["1st Amputation Hist",15]
    psaresults[,298] <- results["2nd Amputation Hist",15]
    psaresults[,299] <- results["Renal Failure Hist",15]
    psaresults[,300] <- results["PVD Hist",15]
    psaresults[,301] <- results["MMALB Hist",15]
    psaresults[,302] <- results["ATFIB Hist",15]
    psaresults[,303] <- results["Breast Cancer Hist",15]
    psaresults[,304] <- results["Colorectal Cancer Hist",15]
    psaresults[,305] <- results["Depression Hist",15]
    psaresults[,306] <- results["Osteoarthritis Hist",15]
    psaresults[,307] <- results["Death", 16]
    psaresults[,308] <- results["Alive_year_end", 16]
    psaresults[,309] <- results["1st MI Hist", 16]
    psaresults[,310] <- results["2nd MI Hist",16]
    psaresults[,311] <- results["1st Stroke Hist",16]
    psaresults[,312] <- results["2nd Stroke Hist",16]
    psaresults[,313] <- results["CHF Hist",16]
    psaresults[,314] <- results["IHD Hist",16]
    psaresults[,315] <- results["Blindness Hist",16]
    psaresults[,316] <- results["Ulcer Hist",16]
    psaresults[,317] <- results["1st Amputation Hist",16]
    psaresults[,318] <- results["2nd Amputation Hist",16]
    psaresults[,319] <- results["Renal Failure Hist",16]
    psaresults[,320] <- results["PVD Hist",16]
    psaresults[,321] <- results["MMALB Hist",16]
    psaresults[,322] <- results["ATFIB Hist",16]
    psaresults[,323] <- results["Breast Cancer Hist",16]
    psaresults[,324] <- results["Colorectal Cancer Hist",16]
    psaresults[,325] <- results["Depression Hist",16]
    psaresults[,326] <- results["Osteoarthritis Hist",16]
    psaresults[,327] <- results["Death", 17]
    psaresults[,328] <- results["Alive_year_end", 17]
    psaresults[,329] <- results["1st MI Hist", 17]
    psaresults[,330] <- results["2nd MI Hist",17]
    psaresults[,331] <- results["1st Stroke Hist",17]
    psaresults[,332] <- results["2nd Stroke Hist",17]
    psaresults[,333] <- results["CHF Hist",17]
    psaresults[,334] <- results["IHD Hist",17]
    psaresults[,335] <- results["Blindness Hist",17]
    psaresults[,336] <- results["Ulcer Hist",17]
    psaresults[,337] <- results["1st Amputation Hist",17]
    psaresults[,338] <- results["2nd Amputation Hist",17]
    psaresults[,339] <- results["Renal Failure Hist",17]
    psaresults[,340] <- results["PVD Hist",17]
    psaresults[,341] <- results["MMALB Hist",17]
    psaresults[,342] <- results["ATFIB Hist",17]
    psaresults[,343] <- results["Breast Cancer Hist",17]
    psaresults[,344] <- results["Colorectal Cancer Hist",17]
    psaresults[,345] <- results["Depression Hist",17]
    psaresults[,346] <- results["Osteoarthritis Hist",17]
    psaresults[,347] <- results["Death", 18]
    psaresults[,348] <- results["Alive_year_end", 18]
    psaresults[,349] <- results["1st MI Hist", 18]
    psaresults[,350] <- results["2nd MI Hist",18]
    psaresults[,351] <- results["1st Stroke Hist",18]
    psaresults[,352] <- results["2nd Stroke Hist",18]
    psaresults[,353] <- results["CHF Hist",18]
    psaresults[,354] <- results["IHD Hist",18]
    psaresults[,355] <- results["Blindness Hist",18]
    psaresults[,356] <- results["Ulcer Hist",18]
    psaresults[,357] <- results["1st Amputation Hist",18]
    psaresults[,358] <- results["2nd Amputation Hist",18]
    psaresults[,359] <- results["Renal Failure Hist",18]
    psaresults[,360] <- results["PVD Hist",18]
    psaresults[,361] <- results["MMALB Hist",18]
    psaresults[,362] <- results["ATFIB Hist",18]
    psaresults[,363] <- results["Breast Cancer Hist",18]
    psaresults[,364] <- results["Colorectal Cancer Hist",18]
    psaresults[,365] <- results["Depression Hist",18]
    psaresults[,366] <- results["Osteoarthritis Hist",18]
    psaresults[,367] <- results["Death", 19]
    psaresults[,368] <- results["Alive_year_end", 19]
    psaresults[,369] <- results["1st MI Hist", 19]
    psaresults[,370] <- results["2nd MI Hist",19]
    psaresults[,371] <- results["1st Stroke Hist",19]
    psaresults[,372] <- results["2nd Stroke Hist",19]
    psaresults[,373] <- results["CHF Hist",19]
    psaresults[,374] <- results["IHD Hist",19]
    psaresults[,375] <- results["Blindness Hist",19]
    psaresults[,376] <- results["Ulcer Hist",19]
    psaresults[,377] <- results["1st Amputation Hist",19]
    psaresults[,378] <- results["2nd Amputation Hist",19]
    psaresults[,379] <- results["Renal Failure Hist",19]
    psaresults[,380] <- results["PVD Hist",19]
    psaresults[,381] <- results["MMALB Hist",19]
    psaresults[,382] <- results["ATFIB Hist",19]
    psaresults[,383] <- results["Breast Cancer Hist",19]
    psaresults[,384] <- results["Colorectal Cancer Hist",19]
    psaresults[,385] <- results["Depression Hist",19]
    psaresults[,386] <- results["Osteoarthritis Hist",19]
    psaresults[,387] <- results["Death", 20]
    psaresults[,388] <- results["Alive_year_end", 20]
    psaresults[,389] <- results["1st MI Hist", 20]
    psaresults[,390] <- results["2nd MI Hist",20]
    psaresults[,391] <- results["1st Stroke Hist",20]
    psaresults[,392] <- results["2nd Stroke Hist",20]
    psaresults[,393] <- results["CHF Hist",20]
    psaresults[,394] <- results["IHD Hist",20]
    psaresults[,395] <- results["Blindness Hist",20]
    psaresults[,396] <- results["Ulcer Hist",20]
    psaresults[,397] <- results["1st Amputation Hist",20]
    psaresults[,398] <- results["2nd Amputation Hist",20]
    psaresults[,399] <- results["Renal Failure Hist",20]
    psaresults[,400] <- results["PVD Hist",20]
    psaresults[,401] <- results["MMALB Hist",20]
    psaresults[,402] <- results["ATFIB Hist",20]
    psaresults[,403] <- results["Breast Cancer Hist",20]
    psaresults[,404] <- results["Colorectal Cancer Hist",20]
    psaresults[,405] <- results["Depression Hist",20]
    psaresults[,406] <- results["Osteoarthritis Hist",20]
    psaresults[,407] <- results["Death", 21]
    psaresults[,408] <- results["Alive_year_end", 21]
    psaresults[,409] <- results["1st MI Hist", 21]
    psaresults[,410] <- results["2nd MI Hist",21]
    psaresults[,411] <- results["1st Stroke Hist",21]
    psaresults[,412] <- results["2nd Stroke Hist",21]
    psaresults[,413] <- results["CHF Hist",21]
    psaresults[,414] <- results["IHD Hist",21]
    psaresults[,415] <- results["Blindness Hist",21]
    psaresults[,416] <- results["Ulcer Hist",21]
    psaresults[,417] <- results["1st Amputation Hist",21]
    psaresults[,418] <- results["2nd Amputation Hist",21]
    psaresults[,419] <- results["Renal Failure Hist",21]
    psaresults[,420] <- results["PVD Hist",21]
    psaresults[,421] <- results["MMALB Hist",21]
    psaresults[,422] <- results["ATFIB Hist",21]
    psaresults[,423] <- results["Breast Cancer Hist",21]
    psaresults[,424] <- results["Colorectal Cancer Hist",21]
    psaresults[,425] <- results["Depression Hist",21]
    psaresults[,426] <- results["Osteoarthritis Hist",21]
    
    #delete the original results matrix
    rm(results)
    return(psaresults)
  }else if(GlobalVars_["Results_output","Value"]=="Patient Level"&
     GlobalVars_["run_psa","Value"]==F){#option to produce the patient 
    #charateristics matrix for checking results stability for number of patients
    return(population_)
  }else{#default is the 
  return(results)
  }
}