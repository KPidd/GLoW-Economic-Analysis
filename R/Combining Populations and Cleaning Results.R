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

####'This page calls in all results from population subsets, combines and cleans output.
####'This page first defines necessary functions, then applies these functions for each analysis

Summarize_results<-function(DE,DEW){
  results<-matrix(NA,nrow=2,ncol=5)
  colnames(results)<-c("Total Discounted NHS Costs per person (£)","Total Discounted QALYs per person","Incremental Costs (£)","Incremental QALYs","Incremental expected net Monetary Benefit (£) – £20000 Threshold")
  rownames(results)<-c("DE","DEW")
  results["DE","Total Discounted NHS Costs per person (£)"]<-mean(DE[,"Discounted Costs"])
  results["DE","Total Discounted QALYs per person"]<-mean(DE[,"Discounted QALYs"])
  results["DEW","Total Discounted NHS Costs per person (£)"]<-mean(DEW[,"Discounted Costs"])
  results["DEW","Total Discounted QALYs per person"]<-mean(DEW[,"Discounted QALYs"])
  results["DEW","Incremental Costs (£)"]<-mean(DEW[,"Discounted Costs"]) - mean(DE[,"Discounted Costs"])
  results["DEW","Incremental QALYs"]<-mean(DEW[,"Discounted QALYs"]) - mean(DE[,"Discounted QALYs"])
  results["DEW","Incremental expected net Monetary Benefit (£) – £20000 Threshold"]<-(results["DEW","Incremental QALYs"]*lambda)-(results["DEW","Incremental Costs (£)"])
  return(results)
}
lambda<-20000
results_column_names <- c("ID","Life Years", 
                          "Discounted Life Years",
                          "QALYs", 
                          "Discounted QALYs", 
                          "Costs",
                          "Discounted Costs",
                          "1 year death", "1 year alive year end","1 year 1st MIs", "1 year 2nd MIs", "1 year 1st Strokes","1 year 2nd Strokes", "1 year CHF", "1 year IHD", "1 year blindness","1 year ulcers", "1 year 1st amputations",
                          "1 year 2nd amputations", "1 year renal failures", "1 year PVD", "1 year MMALB", "1 year Atrial Fibriliation", "1 year breast cancer","1 year colorectal cancer","1 year depression", "1 year oestoarthritis",
                          "2 year death", "2 year alive year end","2 year 1st MIs", "2 year 2nd MIs", "2 year 1st Strokes","2 year 2nd Strokes", "2 year CHF", "2 year IHD", "2 year blindness","2 year ulcers", "2 year 1st amputations",
                          "2 year 2nd amputations", "2 year renal failures", "2 year PVD", "2 year MMALB", "2 year Atrial Fibriliation", "2 year breast cancer","2 year colorectal cancer","2 year depression", "2 year oestoarthritis",
                          "3 year death", "3 year alive year end","3 year 1st MIs", "3 year 2nd MIs", "3 year 1st Strokes","3 year 2nd Strokes", "3 year CHF", "3 year IHD", "3 year blindness","3 year ulcers", "3 year 1st amputations",
                          "3 year 2nd amputations", "3 year renal failures", "3 year PVD", "3 year MMALB", "3 year Atrial Fibriliation", "3 year breast cancer","3 year colorectal cancer","3 year depression", "3 year oestoarthritis",
                          "4 year death", "4 year alive year end","4 year 1st MIs", "4 year 2nd MIs", "4 year 1st Strokes","4 year 2nd Strokes", "4 year CHF", "4 year IHD", "4 year blindness","4 year ulcers", "4 year 1st amputations",
                          "4 year 2nd amputations", "4 year renal failures", "4 year PVD", "4 year MMALB", "4 year Atrial Fibriliation", "4 year breast cancer","4 year colorectal cancer","4 year depression", "4 year oestoarthritis",
                          "5 year death", "5 year alive year end","5 year 1st MIs", "5 year 2nd MIs", "5 year 1st Strokes","5 year 2nd Strokes", "5 year CHF", "5 year IHD", "5 year blindness","5 year ulcers", "5 year 1st amputations",
                          "5 year 2nd amputations", "5 year renal failures", "5 year PVD", "5 year MMALB", "5 year Atrial Fibriliation", "5 year breast cancer","5 year colorectal cancer","5 year depression", "5 year oestoarthritis",
                          "6 year death", "6 year alive year end","6 year 1st MIs", "6 year 2nd MIs", "6 year 1st Strokes","6 year 2nd Strokes", "6 year CHF", "6 year IHD", "6 year blindness","6 year ulcers", "6 year 1st amputations",
                          "6 year 2nd amputations", "6 year renal failures", "6 year PVD", "6 year MMALB", "6 year Atrial Fibriliation", "6 year breast cancer","6 year colorectal cancer","6 year depression", "6 year oestoarthritis",
                          "7 year death", "7 year alive year end","7 year 1st MIs", "7 year 2nd MIs", "7 year 1st Strokes","7 year 2nd Strokes", "7 year CHF", "7 year IHD", "7 year blindness","7 year ulcers", "7 year 1st amputations",
                          "7 year 2nd amputations", "7 year renal failures", "7 year PVD", "7 year MMALB", "7 year Atrial Fibriliation", "7 year breast cancer","7 year colorectal cancer","7 year depression", "7 year oestoarthritis",
                          "8 year death", "8 year alive year end","8 year 1st MIs", "8 year 2nd MIs", "8 year 1st Strokes","8 year 2nd Strokes", "8 year CHF", "8 year IHD", "8 year blindness","8 year ulcers", "8 year 1st amputations",
                          "8 year 2nd amputations", "8 year renal failures", "8 year PVD", "8 year MMALB", "8 year Atrial Fibriliation", "8 year breast cancer","8 year colorectal cancer","8 year depression", "8 year oestoarthritis",
                          "9 year death", "9 year alive year end","9 year 1st MIs", "9 year 2nd MIs", "9 year 1st Strokes","9 year 2nd Strokes", "9 year CHF", "9 year IHD", "9 year blindness","9 year ulcers", "9 year 1st amputations",
                          "9 year 2nd amputations", "9 year renal failures", "9 year PVD", "9 year MMALB", "9 year Atrial Fibriliation", "9 year breast cancer","9 year colorectal cancer","9 year depression", "9 year oestoarthritis",
                          "10 year death", "10 year alive year end","10 year 1st MIs", "10 year 2nd MIs", "10 year 1st Strokes","10 year 2nd Strokes", "10 year CHF", "10 year IHD", "10 year blindness","10 year ulcers", "10 year 1st amputations",
                          "10 year 2nd amputations", "10 year renal failures", "10 year PVD", "10 year MMALB", "10 year Atrial Fibriliation", "10 year breast cancer","10 year colorectal cancer","10 year depression", "10 year oestoarthritis",
                          "11 year death", "11 year alive year end","11 year 1st MIs", "11 year 2nd MIs", "11 year 1st Strokes","11 year 2nd Strokes", "11 year CHF", "11 year IHD", "11 year blindness","11 year ulcers", "11 year 1st amputations",
                          "11 year 2nd amputations", "11 year renal failures", "11 year PVD", "11 year MMALB", "11 year Atrial Fibriliation", "11 year breast cancer","11 year colorectal cancer","11 year depression", "11 year oestoarthritis",
                          "12 year death", "12 year alive year end","12 year 1st MIs", "12 year 2nd MIs", "12 year 1st Strokes","12 year 2nd Strokes", "12 year CHF", "12 year IHD", "12 year blindness","12 year ulcers", "12 year 1st amputations",
                          "12 year 2nd amputations", "12 year renal failures", "12 year PVD", "12 year MMALB", "12 year Atrial Fibriliation", "12 year breast cancer","12 year colorectal cancer","12 year depression", "12 year oestoarthritis",
                          "13 year death", "13 year alive year end","13 year 1st MIs", "13 year 2nd MIs", "13 year 1st Strokes","13 year 2nd Strokes", "13 year CHF", "13 year IHD", "13 year blindness","13 year ulcers", "13 year 1st amputations",
                          "13 year 2nd amputations", "13 year renal failures", "13 year PVD", "13 year MMALB", "13 year Atrial Fibriliation", "13 year breast cancer","13 year colorectal cancer","13 year depression", "13 year oestoarthritis",
                          "14 year death", "14 year alive year end","14 year 1st MIs", "14 year 2nd MIs", "14 year 1st Strokes","14 year 2nd Strokes", "14 year CHF", "14 year IHD", "14 year blindness","14 year ulcers", "14 year 1st amputations",
                          "14 year 2nd amputations", "14 year renal failures", "14 year PVD", "14 year MMALB", "14 year Atrial Fibriliation", "14 year breast cancer","14 year colorectal cancer","14 year depression", "14 year oestoarthritis",
                          "15 year death", "15 year alive year end","15 year 1st MIs", "15 year 2nd MIs", "15 year 1st Strokes","15 year 2nd Strokes", "15 year CHF", "15 year IHD", "15 year blindness","15 year ulcers", "15 year 1st amputations",
                          "15 year 2nd amputations", "15 year renal failures", "15 year PVD", "15 year MMALB", "15 year Atrial Fibriliation", "15 year breast cancer","15 year colorectal cancer","15 year depression", "15 year oestoarthritis",
                          "16 year death", "16 year alive year end","16 year 1st MIs", "16 year 2nd MIs", "16 year 1st Strokes","16 year 2nd Strokes", "16 year CHF", "16 year IHD", "16 year blindness","16 year ulcers", "16 year 1st amputations",
                          "16 year 2nd amputations", "16 year renal failures", "16 year PVD", "16 year MMALB", "16 year Atrial Fibriliation", "16 year breast cancer","16 year colorectal cancer","16 year depression", "16 year oestoarthritis",
                          "17 year death", "17 year alive year end","17 year 1st MIs", "17 year 2nd MIs", "17 year 1st Strokes","17 year 2nd Strokes", "17 year CHF", "17 year IHD", "17 year blindness","17 year ulcers", "17 year 1st amputations",
                          "17 year 2nd amputations", "17 year renal failures", "17 year PVD", "17 year MMALB", "17 year Atrial Fibriliation", "17 year breast cancer","17 year colorectal cancer","17 year depression", "17 year oestoarthritis",
                          "18 year death", "18 year alive year end","18 year 1st MIs", "18 year 2nd MIs", "18 year 1st Strokes","18 year 2nd Strokes", "18 year CHF", "18 year IHD", "18 year blindness","18 year ulcers", "18 year 1st amputations",
                          "18 year 2nd amputations", "18 year renal failures", "18 year PVD", "18 year MMALB", "18 year Atrial Fibriliation", "18 year breast cancer","18 year colorectal cancer","18 year depression", "18 year oestoarthritis",
                          "19 year death", "19 year alive year end","19 year 1st MIs", "19 year 2nd MIs", "19 year 1st Strokes","19 year 2nd Strokes", "19 year CHF", "19 year IHD", "19 year blindness","19 year ulcers", "19 year 1st amputations",
                          "19 year 2nd amputations", "19 year renal failures", "19 year PVD", "19 year MMALB", "19 year Atrial Fibriliation", "19 year breast cancer","19 year colorectal cancer","19 year depression", "19 year oestoarthritis",
                          "20 year death", "20 year alive year end","20 year 1st MIs", "20 year 2nd MIs", "20 year 1st Strokes","20 year 2nd Strokes", "20 year CHF", "20 year IHD", "20 year blindness","20 year ulcers", "20 year 1st amputations",
                          "20 year 2nd amputations", "20 year renal failures", "20 year PVD", "20 year MMALB", "20 year Atrial Fibriliation", "20 year breast cancer","20 year colorectal cancer","20 year depression", "20 year oestoarthritis",
                          "21 year death", "21 year alive year end","21 year 1st MIs", "21 year 2nd MIs", "21 year 1st Strokes","21 year 2nd Strokes", "21 year CHF", "21 year IHD", "21 year blindness","21 year ulcers", "21 year 1st amputations",
                          "21 year 2nd amputations", "21 year renal failures", "21 year PVD", "21 year MMALB", "21 year Atrial Fibriliation", "21 year breast cancer","21 year colorectal cancer","21 year depression", "21 year oestoarthritis")


combining_populations<-function(results1,results2,results3,results4,results5,results6,results7,results8,results9,results10){
  colnames(results1)<-results_column_names
  colnames(results2)<-results_column_names
  colnames(results3)<-results_column_names
  colnames(results4)<-results_column_names
  colnames(results5)<-results_column_names
  colnames(results6)<-results_column_names
  colnames(results7)<-results_column_names
  colnames(results8)<-results_column_names
  colnames(results9)<-results_column_names
  colnames(results10)<-results_column_names
  results_all<-results1+results2+results3+results4+results5+results6+results7+results8+results9+results10
  results_all[,c("Life Years","Discounted Life Years","QALYs","Discounted QALYs","Costs","Discounted Costs")]<- results_all[,c("Life Years","Discounted Life Years","QALYs","Discounted QALYs","Costs","Discounted Costs")]/10
  return(results_all)
}

###########################################################################################################################
###########################        Table 2:Main incremental cost-effectiveness analysis      ################################

###########           assumes mixed F2F and online delivery
resultsDE_Mixed_F2F_A <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationA.csv")
resultsDEW_Mixed_F2F_A <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationA.csv")
Table2_MIXED_F2F_A<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationA.csv")
resultsDE_Mixed_F2F_B <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationB.csv")
resultsDEW_Mixed_F2F_B <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationB.csv")
Table2_MIXED_F2F_B<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationB.csv")
resultsDE_Mixed_F2F_C <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationC.csv")
resultsDEW_Mixed_F2F_C <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationC.csv")
Table2_MIXED_F2F_C<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationC.csv")
resultsDE_Mixed_F2F_D <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationD.csv")
resultsDEW_Mixed_F2F_D <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationD.csv")
Table2_MIXED_F2F_D<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationD.csv")
resultsDE_Mixed_F2F_E <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationE.csv")
resultsDEW_Mixed_F2F_E <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationE.csv")
Table2_MIXED_F2F_E<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationE.csv")
resultsDE_Mixed_F2F_F <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationF.csv")
resultsDEW_Mixed_F2F_F <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationF.csv")
Table2_MIXED_F2F_F<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationF.csv")
resultsDE_Mixed_F2F_G <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationG.csv")
resultsDEW_Mixed_F2F_G <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationG.csv")
Table2_MIXED_F2F_G<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationG.csv")
resultsDE_Mixed_F2F_H <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationH.csv")
resultsDEW_Mixed_F2F_H <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationH.csv")
Table2_MIXED_F2F_H<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationH.csv")
resultsDE_Mixed_F2F_I <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationI.csv")
resultsDEW_Mixed_F2F_I <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationI.csv")
Table2_MIXED_F2F_I<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationI.csv")
resultsDE_Mixed_F2F_J <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationJ.csv")
resultsDEW_Mixed_F2F_J <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationJ.csv")
Table2_MIXED_F2F_J<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationJ.csv")


resultsDE_Mixed_F2F<-combining_populations(resultsDE_Mixed_F2F_A,resultsDE_Mixed_F2F_B,resultsDE_Mixed_F2F_C,resultsDE_Mixed_F2F_D,resultsDE_Mixed_F2F_E,resultsDE_Mixed_F2F_F,resultsDE_Mixed_F2F_G,resultsDE_Mixed_F2F_H,resultsDE_Mixed_F2F_I,resultsDE_Mixed_F2F_J)
resultsDEW_Mixed_F2F<-combining_populations(resultsDEW_Mixed_F2F_A,resultsDEW_Mixed_F2F_B,resultsDEW_Mixed_F2F_C,resultsDEW_Mixed_F2F_D,resultsDEW_Mixed_F2F_E,resultsDEW_Mixed_F2F_F,resultsDEW_Mixed_F2F_G,resultsDEW_Mixed_F2F_H,resultsDEW_Mixed_F2F_I,resultsDEW_Mixed_F2F_J)
Table2_MIXED_F2F<-Summarize_results(resultsDE_Mixed_F2F,resultsDEW_Mixed_F2F)

##################       assumes F2F delivery ONLY
resultsDE_F2F_A <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationA.csv")
resultsDEW_F2F_A <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationA.csv")
Table2_F2F_A<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationA.csv")
resultsDE_F2F_B <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationB.csv")
resultsDEW_F2F_B <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationB.csv")
Table2_F2F_B<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationB.csv")
resultsDE_F2F_C <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationC.csv")
resultsDEW_F2F_C <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationC.csv")
Table2_F2F_C<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationC.csv")
resultsDE_F2F_D <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationD.csv")
resultsDEW_F2F_D <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationD.csv")
Table2_F2F_D<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationD.csv")
resultsDE_F2F_E <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationE.csv")
resultsDEW_F2F_E <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationE.csv")
Table2_F2F_E<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationE.csv")
resultsDE_F2F_F <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationF.csv")
resultsDEW_F2F_F <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationF.csv")
Table2_F2F_F<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationF.csv")
resultsDE_F2F_G <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationG.csv")
resultsDEW_F2F_G <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationG.csv")
Table2_F2F_G<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationG.csv")
resultsDE_F2F_H <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationH.csv")
resultsDEW_F2F_H <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationH.csv")
Table2_F2F_H<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationH.csv")
resultsDE_F2F_I <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationI.csv")
resultsDEW_F2F_I <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationI.csv")
Table2_F2F_I<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationI.csv")
resultsDE_F2F_J <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationJ.csv")
resultsDEW_F2F_J <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationJ.csv")
Table2_F2F_J<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationJ.csv")


resultsDE_F2F<-combining_populations(resultsDE_F2F_A,resultsDE_F2F_B,resultsDE_F2F_C,resultsDE_F2F_D,resultsDE_F2F_E,resultsDE_F2F_F,resultsDE_F2F_G,resultsDE_F2F_H,resultsDE_F2F_I,resultsDE_F2F_J)
resultsDEW_F2F<-combining_populations(resultsDEW_F2F_A,resultsDEW_F2F_B,resultsDEW_F2F_C,resultsDEW_F2F_D,resultsDEW_F2F_E,resultsDEW_F2F_F,resultsDEW_F2F_G,resultsDEW_F2F_H,resultsDEW_F2F_I,resultsDEW_F2F_J)
Table2_F2F<-Summarize_results(resultsDE_F2F,resultsDEW_F2F)

##################       assumes ONLINE delivery ONLY
resultsDE_ONLINE_A <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationA.csv")
resultsDEW_ONLINE_A <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationA.csv")
Table2_ONLINE_A<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationA.csv")
resultsDE_ONLINE_B <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationB.csv")
resultsDEW_ONLINE_B <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationB.csv")
Table2_ONLINE_B<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationB.csv")
resultsDE_ONLINE_C <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationC.csv")
resultsDEW_ONLINE_C <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationC.csv")
Table2_ONLINE_C<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationC.csv")
resultsDE_ONLINE_D <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationD.csv")
resultsDEW_ONLINE_D <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationD.csv")
Table2_ONLINE_D<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationD.csv")
resultsDE_ONLINE_E <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationE.csv")
resultsDEW_ONLINE_E <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationE.csv")
Table2_ONLINE_E<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationE.csv")
resultsDE_ONLINE_F <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationF.csv")
resultsDEW_ONLINE_F <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationF.csv")
Table2_ONLINE_F<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationF.csv")
resultsDE_ONLINE_G <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationG.csv")
resultsDEW_ONLINE_G <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationG.csv")
Table2_ONLINE_G<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationG.csv")
resultsDE_ONLINE_H <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationH.csv")
resultsDEW_ONLINE_H <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationH.csv")
Table2_ONLINE_H<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationH.csv")
resultsDE_ONLINE_I <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationI.csv")
resultsDEW_ONLINE_I <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationI.csv")
Table2_ONLINE_I<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationI.csv")
resultsDE_ONLINE_J <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationJ.csv")
resultsDEW_ONLINE_J <-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationJ.csv")
Table2_ONLINE_J<-read.csv("Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationJ.csv")


resultsDE_ONLINE<-combining_populations(resultsDE_ONLINE_A,resultsDE_ONLINE_B,resultsDE_ONLINE_C,resultsDE_ONLINE_D,resultsDE_ONLINE_E,resultsDE_ONLINE_F,resultsDE_ONLINE_G,resultsDE_ONLINE_H,resultsDE_ONLINE_I,resultsDE_ONLINE_J)
resultsDEW_ONLINE<-combining_populations(resultsDEW_ONLINE_A,resultsDEW_ONLINE_B,resultsDEW_ONLINE_C,resultsDEW_ONLINE_D,resultsDEW_ONLINE_E,resultsDEW_ONLINE_F,resultsDEW_ONLINE_G,resultsDEW_ONLINE_H,resultsDEW_ONLINE_I,resultsDEW_ONLINE_J)
Table2_ONLINE<-Summarize_results(resultsDE_ONLINE,resultsDEW_ONLINE)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
################################             Table 3: Sub-group analysis          ########################################

##################       Disease duration <=1 year
resultsDE_popDM1_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationA.csv")
resultsDEW_popDM1_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationA.csv")
Table3_popDM1_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationA.csv")
resultsDE_popDM1_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationB.csv")
resultsDEW_popDM1_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationB.csv")
Table3_popDM1_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationB.csv")
resultsDE_popDM1_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationC.csv")
resultsDEW_popDM1_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationC.csv")
Table3_popDM1_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationC.csv")
resultsDE_popDM1_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationD.csv")
resultsDEW_popDM1_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationD.csv")
Table3_popDM1_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationD.csv")
resultsDE_popDM1_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationE.csv")
resultsDEW_popDM1_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationE.csv")
Table3_popDM1_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationE.csv")
resultsDE_popDM1_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationF.csv")
resultsDEW_popDM1_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationF.csv")
Table3_popDM1_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationF.csv")
resultsDE_popDM1_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationG.csv")
resultsDEW_popDM1_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationG.csv")
Table3_popDM1_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationG.csv")
resultsDE_popDM1_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationH.csv")
resultsDEW_popDM1_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationH.csv")
Table3_popDM1_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationH.csv")
resultsDE_popDM1_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationI.csv")
resultsDEW_popDM1_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationI.csv")
Table3_popDM1_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationI.csv")
resultsDE_popDM1_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationJ.csv")
resultsDEW_popDM1_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationJ.csv")
Table3_popDM1_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM1_populationJ.csv")

resultsDE_popDM1<-combining_populations(resultsDE_popDM1_A,resultsDE_popDM1_B,resultsDE_popDM1_C,resultsDE_popDM1_D,resultsDE_popDM1_E,resultsDE_popDM1_F,resultsDE_popDM1_G,resultsDE_popDM1_H,resultsDE_popDM1_I,resultsDE_popDM1_J)
resultsDEW_popDM1<-combining_populations(resultsDEW_popDM1_A,resultsDEW_popDM1_B,resultsDEW_popDM1_C,resultsDEW_popDM1_D,resultsDEW_popDM1_E,resultsDEW_popDM1_F,resultsDEW_popDM1_G,resultsDEW_popDM1_H,resultsDEW_popDM1_I,resultsDEW_popDM1_J)
Table3_popDM1<-Summarize_results(resultsDE_popDM1,resultsDEW_popDM1)


##################       Disease duration 2-3 years
resultsDE_popDM2_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationA.csv")
resultsDEW_popDM2_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationA.csv")
Table3_popDM2_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationA.csv")
resultsDE_popDM2_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationB.csv")
resultsDEW_popDM2_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationB.csv")
Table3_popDM2_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationB.csv")
resultsDE_popDM2_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationC.csv")
resultsDEW_popDM2_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationC.csv")
Table3_popDM2_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationC.csv")
resultsDE_popDM2_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationD.csv")
resultsDEW_popDM2_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationD.csv")
Table3_popDM2_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationD.csv")
resultsDE_popDM2_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationE.csv")
resultsDEW_popDM2_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationE.csv")
Table3_popDM2_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationE.csv")
resultsDE_popDM2_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationF.csv")
resultsDEW_popDM2_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationF.csv")
Table3_popDM2_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationF.csv")
resultsDE_popDM2_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationG.csv")
resultsDEW_popDM2_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationG.csv")
Table3_popDM2_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationG.csv")
resultsDE_popDM2_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationH.csv")
resultsDEW_popDM2_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationH.csv")
Table3_popDM2_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationH.csv")
resultsDE_popDM2_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationI.csv")
resultsDEW_popDM2_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationI.csv")
Table3_popDM2_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationI.csv")
resultsDE_popDM2_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationJ.csv")
resultsDEW_popDM2_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationJ.csv")
Table3_popDM2_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popDM2_populationJ.csv")

resultsDE_popDM2<-combining_populations(resultsDE_popDM2_A,resultsDE_popDM2_B,resultsDE_popDM2_C,resultsDE_popDM2_D,resultsDE_popDM2_E,resultsDE_popDM2_F,resultsDE_popDM2_G,resultsDE_popDM2_H,resultsDE_popDM2_I,resultsDE_popDM2_J)
resultsDEW_popDM2<-combining_populations(resultsDEW_popDM2_A,resultsDEW_popDM2_B,resultsDEW_popDM2_C,resultsDEW_popDM2_D,resultsDEW_popDM2_E,resultsDEW_popDM2_F,resultsDEW_popDM2_G,resultsDEW_popDM2_H,resultsDEW_popDM2_I,resultsDEW_popDM2_J)
Table3_popDM2<-Summarize_results(resultsDE_popDM2,resultsDEW_popDM2)

##################       BMI groups: 28-30kg/m2

resultsDE_popBMI2830_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationA.csv")
resultsDEW_popBMI2830_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationA.csv")
Table3_popBMI2830_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationA.csv")
resultsDE_popBMI2830_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationB.csv")
resultsDEW_popBMI2830_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationB.csv")
Table3_popBMI2830_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationB.csv")
resultsDE_popBMI2830_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationC.csv")
resultsDEW_popBMI2830_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationC.csv")
Table3_popBMI2830_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationC.csv")
resultsDE_popBMI2830_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationD.csv")
resultsDEW_popBMI2830_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationD.csv")
Table3_popBMI2830_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationD.csv")
resultsDE_popBMI2830_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationE.csv")
resultsDEW_popBMI2830_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationE.csv")
Table3_popBMI2830_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationE.csv")
resultsDE_popBMI2830_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationF.csv")
resultsDEW_popBMI2830_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationF.csv")
Table3_popBMI2830_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationF.csv")
resultsDE_popBMI2830_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationG.csv")
resultsDEW_popBMI2830_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationG.csv")
Table3_popBMI2830_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationG.csv")
resultsDE_popBMI2830_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationH.csv")
resultsDEW_popBMI2830_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationH.csv")
Table3_popBMI2830_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationH.csv")
resultsDE_popBMI2830_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationI.csv")
resultsDEW_popBMI2830_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationI.csv")
Table3_popBMI2830_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationI.csv")
resultsDE_popBMI2830_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationJ.csv")
resultsDEW_popBMI2830_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationJ.csv")
Table3_popBMI2830_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationJ.csv")

resultsDE_popBMI2830<-combining_populations(resultsDE_popBMI2830_A,resultsDE_popBMI2830_B,resultsDE_popBMI2830_C,resultsDE_popBMI2830_D,resultsDE_popBMI2830_E,resultsDE_popBMI2830_F,resultsDE_popBMI2830_G,resultsDE_popBMI2830_H,resultsDE_popBMI2830_I,resultsDE_popBMI2830_J)
resultsDEW_popBMI2830<-combining_populations(resultsDEW_popBMI2830_A,resultsDEW_popBMI2830_B,resultsDEW_popBMI2830_C,resultsDEW_popBMI2830_D,resultsDEW_popBMI2830_E,resultsDEW_popBMI2830_F,resultsDEW_popBMI2830_G,resultsDEW_popBMI2830_H,resultsDEW_popBMI2830_I,resultsDEW_popBMI2830_J)
Table3_popBMI2830<-Summarize_results(resultsDE_popBMI2830,resultsDEW_popBMI2830)

##################       BMI groups: 30-35kg/m2

resultsDE_popBMI3035_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationA.csv")
resultsDEW_popBMI3035_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationA.csv")
Table3_popBMI3035_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationA.csv")
resultsDE_popBMI3035_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationB.csv")
resultsDEW_popBMI3035_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationB.csv")
Table3_popBMI3035_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationB.csv")
resultsDE_popBMI3035_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationC.csv")
resultsDEW_popBMI3035_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationC.csv")
Table3_popBMI3035_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationC.csv")
resultsDE_popBMI3035_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationD.csv")
resultsDEW_popBMI3035_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationD.csv")
Table3_popBMI3035_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationD.csv")
resultsDE_popBMI3035_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationE.csv")
resultsDEW_popBMI3035_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationE.csv")
Table3_popBMI3035_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationE.csv")
resultsDE_popBMI3035_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationF.csv")
resultsDEW_popBMI3035_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationF.csv")
Table3_popBMI3035_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationF.csv")
resultsDE_popBMI3035_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationG.csv")
resultsDEW_popBMI3035_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationG.csv")
Table3_popBMI3035_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationG.csv")
resultsDE_popBMI3035_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationH.csv")
resultsDEW_popBMI3035_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationH.csv")
Table3_popBMI3035_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationH.csv")
resultsDE_popBMI3035_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationI.csv")
resultsDEW_popBMI3035_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationI.csv")
Table3_popBMI3035_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationI.csv")
resultsDE_popBMI3035_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationJ.csv")
resultsDEW_popBMI3035_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationJ.csv")
Table3_popBMI3035_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationJ.csv")

resultsDE_popBMI3035<-combining_populations(resultsDE_popBMI3035_A,resultsDE_popBMI3035_B,resultsDE_popBMI3035_C,resultsDE_popBMI3035_D,resultsDE_popBMI3035_E,resultsDE_popBMI3035_F,resultsDE_popBMI3035_G,resultsDE_popBMI3035_H,resultsDE_popBMI3035_I,resultsDE_popBMI3035_J)
resultsDEW_popBMI3035<-combining_populations(resultsDEW_popBMI3035_A,resultsDEW_popBMI3035_B,resultsDEW_popBMI3035_C,resultsDEW_popBMI3035_D,resultsDEW_popBMI3035_E,resultsDEW_popBMI3035_F,resultsDEW_popBMI3035_G,resultsDEW_popBMI3035_H,resultsDEW_popBMI3035_I,resultsDEW_popBMI3035_J)
Table3_popBMI3035<-Summarize_results(resultsDE_popBMI3035,resultsDEW_popBMI3035)


##################       BMI groups: 35-40kg/m2

resultsDE_popBMI3540_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationA.csv")
resultsDEW_popBMI3540_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationA.csv")
Table3_popBMI3540_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationA.csv")
resultsDE_popBMI3540_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationB.csv")
resultsDEW_popBMI3540_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationB.csv")
Table3_popBMI3540_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationB.csv")
resultsDE_popBMI3540_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationC.csv")
resultsDEW_popBMI3540_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationC.csv")
Table3_popBMI3540_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationC.csv")
resultsDE_popBMI3540_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationD.csv")
resultsDEW_popBMI3540_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationD.csv")
Table3_popBMI3540_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationD.csv")
resultsDE_popBMI3540_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationE.csv")
resultsDEW_popBMI3540_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationE.csv")
Table3_popBMI3540_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationE.csv")
resultsDE_popBMI3540_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationF.csv")
resultsDEW_popBMI3540_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationF.csv")
Table3_popBMI3540_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationF.csv")
resultsDE_popBMI3540_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationG.csv")
resultsDEW_popBMI3540_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationG.csv")
Table3_popBMI3540_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationG.csv")
resultsDE_popBMI3540_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationH.csv")
resultsDEW_popBMI3540_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationH.csv")
Table3_popBMI3540_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationH.csv")
resultsDE_popBMI3540_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationI.csv")
resultsDEW_popBMI3540_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationI.csv")
Table3_popBMI3540_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationI.csv")
resultsDE_popBMI3540_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationJ.csv")
resultsDEW_popBMI3540_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationJ.csv")
Table3_popBMI3540_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationJ.csv")

resultsDE_popBMI3540<-combining_populations(resultsDE_popBMI3540_A,resultsDE_popBMI3540_B,resultsDE_popBMI3540_C,resultsDE_popBMI3540_D,resultsDE_popBMI3540_E,resultsDE_popBMI3540_F,resultsDE_popBMI3540_G,resultsDE_popBMI3540_H,resultsDE_popBMI3540_I,resultsDE_popBMI3540_J)
resultsDEW_popBMI3540<-combining_populations(resultsDEW_popBMI3540_A,resultsDEW_popBMI3540_B,resultsDEW_popBMI3540_C,resultsDEW_popBMI3540_D,resultsDEW_popBMI3540_E,resultsDEW_popBMI3540_F,resultsDEW_popBMI3540_G,resultsDEW_popBMI3540_H,resultsDEW_popBMI3540_I,resultsDEW_popBMI3540_J)
Table3_popBMI3540<-Summarize_results(resultsDE_popBMI3540,resultsDEW_popBMI3540)
##################       BMI groups: 40+kg/m2

resultsDE_popBMI40_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationA.csv")
resultsDEW_popBMI40_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationA.csv")
Table3_popBMI40_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationA.csv")
resultsDE_popBMI40_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationB.csv")
resultsDEW_popBMI40_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationB.csv")
Table3_popBMI40_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationB.csv")
resultsDE_popBMI40_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationC.csv")
resultsDEW_popBMI40_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationC.csv")
Table3_popBMI40_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationC.csv")
resultsDE_popBMI40_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationD.csv")
resultsDEW_popBMI40_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationD.csv")
Table3_popBMI40_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationD.csv")
resultsDE_popBMI40_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationE.csv")
resultsDEW_popBMI40_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationE.csv")
Table3_popBMI40_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationE.csv")
resultsDE_popBMI40_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationF.csv")
resultsDEW_popBMI40_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationF.csv")
Table3_popBMI40_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationF.csv")
resultsDE_popBMI40_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationG.csv")
resultsDEW_popBMI40_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationG.csv")
Table3_popBMI40_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationG.csv")
resultsDE_popBMI40_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationH.csv")
resultsDEW_popBMI40_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationH.csv")
Table3_popBMI40_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationH.csv")
resultsDE_popBMI40_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationI.csv")
resultsDEW_popBMI40_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationI.csv")
Table3_popBMI40_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationI.csv")
resultsDE_popBMI40_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationJ.csv")
resultsDEW_popBMI40_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationJ.csv")
Table3_popBMI40_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationJ.csv")

resultsDE_popBMI40<-combining_populations(resultsDE_popBMI40_A,resultsDE_popBMI40_B,resultsDE_popBMI40_C,resultsDE_popBMI40_D,resultsDE_popBMI40_E,resultsDE_popBMI40_F,resultsDE_popBMI40_G,resultsDE_popBMI40_H,resultsDE_popBMI40_I,resultsDE_popBMI40_J)
resultsDEW_popBMI40<-combining_populations(resultsDEW_popBMI40_A,resultsDEW_popBMI40_B,resultsDEW_popBMI40_C,resultsDEW_popBMI40_D,resultsDEW_popBMI40_E,resultsDEW_popBMI40_F,resultsDEW_popBMI40_G,resultsDEW_popBMI40_H,resultsDEW_popBMI40_I,resultsDEW_popBMI40_J)
Table3_popBMI40<-Summarize_results(resultsDE_popBMI40,resultsDEW_popBMI40)
##################      IMD Socioeconomic quantiles (1: least deprived)

resultsDE_popIMD1_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationA.csv")
resultsDEW_popIMD1_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationA.csv")
Table3_popIMD1_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationA.csv")
resultsDE_popIMD1_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationB.csv")
resultsDEW_popIMD1_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationB.csv")
Table3_popIMD1_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationB.csv")
resultsDE_popIMD1_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationC.csv")
resultsDEW_popIMD1_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationC.csv")
Table3_popIMD1_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationC.csv")
resultsDE_popIMD1_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationD.csv")
resultsDEW_popIMD1_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationD.csv")
Table3_popIMD1_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationD.csv")
resultsDE_popIMD1_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationE.csv")
resultsDEW_popIMD1_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationE.csv")
Table3_popIMD1_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationE.csv")
resultsDE_popIMD1_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationF.csv")
resultsDEW_popIMD1_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationF.csv")
Table3_popIMD1_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationF.csv")
resultsDE_popIMD1_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationG.csv")
resultsDEW_popIMD1_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationG.csv")
Table3_popIMD1_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationG.csv")
resultsDE_popIMD1_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationH.csv")
resultsDEW_popIMD1_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationH.csv")
Table3_popIMD1_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationH.csv")
resultsDE_popIMD1_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationI.csv")
resultsDEW_popIMD1_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationI.csv")
Table3_popIMD1_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationI.csv")
resultsDE_popIMD1_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationJ.csv")
resultsDEW_popIMD1_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationJ.csv")
Table3_popIMD1_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationJ.csv")

resultsDE_popIMD1<-combining_populations(resultsDE_popIMD1_A,resultsDE_popIMD1_B,resultsDE_popIMD1_C,resultsDE_popIMD1_D,resultsDE_popIMD1_E,resultsDE_popIMD1_F,resultsDE_popIMD1_G,resultsDE_popIMD1_H,resultsDE_popIMD1_I,resultsDE_popIMD1_J)
resultsDEW_popIMD1<-combining_populations(resultsDEW_popIMD1_A,resultsDEW_popIMD1_B,resultsDEW_popIMD1_C,resultsDEW_popIMD1_D,resultsDEW_popIMD1_E,resultsDEW_popIMD1_F,resultsDEW_popIMD1_G,resultsDEW_popIMD1_H,resultsDEW_popIMD1_I,resultsDEW_popIMD1_J)
Table3_popIMD1<-Summarize_results(resultsDE_popIMD1,resultsDEW_popIMD1)
##################      IMD Socioeconomic quantiles 2

resultsDE_popIMD2_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationA.csv")
resultsDEW_popIMD2_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationA.csv")
Table3_popIMD2_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationA.csv")
resultsDE_popIMD2_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationB.csv")
resultsDEW_popIMD2_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationB.csv")
Table3_popIMD2_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationB.csv")
resultsDE_popIMD2_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationC.csv")
resultsDEW_popIMD2_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationC.csv")
Table3_popIMD2_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationC.csv")
resultsDE_popIMD2_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationD.csv")
resultsDEW_popIMD2_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationD.csv")
Table3_popIMD2_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationD.csv")
resultsDE_popIMD2_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationE.csv")
resultsDEW_popIMD2_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationE.csv")
Table3_popIMD2_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationE.csv")
resultsDE_popIMD2_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationF.csv")
resultsDEW_popIMD2_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationF.csv")
Table3_popIMD2_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationF.csv")
resultsDE_popIMD2_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationG.csv")
resultsDEW_popIMD2_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationG.csv")
Table3_popIMD2_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationG.csv")
resultsDE_popIMD2_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationH.csv")
resultsDEW_popIMD2_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationH.csv")
Table3_popIMD2_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationH.csv")
resultsDE_popIMD2_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationI.csv")
resultsDEW_popIMD2_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationI.csv")
Table3_popIMD2_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationI.csv")
resultsDE_popIMD2_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationJ.csv")
resultsDEW_popIMD2_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationJ.csv")
Table3_popIMD2_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationJ.csv")

resultsDE_popIMD2<-combining_populations(resultsDE_popIMD2_A,resultsDE_popIMD2_B,resultsDE_popIMD2_C,resultsDE_popIMD2_D,resultsDE_popIMD2_E,resultsDE_popIMD2_F,resultsDE_popIMD2_G,resultsDE_popIMD2_H,resultsDE_popIMD2_I,resultsDE_popIMD2_J)
resultsDEW_popIMD2<-combining_populations(resultsDEW_popIMD2_A,resultsDEW_popIMD2_B,resultsDEW_popIMD2_C,resultsDEW_popIMD2_D,resultsDEW_popIMD2_E,resultsDEW_popIMD2_F,resultsDEW_popIMD2_G,resultsDEW_popIMD2_H,resultsDEW_popIMD2_I,resultsDEW_popIMD2_J)
Table3_popIMD2<-Summarize_results(resultsDE_popIMD2,resultsDEW_popIMD2)
##################      IMD Socioeconomic quantiles 3

resultsDE_popIMD3_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationA.csv")
resultsDEW_popIMD3_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationA.csv")
Table3_popIMD3_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationA.csv")
resultsDE_popIMD3_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationB.csv")
resultsDEW_popIMD3_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationB.csv")
Table3_popIMD3_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationB.csv")
resultsDE_popIMD3_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationC.csv")
resultsDEW_popIMD3_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationC.csv")
Table3_popIMD3_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationC.csv")
resultsDE_popIMD3_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationD.csv")
resultsDEW_popIMD3_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationD.csv")
Table3_popIMD3_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationD.csv")
resultsDE_popIMD3_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationE.csv")
resultsDEW_popIMD3_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationE.csv")
Table3_popIMD3_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationE.csv")
resultsDE_popIMD3_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationF.csv")
resultsDEW_popIMD3_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationF.csv")
Table3_popIMD3_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationF.csv")
resultsDE_popIMD3_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationG.csv")
resultsDEW_popIMD3_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationG.csv")
Table3_popIMD3_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationG.csv")
resultsDE_popIMD3_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationH.csv")
resultsDEW_popIMD3_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationH.csv")
Table3_popIMD3_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationH.csv")
resultsDE_popIMD3_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationI.csv")
resultsDEW_popIMD3_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationI.csv")
Table3_popIMD3_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationI.csv")
resultsDE_popIMD3_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationJ.csv")
resultsDEW_popIMD3_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationJ.csv")
Table3_popIMD3_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationJ.csv")

resultsDE_popIMD3<-combining_populations(resultsDE_popIMD3_A,resultsDE_popIMD3_B,resultsDE_popIMD3_C,resultsDE_popIMD3_D,resultsDE_popIMD3_E,resultsDE_popIMD3_F,resultsDE_popIMD3_G,resultsDE_popIMD3_H,resultsDE_popIMD3_I,resultsDE_popIMD3_J)
resultsDEW_popIMD3<-combining_populations(resultsDEW_popIMD3_A,resultsDEW_popIMD3_B,resultsDEW_popIMD3_C,resultsDEW_popIMD3_D,resultsDEW_popIMD3_E,resultsDEW_popIMD3_F,resultsDEW_popIMD3_G,resultsDEW_popIMD3_H,resultsDEW_popIMD3_I,resultsDEW_popIMD3_J)
Table3_popIMD3<-Summarize_results(resultsDE_popIMD3,resultsDEW_popIMD3)
##################      IMD Socioeconomic quantiles 4

resultsDE_popIMD4_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationA.csv")
resultsDEW_popIMD4_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationA.csv")
Table3_popIMD4_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationA.csv")
resultsDE_popIMD4_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationB.csv")
resultsDEW_popIMD4_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationB.csv")
Table3_popIMD4_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationB.csv")
resultsDE_popIMD4_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationC.csv")
resultsDEW_popIMD4_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationC.csv")
Table3_popIMD4_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationC.csv")
resultsDE_popIMD4_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationD.csv")
resultsDEW_popIMD4_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationD.csv")
Table3_popIMD4_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationD.csv")
resultsDE_popIMD4_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationE.csv")
resultsDEW_popIMD4_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationE.csv")
Table3_popIMD4_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationE.csv")
resultsDE_popIMD4_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationF.csv")
resultsDEW_popIMD4_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationF.csv")
Table3_popIMD4_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationF.csv")
resultsDE_popIMD4_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationG.csv")
resultsDEW_popIMD4_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationG.csv")
Table3_popIMD4_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationG.csv")
resultsDE_popIMD4_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationH.csv")
resultsDEW_popIMD4_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationH.csv")
Table3_popIMD4_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationH.csv")
resultsDE_popIMD4_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationI.csv")
resultsDEW_popIMD4_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationI.csv")
Table3_popIMD4_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationI.csv")
resultsDE_popIMD4_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationJ.csv")
resultsDEW_popIMD4_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationJ.csv")
Table3_popIMD4_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationJ.csv")

resultsDE_popIMD4<-combining_populations(resultsDE_popIMD4_A,resultsDE_popIMD4_B,resultsDE_popIMD4_C,resultsDE_popIMD4_D,resultsDE_popIMD4_E,resultsDE_popIMD4_F,resultsDE_popIMD4_G,resultsDE_popIMD4_H,resultsDE_popIMD4_I,resultsDE_popIMD4_J)
resultsDEW_popIMD4<-combining_populations(resultsDEW_popIMD4_A,resultsDEW_popIMD4_B,resultsDEW_popIMD4_C,resultsDEW_popIMD4_D,resultsDEW_popIMD4_E,resultsDEW_popIMD4_F,resultsDEW_popIMD4_G,resultsDEW_popIMD4_H,resultsDEW_popIMD4_I,resultsDEW_popIMD4_J)
Table3_popIMD4<-Summarize_results(resultsDE_popIMD4,resultsDEW_popIMD4)
##################      IMD Socioeconomic quantiles (5: most deprived)

resultsDE_popIMD5_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationA.csv")
resultsDEW_popIMD5_A <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationA.csv")
Table3_popIMD5_A<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationA.csv")
resultsDE_popIMD5_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationB.csv")
resultsDEW_popIMD5_B <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationB.csv")
Table3_popIMD5_B<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationB.csv")
resultsDE_popIMD5_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationC.csv")
resultsDEW_popIMD5_C <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationC.csv")
Table3_popIMD5_C<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationC.csv")
resultsDE_popIMD5_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationD.csv")
resultsDEW_popIMD5_D <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationD.csv")
Table3_popIMD5_D<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationD.csv")
resultsDE_popIMD5_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationE.csv")
resultsDEW_popIMD5_E <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationE.csv")
Table3_popIMD5_E<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationE.csv")
resultsDE_popIMD5_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationF.csv")
resultsDEW_popIMD5_F <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationF.csv")
Table3_popIMD5_F<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationF.csv")
resultsDE_popIMD5_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationG.csv")
resultsDEW_popIMD5_G <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationG.csv")
Table3_popIMD5_G<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationG.csv")
resultsDE_popIMD5_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationH.csv")
resultsDEW_popIMD5_H <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationH.csv")
Table3_popIMD5_H<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationH.csv")
resultsDE_popIMD5_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationI.csv")
resultsDEW_popIMD5_I <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationI.csv")
Table3_popIMD5_I<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationI.csv")
resultsDE_popIMD5_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationJ.csv")
resultsDEW_popIMD5_J <-read.csv("Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationJ.csv")
Table3_popIMD5_J<-read.csv("Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationJ.csv")

resultsDE_popIMD5<-combining_populations(resultsDE_popIMD5_A,resultsDE_popIMD5_B,resultsDE_popIMD5_C,resultsDE_popIMD5_D,resultsDE_popIMD5_E,resultsDE_popIMD5_F,resultsDE_popIMD5_G,resultsDE_popIMD5_H,resultsDE_popIMD5_I,resultsDE_popIMD5_J)
resultsDEW_popIMD5<-combining_populations(resultsDEW_popIMD5_A,resultsDEW_popIMD5_B,resultsDEW_popIMD5_C,resultsDEW_popIMD5_D,resultsDEW_popIMD5_E,resultsDEW_popIMD5_F,resultsDEW_popIMD5_G,resultsDEW_popIMD5_H,resultsDEW_popIMD5_I,resultsDEW_popIMD5_J)
Table3_popIMD5<-Summarize_results(resultsDE_popIMD5,resultsDEW_popIMD5)




###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
################################             Table 4: Sensitivity analysis          ########################################

##################       Treatment Effect Maintenance 3 years

resultsDE_TreatDur3_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationA.csv")
resultsDEW_TreatDur3_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationA.csv")
Table4_TreatDur3_A<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationA.csv")
resultsDE_TreatDur3_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationB.csv")
resultsDEW_TreatDur3_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationB.csv")
Table4_TreatDur3_B<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationB.csv")
resultsDE_TreatDur3_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationC.csv")
resultsDEW_TreatDur3_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationC.csv")
Table4_TreatDur3_C<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationC.csv")
resultsDE_TreatDur3_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationD.csv")
resultsDEW_TreatDur3_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationD.csv")
Table4_TreatDur3_D<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationD.csv")
resultsDE_TreatDur3_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationE.csv")
resultsDEW_TreatDur3_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationE.csv")
Table4_TreatDur3_E<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationE.csv")
resultsDE_TreatDur3_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationF.csv")
resultsDEW_TreatDur3_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationF.csv")
Table4_TreatDur3_F<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationF.csv")
resultsDE_TreatDur3_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationG.csv")
resultsDEW_TreatDur3_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationG.csv")
Table4_TreatDur3_G<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationG.csv")
resultsDE_TreatDur3_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationH.csv")
resultsDEW_TreatDur3_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationH.csv")
Table4_TreatDur3_H<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationH.csv")
resultsDE_TreatDur3_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationI.csv")
resultsDEW_TreatDur3_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationI.csv")
Table4_TreatDur3_I<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationI.csv")
resultsDE_TreatDur3_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationJ.csv")
resultsDEW_TreatDur3_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationJ.csv")
Table4_TreatDur3_J<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationJ.csv")

resultsDE_TreatDur3<-combining_populations(resultsDE_TreatDur3_A,resultsDE_TreatDur3_B,resultsDE_TreatDur3_C,resultsDE_TreatDur3_D,resultsDE_TreatDur3_E,resultsDE_TreatDur3_F,resultsDE_TreatDur3_G,resultsDE_TreatDur3_H,resultsDE_TreatDur3_I,resultsDE_TreatDur3_J)
resultsDEW_TreatDur3<-combining_populations(resultsDEW_TreatDur3_A,resultsDEW_TreatDur3_B,resultsDEW_TreatDur3_C,resultsDEW_TreatDur3_D,resultsDEW_TreatDur3_E,resultsDEW_TreatDur3_F,resultsDEW_TreatDur3_G,resultsDEW_TreatDur3_H,resultsDEW_TreatDur3_I,resultsDEW_TreatDur3_J)
Table4_TreatDur3<-Summarize_results(resultsDE_TreatDur3,resultsDEW_TreatDur3)


##################       Treatment Effect Maintenance 13 years

resultsDE_TreatDur13_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationA.csv")
resultsDEW_TreatDur13_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationA.csv")
Table4_TreatDur13_A<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationA.csv")
resultsDE_TreatDur13_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationB.csv")
resultsDEW_TreatDur13_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationB.csv")
Table4_TreatDur13_B<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationB.csv")
resultsDE_TreatDur13_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationC.csv")
resultsDEW_TreatDur13_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationC.csv")
Table4_TreatDur13_C<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationC.csv")
resultsDE_TreatDur13_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationD.csv")
resultsDEW_TreatDur13_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationD.csv")
Table4_TreatDur13_D<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationD.csv")
resultsDE_TreatDur13_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationE.csv")
resultsDEW_TreatDur13_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationE.csv")
Table4_TreatDur13_E<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationE.csv")
resultsDE_TreatDur13_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationF.csv")
resultsDEW_TreatDur13_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationF.csv")
Table4_TreatDur13_F<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationF.csv")
resultsDE_TreatDur13_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationG.csv")
resultsDEW_TreatDur13_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationG.csv")
Table4_TreatDur13_G<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationG.csv")
resultsDE_TreatDur13_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationH.csv")
resultsDEW_TreatDur13_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationH.csv")
Table4_TreatDur13_H<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationH.csv")
resultsDE_TreatDur13_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationI.csv")
resultsDEW_TreatDur13_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationI.csv")
Table4_TreatDur13_I<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationI.csv")
resultsDE_TreatDur13_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationJ.csv")
resultsDEW_TreatDur13_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationJ.csv")
Table4_TreatDur13_J<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationJ.csv")

resultsDE_TreatDur13<-combining_populations(resultsDE_TreatDur13_A,resultsDE_TreatDur13_B,resultsDE_TreatDur13_C,resultsDE_TreatDur13_D,resultsDE_TreatDur13_E,resultsDE_TreatDur13_F,resultsDE_TreatDur13_G,resultsDE_TreatDur13_H,resultsDE_TreatDur13_I,resultsDE_TreatDur13_J)
resultsDEW_TreatDur13<-combining_populations(resultsDEW_TreatDur13_A,resultsDEW_TreatDur13_B,resultsDEW_TreatDur13_C,resultsDEW_TreatDur13_D,resultsDEW_TreatDur13_E,resultsDEW_TreatDur13_F,resultsDEW_TreatDur13_G,resultsDEW_TreatDur13_H,resultsDEW_TreatDur13_I,resultsDEW_TreatDur13_J)
Table4_TreatDur13<-Summarize_results(resultsDE_TreatDur13,resultsDEW_TreatDur13)


##################       Cost of DESMOND face to face assuming optimistic utilisation (£95.05: £159 F2F and £7.47 online)

resultsDE_OptimisticUtil_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationA.csv")
resultsDEW_OptimisticUtil_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationA.csv")
Table4_OptimisticUtil_A<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationA.csv")
resultsDE_OptimisticUtil_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationB.csv")
resultsDEW_OptimisticUtil_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationB.csv")
Table4_OptimisticUtil_B<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationB.csv")
resultsDE_OptimisticUtil_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationC.csv")
resultsDEW_OptimisticUtil_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationC.csv")
Table4_OptimisticUtil_C<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationC.csv")
resultsDE_OptimisticUtil_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationD.csv")
resultsDEW_OptimisticUtil_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationD.csv")
Table4_OptimisticUtil_D<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationD.csv")
resultsDE_OptimisticUtil_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationE.csv")
resultsDEW_OptimisticUtil_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationE.csv")
Table4_OptimisticUtil_E<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationE.csv")
resultsDE_OptimisticUtil_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationF.csv")
resultsDEW_OptimisticUtil_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationF.csv")
Table4_OptimisticUtil_F<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationF.csv")
resultsDE_OptimisticUtil_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationG.csv")
resultsDEW_OptimisticUtil_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationG.csv")
Table4_OptimisticUtil_G<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationG.csv")
resultsDE_OptimisticUtil_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationH.csv")
resultsDEW_OptimisticUtil_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationH.csv")
Table4_OptimisticUtil_H<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationH.csv")
resultsDE_OptimisticUtil_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationI.csv")
resultsDEW_OptimisticUtil_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationI.csv")
Table4_OptimisticUtil_I<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationI.csv")
resultsDE_OptimisticUtil_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationJ.csv")
resultsDEW_OptimisticUtil_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationJ.csv")
Table4_OptimisticUtil_J<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationJ.csv")

resultsDE_OptimisticUtil<-combining_populations(resultsDE_OptimisticUtil_A,resultsDE_OptimisticUtil_B,resultsDE_OptimisticUtil_C,resultsDE_OptimisticUtil_D,resultsDE_OptimisticUtil_E,resultsDE_OptimisticUtil_F,resultsDE_OptimisticUtil_G,resultsDE_OptimisticUtil_H,resultsDE_OptimisticUtil_I,resultsDE_OptimisticUtil_J)
resultsDEW_OptimisticUtil<-combining_populations(resultsDEW_OptimisticUtil_A,resultsDEW_OptimisticUtil_B,resultsDEW_OptimisticUtil_C,resultsDEW_OptimisticUtil_D,resultsDEW_OptimisticUtil_E,resultsDEW_OptimisticUtil_F,resultsDEW_OptimisticUtil_G,resultsDEW_OptimisticUtil_H,resultsDEW_OptimisticUtil_I,resultsDEW_OptimisticUtil_J)
Table4_OptimisticUtil<-Summarize_results(resultsDE_OptimisticUtil,resultsDEW_OptimisticUtil)


######removing BMi utility decrements

resultsDE_zerobmidec_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationA.csv")
resultsDEW_zerobmidec_A <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationA.csv")
Table4_zerobmidec_A<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationA.csv")
resultsDE_zerobmidec_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationB.csv")
resultsDEW_zerobmidec_B <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationB.csv")
Table4_zerobmidec_B<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationB.csv")
resultsDE_zerobmidec_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationC.csv")
resultsDEW_zerobmidec_C <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationC.csv")
Table4_zerobmidec_C<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationC.csv")
resultsDE_zerobmidec_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationD.csv")
resultsDEW_zerobmidec_D <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationD.csv")
Table4_zerobmidec_D<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationD.csv")
resultsDE_zerobmidec_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationE.csv")
resultsDEW_zerobmidec_E <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationE.csv")
Table4_zerobmidec_E<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationE.csv")
resultsDE_zerobmidec_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationF.csv")
resultsDEW_zerobmidec_F <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationF.csv")
Table4_zerobmidec_F<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationF.csv")
resultsDE_zerobmidec_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationG.csv")
resultsDEW_zerobmidec_G <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationG.csv")
Table4_zerobmidec_G<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationG.csv")
resultsDE_zerobmidec_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationH.csv")
resultsDEW_zerobmidec_H <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationH.csv")
Table4_zerobmidec_H<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationH.csv")
resultsDE_zerobmidec_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationI.csv")
resultsDEW_zerobmidec_I <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationI.csv")
Table4_zerobmidec_I<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationI.csv")
resultsDE_zerobmidec_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationJ.csv")
resultsDEW_zerobmidec_J <-read.csv("Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationJ.csv")
Table4_zerobmidec_J<-read.csv("Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationJ.csv")

resultsDE_zerobmidec<-combining_populations(resultsDE_zerobmidec_A,resultsDE_zerobmidec_B,resultsDE_zerobmidec_C,resultsDE_zerobmidec_D,resultsDE_zerobmidec_E,resultsDE_zerobmidec_F,resultsDE_zerobmidec_G,resultsDE_zerobmidec_H,resultsDE_zerobmidec_I,resultsDE_zerobmidec_J)
resultsDEW_zerobmidec<-combining_populations(resultsDEW_zerobmidec_A,resultsDEW_zerobmidec_B,resultsDEW_zerobmidec_C,resultsDEW_zerobmidec_D,resultsDEW_zerobmidec_E,resultsDEW_zerobmidec_F,resultsDEW_zerobmidec_G,resultsDEW_zerobmidec_H,resultsDEW_zerobmidec_I,resultsDEW_zerobmidec_J)
Table4_zerobmidec<-Summarize_results(resultsDE_zerobmidec,resultsDEW_zerobmidec)



######removing calibration

resultsDE_nocalibration_A <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationA.csv")
resultsDEW_nocalibration_A <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationA.csv")
Table4_nocalibration_A<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationA.csv")
resultsDE_nocalibration_B <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationB.csv")
resultsDEW_nocalibration_B <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationB.csv")
Table4_nocalibration_B<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationB.csv")
resultsDE_nocalibration_C <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationC.csv")
resultsDEW_nocalibration_C <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationC.csv")
Table4_nocalibration_C<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationC.csv")
resultsDE_nocalibration_D <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationD.csv")
resultsDEW_nocalibration_D <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationD.csv")
Table4_nocalibration_D<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationD.csv")
resultsDE_nocalibration_E <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationE.csv")
resultsDEW_nocalibration_E <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationE.csv")
Table4_nocalibration_E<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationE.csv")
resultsDE_nocalibration_F <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationF.csv")
resultsDEW_nocalibration_F <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationF.csv")
Table4_nocalibration_F<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationF.csv")
resultsDE_nocalibration_G <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationG.csv")
resultsDEW_nocalibration_G <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationG.csv")
Table4_nocalibration_G<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationG.csv")
resultsDE_nocalibration_H <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationH.csv")
resultsDEW_nocalibration_H <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationH.csv")
Table4_nocalibration_H<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationH.csv")
resultsDE_nocalibration_I <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationI.csv")
resultsDEW_nocalibration_I <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationI.csv")
Table4_nocalibration_I<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationI.csv")
resultsDE_nocalibration_J <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_populationJ.csv")
resultsDEW_nocalibration_J <-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_populationJ.csv")
Table4_nocalibration_J<-read.csv("Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_populationJ.csv")

resultsDE_nocalibration<-combining_populations(resultsDE_nocalibration_A,resultsDE_nocalibration_B,resultsDE_nocalibration_C,resultsDE_nocalibration_D,resultsDE_nocalibration_E,resultsDE_nocalibration_F,resultsDE_nocalibration_G,resultsDE_nocalibration_H,resultsDE_nocalibration_I,resultsDE_nocalibration_J)
resultsDEW_nocalibration<-combining_populations(resultsDEW_nocalibration_A,resultsDEW_nocalibration_B,resultsDEW_nocalibration_C,resultsDEW_nocalibration_D,resultsDEW_nocalibration_E,resultsDEW_nocalibration_F,resultsDEW_nocalibration_G,resultsDEW_nocalibration_H,resultsDEW_nocalibration_I,resultsDEW_nocalibration_J)
Table4_nocalibration<-Summarize_results(resultsDE_nocalibration,resultsDEW_nocalibration)




### Save all combined outputs
write.csv(resultsDE_Mixed_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_MIXED_F2F_ALL.csv")
write.csv(resultsDEW_Mixed_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_MIXED_F2F_ALL.csv")
write.csv(Table2_MIXED_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_ALL.csv")

write.csv(resultsDE_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_ALL.csv")
write.csv(resultsDEW_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_ALL.csv")
write.csv(Table2_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_ALL.csv")

write.csv(resultsDE_ONLINE,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_ALL.csv")
write.csv(resultsDEW_ONLINE,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_ALL.csv")
write.csv(Table2_ONLINE,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_ALL.csv")

write.csv(resultsDE_popDM1,"Results/Table 3 - Sub-group analysis/resultsDE_popDM1_ALL.csv")
write.csv(resultsDEW_popDM1,"Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_ALL.csv")
write.csv(Table3_popDM1,"Results/Table 3 - Sub-group analysis/Table3_popDM1_ALL.csv")

write.csv(resultsDE_popDM2,"Results/Table 3 - Sub-group analysis/resultsDE_popDM2_ALL.csv")
write.csv(resultsDEW_popDM2,"Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_ALL.csv")
write.csv(Table3_popDM2,"Results/Table 3 - Sub-group analysis/Table3_popDM2_ALL.csv")

write.csv(resultsDE_popBMI2830,"Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_ALL.csv")
write.csv(resultsDEW_popBMI2830,"Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_ALL.csv")
write.csv(Table3_popBMI2830,"Results/Table 3 - Sub-group analysis/Table3_popBMI2830_ALL.csv")

write.csv(resultsDE_popBMI3035,"Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_ALL.csv")
write.csv(resultsDEW_popBMI3035,"Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_ALL.csv")
write.csv(Table3_popBMI3035,"Results/Table 3 - Sub-group analysis/Table3_popBMI3035_ALL.csv")

write.csv(resultsDE_popBMI3540,"Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_ALL.csv")
write.csv(resultsDEW_popBMI3540,"Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_ALL.csv")
write.csv(Table3_popBMI3540,"Results/Table 3 - Sub-group analysis/Table3_popBMI3540_ALL.csv")

write.csv(resultsDE_popBMI40,"Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_ALL.csv")
write.csv(resultsDEW_popBMI40,"Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_ALL.csv")
write.csv(Table3_popBMI40,"Results/Table 3 - Sub-group analysis/Table3_popBMI40_ALL.csv")

write.csv(resultsDE_popIMD1,"Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_ALL.csv")
write.csv(resultsDEW_popIMD1,"Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_ALL.csv")
write.csv(Table3_popIMD1,"Results/Table 3 - Sub-group analysis/Table3_popIMD1_ALL.csv")

write.csv(resultsDE_popIMD2,"Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_ALL.csv")
write.csv(resultsDEW_popIMD2,"Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_ALL.csv")
write.csv(Table3_popIMD2,"Results/Table 3 - Sub-group analysis/Table3_popIMD2_ALL.csv")

write.csv(resultsDE_popIMD3,"Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_ALL.csv")
write.csv(resultsDEW_popIMD3,"Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_ALL.csv")
write.csv(Table3_popIMD3,"Results/Table 3 - Sub-group analysis/Table3_popIMD3_ALL.csv")

write.csv(resultsDE_popIMD4,"Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_ALL.csv")
write.csv(resultsDEW_popIMD4,"Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_ALL.csv")
write.csv(Table3_popIMD4,"Results/Table 3 - Sub-group analysis/Table3_popIMD4_ALL.csv")

write.csv(resultsDE_popIMD5,"Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_ALL.csv")
write.csv(resultsDEW_popIMD5,"Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_ALL.csv")
write.csv(Table3_popIMD5,"Results/Table 3 - Sub-group analysis/Table3_popIMD5_ALL.csv")

write.csv(resultsDE_TreatDur3,"Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_ALL.csv")
write.csv(resultsDEW_TreatDur3,"Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_ALL.csv")
write.csv(Table4_TreatDur3,"Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_ALL.csv")

write.csv(resultsDE_TreatDur13,"Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_ALL.csv")
write.csv(resultsDEW_TreatDur13,"Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_ALL.csv")
write.csv(Table4_TreatDur13,"Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_ALL.csv")

write.csv(resultsDE_OptimisticUtil,"Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_ALL.csv")
write.csv(resultsDEW_OptimisticUtil,"Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_ALL.csv")
write.csv(Table4_OptimisticUtil,"Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_ALL.csv")

write.csv(resultsDE_zerobmidec,"Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_ALL.csv")
write.csv(resultsDEW_zerobmidec,"Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_ALL.csv")
write.csv(Table4_zerobmidec,"Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_ALL.csv")

write.csv(resultsDE_nocalibration,"Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDE_nocalibration_ALL.csv")
write.csv(resultsDEW_nocalibration,"Results/Table 4 - Sensitivity Analysis/No Calibration/resultsDEW_nocalibration_ALL.csv")
write.csv(Table4_nocalibration,"Results/Table 4 - Sensitivity Analysis/No Calibration/Table4_nocalibration_ALL.csv")


### Generate table of results
Table3<-matrix(NA,22,6)
colnames(Table3)<-c("Arm",colnames(Table3_popBMI2830))
Table3[,1]<-c("DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW","DE","DEW")
rownames(Table3)<-c("DM1_DE","DM1_DEW","DM2_DE","DM2_DEW","BMI2830_DE","BMI2830_DEW","BMI3035_DE","BMI3035_DEW","BMI3540_DE","BMI3540_DEW","BMI40_DE","BMI40_DEW","IMD1_DE","IMD1_DEW","IMD2_DE","IMD2_DEW","IMD3_DE","IMD3_DEW","IMD4_DE","IMD4_DEW","IMD5_DE","IMD5_DEW")
Table3[c("DM1_DE","DM1_DEW"),2:6]<-Table3_popDM1
Table3[c("DM2_DE","DM2_DEW"),2:6]<-Table3_popDM2
Table3[c("BMI2830_DE","BMI2830_DEW"),2:6]<-Table3_popBMI2830
Table3[c("BMI3035_DE","BMI3035_DEW"),2:6]<-Table3_popBMI3035
Table3[c("BMI3540_DE","BMI3540_DEW"),2:6]<-Table3_popBMI3540
Table3[c("BMI40_DE","BMI40_DEW"),2:6]<-Table3_popBMI40
Table3[c("IMD1_DE","IMD1_DEW"),2:6]<-Table3_popIMD1
Table3[c("IMD2_DE","IMD2_DEW"),2:6]<-Table3_popIMD2
Table3[c("IMD3_DE","IMD3_DEW"),2:6]<-Table3_popIMD3
Table3[c("IMD4_DE","IMD4_DEW"),2:6]<-Table3_popIMD4
Table3[c("IMD5_DE","IMD5_DEW"),2:6]<-Table3_popIMD5

Table3<-replace(Table3,is.na(Table3),"-")
write.csv(Table3,"Results/Table 3 - Sub-group analysis/Table3.csv")


Table4<-matrix(NA,6,6)
colnames(Table4)<-c("Arm",colnames(Table4_TreatDur13))
Table4[,1]<-c("DE","DEW","DE","DEW","DE","DEW")
rownames(Table4)<-c("TreatDur3_DE","TreatDur3_DEW","TreatDur13_DE","TreatDur13_DEW","zerobmidec_DE","zerobmidec_DEW")
Table4[c("TreatDur3_DE","TreatDur3_DEW"),2:6]<-Table4_TreatDur3
Table4[c("TreatDur13_DE","TreatDur13_DEW"),2:6]<-Table4_TreatDur13
Table4[c("zerobmidec_DE","zerobmidec_DEW"),2:6]<-Table4_zerobmidec
Table4<-replace(Table4,is.na(Table4),"-")
write.csv(Table4,"Results/Table 4 - Sensitivity Analysis/Table4.csv")


