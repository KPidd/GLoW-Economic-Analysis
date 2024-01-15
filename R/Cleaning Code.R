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
###Results Cleaning
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

