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

