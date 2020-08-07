require(corHMM)

#get ages of transitions
make.trans.table<-function(corhmm.tree,corhmm.node.states,corhmm.tip.states){
	node.ages<-NULL
	node.df<-NULL
	#print(corhmm.object["phy"])
	node.ages<-branching.times(corhmm.tree)
	names(node.ages)<-(length(corhmm.tree$tip.label)+1):(nrow(corhmm.node.states) + length(corhmm.tree$tip.label))
	#print(node.ages)
	node.df<-as.data.frame(corhmm.tree$edge)
	#print(node.df)
	colnames(node.df)<-c("From","To")
	node.df[,c("FromAge","FromState","ToAge","ToState","MicroGain","MacroGain","NonGain","TrentGain","TrebGain","CyanoGain","CephalodiaGain","AnyCyanoGain","AnyCyanoLoss")]<-NA	
	for(x in 1:nrow(node.df)){
		node.df$FromAge[x]<-node.ages[names(node.ages) %in% node.df$From[x]][[1]]
		#subtracted by 1 here bc corhmm numbers them from 1-6 instead of 0-5
		node.df$FromState[x]<-corhmm.tree$node.label[node.df$From[x]-length(corhmm.tree$tip.label)]-1
		if(node.df$To[x]<=length(corhmm.tree$tip.label)){
			#print(node.df$To[x])
			node.df$ToAge[x]<-0
			#print(dx$DietModslashboth[dx$Taxon %in% ard.recon$phy$tip.label[node.df$To[x]]])
			node.df$ToState[x]<-corhmm.tip.states[corhmm.tip.states$Taxon %in% corhmm.tree$tip.label[node.df$To[x]],2]
		}
		if(node.df$To[x]>length(corhmm.tree$tip.label)){
			node.df$ToAge[x]<-node.ages[names(node.ages) %in% node.df$To[x]][[1]]
			#subtracted by 1 here bc corhmm numbers them from 1-6 instead of 0-5
			node.df$ToState[x]<-corhmm.tree$node.label[node.df$To[x]-length(corhmm.tree$tip.label)]-1
		}
		if(node.df$FromState[x] %in% c(0,1,2,4,6) && node.df$ToState[x] %in% c(0,1,2,4,6,9) || is.na(node.df$ToState[x])){
			node.df$MicroGain[x]<-FALSE
			node.df$MacroGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(3,5,7) && node.df$ToState[x] %in% c(3,5,7,8) || is.na(node.df$ToState[x])){
			node.df$MicroGain[x]<-FALSE
			node.df$MacroGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0,1,2,4,6) && !node.df$ToState[x] %in% c(0,1,2,4,6,9) && !is.na(node.df$ToState[x])){
			node.df$MicroGain[x]<-FALSE
			node.df$MacroGain[x]<-TRUE
		}
		if(node.df$FromState[x] %in% c(3,5,7) && !node.df$ToState[x] %in% c(3,5,7,8) && !is.na(node.df$ToState[x])){
			node.df$MicroGain[x]<-TRUE
			node.df$MacroGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% 0 && node.df$ToState[x] %in% 0 || is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% 1 && node.df$ToState[x] %in% 1 || is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(2,3,6,7) && node.df$ToState[x] %in% c(2,3,6,7) || is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(4,5) && node.df$ToState[x] %in% c(4,5) || is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(1:7) && node.df$ToState[x] %in% 0 && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-TRUE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0,2:9) && node.df$ToState[x] %in% 1 && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-TRUE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0,1,4,5) && node.df$ToState[x] %in% c(2,3,6,7) && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-TRUE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0,1,2,3,6,7) && node.df$ToState[x] %in% c(4,5) && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-TRUE
		}
		if(node.df$FromState[x] %in% 0 && node.df$ToState[x] %in% 9 && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-TRUE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(2,3,6,7) && node.df$ToState[x] %in% 9 && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-TRUE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(1,4,5) && node.df$ToState[x] %in% 9 && !is.na(node.df$ToState[x])){
			print("double-gain state145 to 9")
			node.df$NonGain[x]<-TRUE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-TRUE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(2,3,6,7) && node.df$ToState[x] %in% 8 && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-FALSE
			node.df$CyanoGain[x]<-TRUE
		}
		if(node.df$FromState[x] %in% c(4,5) && node.df$ToState[x] %in% 8 && !is.na(node.df$ToState[x])){
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-TRUE
			node.df$CyanoGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0,1) && node.df$ToState[x] %in% 8 && !is.na(node.df$ToState[x])){
			print("double-gain state01 to 8")
			node.df$NonGain[x]<-FALSE
			node.df$TrentGain[x]<-FALSE
			node.df$TrebGain[x]<-TRUE
			node.df$CyanoGain[x]<-TRUE
		}
		if(node.df$FromState[x] %in% c(0:5) && node.df$ToState[x] %in% c(0:5,9) || is.na(node.df$ToState[x])){
			node.df$CephalodiaGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(6:7) && node.df$ToState[x] %in% c(6:7) || is.na(node.df$ToState[x])){
			node.df$CephalodiaGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0:5) && node.df$ToState[x] %in% 8 || is.na(node.df$ToState[x])){
			node.df$CephalodiaGain[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(0:5,9) && !node.df$ToState[x] %in% c(0:5,9) || is.na(node.df$ToState[x])){
			node.df$CephalodiaGain[x]<-TRUE
		}
		if(node.df$FromState[x] %in% c(0:3,9) && node.df$ToState[x] %in% c(0:3,9) || is.na(node.df$ToState[x])){
			node.df$AnyCyanoGain[x]<-FALSE
			node.df$AnyCyanoLoss[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(4:8) && node.df$ToState[x] %in% c(4:8) || is.na(node.df$ToState[x])){
			node.df$AnyCyanoGain[x]<-FALSE
			node.df$AnyCyanoLoss[x]<-FALSE
		}
		if(node.df$FromState[x] %in% c(4:8) && !node.df$ToState[x] %in% c(4:8) && !is.na(node.df$ToState[x])){
			node.df$AnyCyanoGain[x]<-FALSE
			node.df$AnyCyanoLoss[x]<-TRUE
		}
		if(node.df$FromState[x] %in% c(0:3,9) && !node.df$ToState[x] %in% c(0:3,9) && !is.na(node.df$ToState[x])){
			node.df$AnyCyanoGain[x]<-TRUE
			node.df$AnyCyanoLoss[x]<-FALSE
		}
	}
return(node.df)
}

file.path<-"/home/mpnelsen/examl_bootstrapping/"
mods<-read.csv(file="/home/mpnelsen/examl_bootstrapping/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_MULTISTATE_cephalodia.anycyano.csv",stringsAsFactors=FALSE)
rownames(mods)<-mods$Taxon

for(q in 1:100){
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	tr<-my.rate.mat.corr.cleanerer.fitz.no.diagn$phy
	mods<-mods[match(tr$tip.label,rownames(mods)),]
	tr.tab.correlated<-	make.trans.table(corhmm.tree=my.rate.mat.corr.cleanerer.fitz.no.diagn$phy,corhmm.node.states=my.rate.mat.corr.cleanerer.fitz.no.diagn$states,corhmm.tip.states=mods)
	write.csv(tr.tab.correlated,file=paste(file.path,"bootstrap_",q,"_node.info.csv",sep=""),row.names=FALSE)
}




recons<-c("MicroGain","MacroGain","NonGain","TrentGain","TrebGain","CyanoGain","CephalodiaGain","AnyCyanoGain","AnyCyanoLoss")
path<-"/home/mpnelsen/examl_bootstrapping/"
cum.df<-data.frame(matrix(vector(),0,9,dimnames=list(c(),c("From","To","FromAge","FromState","ToAge","ToState","Change","Reconstruction","Bootstrap"))),stringsAsFactors=FALSE)
for(q in 1:100){
	for(x in 1:length(recons)){
		df<-read.csv(file=paste(path,"bootstrap_",q,"_node.info.csv",sep=""),stringsAsFactors=FALSE)
		if(recons[x] %in% "MicroGain"){
			df<-df[,c(1:7)]	
		}
		if(recons[x] %in% "MacroGain"){
			df<-df[,c(1:6,8)]	
		}
		if(recons[x] %in% "NonGain"){
			df<-df[,c(1:6,9)]	
		}
		if(recons[x] %in% "TrentGain"){
			df<-df[,c(1:6,10)]	
		}
		if(recons[x] %in% "TrebGain"){
			df<-df[,c(1:6,11)]	
		}
		if(recons[x] %in% "CyanoGain"){
			df<-df[,c(1:6,12)]	
		}
		if(recons[x] %in% "CephalodiaGain"){
			#print("hey")
			df<-df[,c(1:6,13)]
			#print(head(df))	
		}
		if(recons[x] %in% "AnyCyanoGain"){
			df<-df[,c(1:6,14)]	
		}
		if(recons[x] %in% "AnyCyanoLoss"){
			df<-df[,c(1:6,15)]	
		}
		df$Reconstruction<-recons[x]
		df$Bootstrap<-q
		if(recons[x] %in% "MicroGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "MacroGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "NonGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "TrentGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "TrebGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "CyanoGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "CephalodiaGain"){
			#print("double hey")
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")
			#print(head(df))	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
			#print(tail(cum.df))
		}
		if(recons[x] %in% "AnyCyanoGain"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
		if(recons[x] %in% "AnyCyanoLoss"){
			colnames(df)[7:9]<-c("Change","Reconstruction","Bootstrap")	
			cum.df<-rbind(cum.df,df[df$Change %in% "TRUE" & df$ToAge>0,])
		}
	}
}

cum.df$Reconstruction[cum.df$Reconstruction %in% "MicroGain"]<-"Form: Microlichen"
cum.df$Reconstruction[cum.df$Reconstruction %in% "MacroGain"]<-"Form: Macrolichen"
cum.df$Reconstruction[cum.df$Reconstruction %in% "NonGain"]<-"Photobiont: None"
cum.df$Reconstruction[cum.df$Reconstruction %in% "TrentGain"]<-"Photobiont: Trentepohliales"
cum.df$Reconstruction[cum.df$Reconstruction %in% "TrebGain"]<-"Photobiont: Trebouxiophyceae"
cum.df$Reconstruction[cum.df$Reconstruction %in% "CyanoGain"]<-"Photobiont: Cyanobacteria"
cum.df$Reconstruction[cum.df$Reconstruction %in% "CephalodiaGain"]<-"Cephalodia: Present"
cum.df$Reconstruction[cum.df$Reconstruction %in% "AnyCyanoGain"]<-"Cyanobacteria: Present"

colnames(cum.df)[8]<-"Trait"

cum.df<-cum.df[!cum.df$Trait %in% "Form: Microlichen",]
write.csv(cum.df,file=paste(path,"bootstrap_cum.df.csv",sep=""),row.names=FALSE)



require(ggplot2)
#copy over to local
cum.df<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/bootstrap_cum.df.csv",stringsAsFactors=FALSE)
ml.cum.df<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ml_cum.df.csv",stringsAsFactors=FALSE)

#make ML tree bootstrap replicate 101
ml.cum.df$Bootstrap<-101
cum.df<-rbind(cum.df,ml.cum.df)
cum.df$Kind[cum.df$Bootstrap==101]<-"ML Tree"
cum.df$Kind[cum.df$Bootstrap!=101]<-"Bootstrap Replicate"

timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)
for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255)
}
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$Start<275,]
require(ggplot2)
require(gridExtra)
mycols<-c("lavenderblush4","paleturquoise1","lemonchiffon4","blue","tan","forestgreen","goldenrod2")
timey$Name<-c(NA,"Ng","Pg","K","J","T")

timeplot.cephalodia<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[1],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Gains",title="Cephalodia")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

timeplot.anycyano<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[2],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Gains",title="Cyanobacteria (Photobiont or Cephalodia)")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

timeplot.macro<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[3],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Gains",title="Macrolichen Growth Form")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

timeplot.cyanopb<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[4],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Gains",title="Cyanobacteria Photobiont")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

timeplot.nopb<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[5],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Losses",title="Photobiont Losses")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

timeplot.trebpb<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[6],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Gains",title="Trebouxiophyceae Photobiont")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

timeplot.trentpb<-ggplot(cum.df,aes(x=ToAge,color=Kind))+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect",xmin=timey$End[6],xmax=250,ymin=-Inf,ymax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+theme_bw()+scale_x_reverse(expand=c(0,0),limits=c(250,0))+scale_y_continuous(expand=c(0,0),limits=c(-2,40))+annotate("text",x=timey$Midpoint[6],y=38,label=timey$Name[6],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[5],y=38,label=timey$Name[5],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[4],y=38,label=timey$Name[4],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[3],y=38,label=timey$Name[3],size=3,colour="grey35")+annotate("text",x=timey$Midpoint[2],y=38,label=timey$Name[2],size=3,colour="grey35")+annotate("text",x=1.3,y=38,label=timey$Name[1],size=3,colour="grey35")+scale_colour_manual(values=c(mycols[7],"black"))+labs(x="Node Change Age (Ma)",y="Cumulative Gains",title="Trentepohliales Photobiont")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.23,0.65),legend.title=element_blank())+guides(colour=guide_legend(override.aes=list(lty=c(1,1),alpha=c(0.25,1))))

for(x in 1:101){
	if(x==1){
		cephalodia<-timeplot.cephalodia+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Cephalodia: Present",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)		
		anycyano<-timeplot.anycyano+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Cyanobacteria: Present",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)
		macro<-timeplot.macro+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Form: Macrolichen",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)
		cyanopb<-timeplot.cyanopb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Photobiont: Cyanobacteria",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)
		nopb<-timeplot.nopb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Photobiont: None",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)
		trebpb<-timeplot.trebpb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Photobiont: Trebouxiophyceae",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)
		trentpb<-timeplot.trentpb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==c(x,101),],Trait=="Photobiont: Trentepohliales",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0)
	}
	if(x<101){
		cephalodia<-cephalodia+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Cephalodia: Present",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
		anycyano<-anycyano+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Cyanobacteria: Present",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
		macro<-macro+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Form: Macrolichen",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
		cyanopb<-cyanopb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: Cyanobacteria",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
		nopb<-nopb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: None",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
		trebpb<-trebpb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: Trebouxiophyceae",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
		trentpb<-trentpb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: Trentepohliales",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=0.25)
	}
	if(x==101){
		cephalodia<-cephalodia+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Cephalodia: Present",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)
		anycyano<-anycyano+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Cyanobacteria: Present",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)
		macro<-macro+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Form: Macrolichen",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)
		cyanopb<-cyanopb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: Cyanobacteria",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)
		nopb<-nopb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: None",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)
		trebpb<-trebpb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: Trebouxiophyceae",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)
		trentpb<-trentpb+stat_bin(binwidth=1,data=subset(cum.df[cum.df$Bootstrap==x,],Trait=="Photobiont: Trentepohliales",breaks=seq(-round_any(max(cum.df$ToAge),10,f=ceiling),0,by=1)),aes(y=cumsum(..count..)),geom="step",linetype=1,alpha=1)		
	}
}

#pdf(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/lecanoromycetes_cumulative_trait_changes_bootstraps_w_ml.pdf")
#grid.arrange(cephalodia,anycyano,macro,cyanopb,nopb,trebpb,trentpb,ncol=2)
grid.arrange(nopb,trebpb,trentpb,cyanopb,anycyano,cephalodia,macro,ncol=2)
#dev.off()