require(corHMM)

tr<-read.tree(file="/home/mpnelsen/lecanoromycetes_correlated_evolution/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")

dat<-read.csv(file="/home/mpnelsen/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[,c("edited_name","multistate")]
dat<-dat[dat$edited_name %in% tr$tip.label,]

dat$multistate<-NA
for(x in 1:nrow(dat)){
	#micro, non-lichenized	
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==0){
		dat$multistate[x]<-0
	}
	#micro, Trent
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==1){
		dat$multistate[x]<-1
	}
	#micro, Treb
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==2){
		dat$multistate[x]<-2
	}
	#macro, Treb
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]==2){
		dat$multistate[x]<-3
	}
	#micro, Cyano
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==3){
		dat$multistate[x]<-4
	}
	#macro, Cyano
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]==3){
		dat$multistate[x]<-5
	}
	#macro, Treb and Cyano
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]=="2&3"){
		dat$multistate[x]<-"3&5"
	}
	#micro, non-lichenized and Treb
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]=="0&2"){
		dat$multistate[x]<-"0&2"
	}
}


dat<-dat[,c("edited_name","multistate")]

thallus<-c("micro","macro")
pb<-c("non","trent","treb","cyano")

comb<-as.vector(outer(thallus,pb,paste,sep="_"))

#retain only states that are observed
combs.pres<-comb[c(1,3,5,6,7,8)]

#makes correlated
my.rate.mat.corr<-rate.mat.maker(rate.cat=2,hrm=FALSE,ntraits=1,nstates=length(combs.pres),model="ARD")
my.rate.mat.corr<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(16,26,17,27,28,3,8,24,19,5,10,15))
rownames(my.rate.mat.corr)<-0:5
colnames(my.rate.mat.corr)<-0:5
my.rate.mat.corr
my.rate.mat.corr.fitz<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
my.rate.mat.corr.no15.no<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(1,5))
my.rate.mat.corr.fitzno15.no<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat=my.rate.mat.corr.no15.no,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

my.rate.mat.corr.no135.no<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(1,3,5))
my.rate.mat.corr.fitzno135.no<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat=my.rate.mat.corr.no135.no,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

my.rate.mat.corr.no1356.no<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(1,3,5,6))
my.rate.mat.corr.fitzno1356.no<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat=my.rate.mat.corr.no1356.no,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

my.rate.mat.corr.no136.no<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(1,3,6))
my.rate.mat.corr.fitzno136.no<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat=my.rate.mat.corr.no136.no,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

my.rate.mat.corr.no13.no<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(1,3))
my.rate.mat.corr.fitzno13.no<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat=my.rate.mat.corr.no13.no,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")














#cephalodia

tr<-read.tree(file="/home/mpnelsen/lecanoromycetes_correlated_evolution/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")

dat<-read.csv(file="/home/mpnelsen/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[,c("edited_name","multistate")]
dat<-dat[dat$edited_name %in% tr$tip.label,]

dat$multistate<-NA
for(x in 1:nrow(dat)){
	#micro, non-lichenized, no ceph	
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==0 & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-0
	}
	#micro, Trent, no ceph
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==1 & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-1
	}
	#micro, Treb, no ceph
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==2 & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-2
	}
	#micro, Treb, ceph
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==2 & dat$TripartitePossible[x]==1){
		dat$multistate[x]<-6
	}
	#macro, Treb, no ceph
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]==2 & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-3
	}
	#macro, Treb, ceph
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]==2 & dat$TripartitePossible[x]==1){
		dat$multistate[x]<-7
	}
	#micro, Cyano, no ceph
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]==3 & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-4
	}
	#macro, Cyano, no ceph
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]==3 & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-5
	}
	#macro, Treb no ceph, and Cyano
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]=="2&3" & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-"3&5"
	}
	#macro, Treb ceph, and Cyano
	if(dat$Micro[x]==0 & dat$Macro[x]==1 & dat$NutMult[x]=="2&3" & dat$TripartitePossible[x]==1){
		dat$multistate[x]<-"5&7"
	}
	#micro, non-lichenized and Treb no cepph
	if(dat$Micro[x]==1 & dat$Macro[x]==0 & dat$NutMult[x]=="0&2" & dat$TripartitePossible[x]==0){
		dat$multistate[x]<-"0&2"
	}
}



dat<-dat[,c("edited_name","multistate")]

#thallus<-c("micro","macro")
#pb<-c("non","trent","treb","cyano")

#comb<-as.vector(outer(thallus,pb,paste,sep="_"))

#retain only states that are observed
#combs.pres<-comb[c(1,3,5,6,7,8)]

combs.pres<-c("micro_non_noceph","micro_trent_noceph","micro_treb_noceph","macro_treb_noceph","micro_cyano_noceph","macro_cyano_noceph","micro_treb_ceph","macro_treb_ceph")


#makes correlated
my.rate.mat.corr<-rate.mat.maker(rate.cat=2,hrm=FALSE,ntraits=1,nstates=length(combs.pres),model="ARD")
rownames(my.rate.mat.corr)<-0:7
colnames(my.rate.mat.corr)<-0:7
rownames(my.rate.mat.corr)<-combs.pres
colnames(my.rate.mat.corr)<-combs.pres
rownames(my.rate.mat.corr)<-0:7
colnames(my.rate.mat.corr)<-0:7

#47,55,34,42??
my.rate.mat.corr.cleaner<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(22,36,43,50,23,37,44,51,38,52,3,10,32,46,25,47,54,5,12,19,48,55,6,13,27,34,41,7,14,21,35,42))
rownames(my.rate.mat.corr.cleaner)<-0:7
colnames(my.rate.mat.corr.cleaner)<-0:7
my.rate.mat.corr.cleaner.fitz<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleaner,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")


#include 47,55,34,42
my.rate.mat.corr.cleanerer<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(22,36,43,50,23,37,44,51,38,52,3,10,32,46,25,54,5,12,19,48,6,13,27,41,7,14,21,35))
rownames(my.rate.mat.corr.cleanerer)<-0:7
colnames(my.rate.mat.corr.cleanerer)<-0:7
my.rate.mat.corr.cleanerer.fitz<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer,node.states="marginal",diagn=TRUE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
my.rate.mat.corr.cleanerer.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")



#make micro to macro equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer,eq.par=c(12,21,28))

#make macro to micro equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(9,18,24))

#make treb no ceph to cyano equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(17,19))

#make cyano to treb no ceph equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(10,13))

#make treb ceph to cyano equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(17,18))

#make cyano to treb ceph equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(19,21))

#make treb ceph to treb no ceph equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(11,13))

#make treb no ceph to treb ceph equal
my.rate.mat.corr.cleanerer.more<-rate.par.eq(rate.mat.index= my.rate.mat.corr.cleanerer.more,eq.par=c(17,19))

my.rate.mat.corr.cleanerer.constr<-my.rate.mat.corr.cleanerer.more

rownames(my.rate.mat.corr.cleanerer.constr)<-0:7
colnames(my.rate.mat.corr.cleanerer.constr)<-0:7
my.rate.mat.corr.cleanerer.constr.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.constr,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

	#make macro/micro equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.more,eq.par=c(12,9))
	#make non/trent equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(4,1))
	#make treb no ceph/non equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(6,2))
	#make cyano/non equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(10,3))
	#make cyano/trent equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(10,5))
	#make treb/trent equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(6,4))
	#make treb no ceph/cyano equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(9,7))
	#make treb ceph/cyano equal
	my.rate.mat.corr.cleanerer.eq.notceph<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.notceph,eq.par=c(9,11))
	my.rate.mat.corr.cleanerer.eq.notceph.more<-my.rate.mat.corr.cleanerer.eq.notceph
	rownames(my.rate.mat.corr.cleanerer.eq.notceph.more)<-0:7
	colnames(my.rate.mat.corr.cleanerer.eq.notceph.more)<-0:7
	my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.eq.notceph.more,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

my.rate.mat.corr.cleanerer.more.more<-my.rate.mat.corr.cleanerer.more
#make all shifts from treb to cyano equal (regardless of ceph)
my.rate.mat.corr.cleanerer.more.more<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.more.more,eq.par=c(15,16))

#make all shifts from cyano to treb equal (regardless of ceph)
my.rate.mat.corr.cleanerer.more.more<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.more.more,eq.par=c(10,17))

my.rate.mat.corr.cleanerer.constr.constr<-my.rate.mat.corr.cleanerer.more.more
rownames(my.rate.mat.corr.cleanerer.constr.constr)<-0:7
colnames(my.rate.mat.corr.cleanerer.constr.constr)<-0:7
my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.constr.constr,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

	#make macro/micro equal
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.more.more,eq.par=c(12,9))
	#make non/trent equal
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.cephnocephsame,eq.par=c(4,1))
	#make treb no ceph/non equal
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.cephnocephsame,eq.par=c(2,6))
	#make cyano/non equal
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.cephnocephsame,eq.par=c(3,10))
	#make cyano/trent equal
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.cephnocephsame,eq.par=c(5,10))
	#make treb/trent equal
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.cephnocephsame,eq.par=c(6,4))
	#make treb to cyano equal regardless of ceph
	my.rate.mat.corr.cleanerer.eq.cephnocephsame<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.eq.cephnocephsame,eq.par=c(9,7))

	my.rate.mat.corr.cleanerer.eq.cephnocephsame.more<-my.rate.mat.corr.cleanerer.eq.cephnocephsame
	rownames(my.rate.mat.corr.cleanerer.eq.cephnocephsame.more)<-0:7
	colnames(my.rate.mat.corr.cleanerer.eq.cephnocephsame.more)<-0:7
	my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.eq.cephnocephsame.more,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")

my.rate.mat.corr.cleanerer.fitz.no.diagn
my.rate.mat.corr.cleanerer.constr.fitz.no.diagn
my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn
my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn
my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn


my.rate.mat.corr.cleanerer.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.constr.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn$AICc

save(my.rate.mat.corr.cleanerer.fitz.no.diagn,file="~/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave")
save(my.rate.mat.corr.cleanerer.constr.fitz.no.diagn,file="~/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.constr.fitz.no.diagn.Rsave")
save(my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn,file="~/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn.Rsave")
save(my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn,file="~/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn.Rsave")
save(my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn,file="~/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn.Rsave")

