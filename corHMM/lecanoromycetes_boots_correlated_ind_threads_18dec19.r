require(ape)
file.path<-"/home/mpnelsen/examl_bootstrapping/"

#ladderize trees
for(q in 31:100){
	tr<-read.tree(file=paste(file.path,"bs_",q,"_trimmed_treepl.tre",sep=""))
	tr<-ladderize(tr,FALSE)
	write.tree(tr,file=paste(file.path,"bs_",q,"_trimmed_treepl_ladderized.tre",sep=""))
}



#ind threads

dat<-read.csv(file="/home/mpnelsen/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)

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

dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[,c("edited_name","multistate")]

testtr<-read.tree(file=paste(file.path,"bs_1_trimmed_treepl_ladderized.tre",sep=""))

dat<-dat[dat$edited_name %in% testtr$tip.label,]

#thallus<-c("micro","macro")
#pb<-c("non","trent","treb","cyano")

#comb<-as.vector(outer(thallus,pb,paste,sep="_"))

#retain only states that are observed
#combs.pres<-comb[c(1,3,5,6,7,8)]

require(corHMM)
combs.pres<-c("micro_non_noceph","micro_trent_noceph","micro_treb_noceph","macro_treb_noceph","micro_cyano_noceph","macro_cyano_noceph","micro_treb_ceph","macro_treb_ceph")


#makes correlated
my.rate.mat.corr<-rate.mat.maker(rate.cat=2,hrm=FALSE,ntraits=1,nstates=length(combs.pres),model="ARD")
rownames(my.rate.mat.corr)<-0:7
colnames(my.rate.mat.corr)<-0:7
rownames(my.rate.mat.corr)<-combs.pres
colnames(my.rate.mat.corr)<-combs.pres
rownames(my.rate.mat.corr)<-0:7
colnames(my.rate.mat.corr)<-0:7

#include 47,55,34,42
my.rate.mat.corr.cleanerer<-rate.par.drop(rate.mat.index=my.rate.mat.corr,drop.par=c(22,36,43,50,23,37,44,51,38,52,3,10,32,46,25,54,5,12,19,48,6,13,27,41,7,14,21,35))
rownames(my.rate.mat.corr.cleanerer)<-0:7
colnames(my.rate.mat.corr.cleanerer)<-0:7

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

my.rate.mat.corr.cleanerer.more.more<-my.rate.mat.corr.cleanerer.more
#make all shifts from treb to cyano equal (regardless of ceph)
my.rate.mat.corr.cleanerer.more.more<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.more.more,eq.par=c(15,16))

#make all shifts from cyano to treb equal (regardless of ceph)
my.rate.mat.corr.cleanerer.more.more<-rate.par.eq(rate.mat.index=my.rate.mat.corr.cleanerer.more.more,eq.par=c(10,17))

my.rate.mat.corr.cleanerer.constr.constr<-my.rate.mat.corr.cleanerer.more.more
rownames(my.rate.mat.corr.cleanerer.constr.constr)<-0:7
colnames(my.rate.mat.corr.cleanerer.constr.constr)<-0:7

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

#First...slow
for(q in 65:100){
	tr<-read.tree(file=paste(file.path,"bs_",q,"_trimmed_treepl_ladderized.tre",sep=""))
	my.rate.mat.corr.cleanerer.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
	save(my.rate.mat.corr.cleanerer.fitz.no.diagn,file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
}

#Second...slow
for(q in 65:100){
	tr<-read.tree(file=paste(file.path,"bs_",q,"_trimmed_treepl_ladderized.tre",sep=""))
	my.rate.mat.corr.cleanerer.constr.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.constr,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
	save(my.rate.mat.corr.cleanerer.constr.fitz.no.diagn,file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.constr.fitz.no.diagn.Rsave",sep=""))
}

#Third...fast
for(q in 30:100){
	tr<-read.tree(file=paste(file.path,"bs_",q,"_trimmed_treepl_ladderized.tre",sep=""))
	my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.eq.notceph.more,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
	save(my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn,file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn.Rsave",sep=""))
}

#Fourth...slow
for(q in 65:100){
	tr<-read.tree(file=paste(file.path,"bs_",q,"_trimmed_treepl_ladderized.tre",sep=""))
	my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.constr.constr,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
	save(my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn,file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn.Rsave",sep=""))
}

#Fifth...fast
for(q in 30:100){
	tr<-read.tree(file=paste(file.path,"bs_",q,"_trimmed_treepl_ladderized.tre",sep=""))
	my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn<-rayDISC(phy=tr,data=dat,ntraits=1,rate.mat= my.rate.mat.corr.cleanerer.eq.cephnocephsame.more,node.states="marginal",diagn=FALSE,root.p="maddfitz",ub=1000000,lb=-1000000,model="ARD")
	save(my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn,file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn.Rsave",sep=""))
}




#Summary
cns<-c("Tree","my_rate_mat_corr_cleanerer_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_constr_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_eq_notceph_more_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_constr_constr_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_eq_cephnocephsame_more_fitz_no_diagn_AICc","Best")
sum.df<-as.data.frame(matrix(nrow=100,ncol=length(cns)))
colnames(sum.df)<-cns
sum.df$Tree<-1:100

for(q in 1:100){
	#my.rate.mat.corr.cleanerer.fitz.no.diagn<-load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	#my.rate.mat.corr.cleanerer.constr.fitz.no.diagn<-load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.constr.fitz.no.diagn.Rsave",sep=""))
	#my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn<-load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn.Rsave",sep=""))
	#my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn<-load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn.Rsave",sep=""))
	#my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn<-load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn.Rsave",sep=""))
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.constr.fitz.no.diagn.Rsave",sep=""))
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn.Rsave",sep=""))
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn.Rsave",sep=""))
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn.Rsave",sep=""))
	sum.df$my_rate_mat_corr_cleanerer_fitz_no_diagn_AICc[q]<-my.rate.mat.corr.cleanerer.fitz.no.diagn$AICc
	sum.df$my_rate_mat_corr_cleanerer_constr_fitz_no_diagn_AICc[q]<-my.rate.mat.corr.cleanerer.constr.fitz.no.diagn$AICc
	sum.df$my_rate_mat_corr_cleanerer_eq_notceph_more_fitz_no_diagn_AICc[q]<-my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn$AICc
	sum.df$my_rate_mat_corr_cleanerer_constr_constr_fitz_no_diagn_AICc[q]<-my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn$AICc
	sum.df$my_rate_mat_corr_cleanerer_eq_cephnocephsame_more_fitz_no_diagn_AICc[q]<-my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn$AICc
	sum.df[q,"Best"]<-c("my_rate_mat_corr_cleanerer_fitz_no_diagn","my_rate_mat_corr_cleanerer_constr_fitz_no_diagn","my_rate_mat_corr_cleanerer_eq_notceph_more_fitz_no_diagn","my_rate_mat_corr_cleanerer_constr_constr_fitz_no_diagn","my_rate_mat_corr_cleanerer_eq_cephnocephsame_more_fitz_no_diagn")[sum.df[q,c("my_rate_mat_corr_cleanerer_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_constr_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_eq_notceph_more_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_constr_constr_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_eq_cephnocephsame_more_fitz_no_diagn_AICc")] %in% min(sum.df[q,c("my_rate_mat_corr_cleanerer_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_constr_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_eq_notceph_more_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_constr_constr_fitz_no_diagn_AICc","my_rate_mat_corr_cleanerer_eq_cephnocephsame_more_fitz_no_diagn_AICc")])]
}

table(sum.df$Best)
write.csv(sum.df,file=paste(file.path,"_boot_models_summary.csv",sep=""))

#my.rate.mat.corr.cleanerer.fitz.no.diagn




require(corHMM)
file.path<-"/home/mpnelsen/examl_bootstrapping/"

my10<-NULL
my20<-NULL
my40<-NULL
my01<-NULL
my21<-NULL
my41<-NULL
my02<-NULL
my12<-NULL
my32<-NULL
my42<-NULL
my62<-NULL
my23<-NULL
my53<-NULL
my73<-NULL
my04<-NULL
my14<-NULL
my24<-NULL
my54<-NULL
my64<-NULL
my35<-NULL
my45<-NULL
my75<-NULL
my26<-NULL
my46<-NULL
my76<-NULL
my37<-NULL
my57<-NULL
my67<-NULL

for(q in 1:100){
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	my10<-c(my10,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[2,1])
	my20<-c(my20,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[3,1])
	my40<-c(my40,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[5,1])
	my01<-c(my01,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[1,2])
	my21<-c(my21,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[3,2])
	my41<-c(my41,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[5,2])
	my02<-c(my02,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[1,3])
	my12<-c(my12,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[2,3])
	my32<-c(my32,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[4,3])
	my42<-c(my42,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[5,3])
	my62<-c(my62,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[7,3])
	my23<-c(my23,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[3,4])
	my53<-c(my53,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[6,4])
	my73<-c(my73,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[8,4])
	my04<-c(my04,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[1,5])
	my14<-c(my14,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[2,5])
	my24<-c(my24,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[3,5])
	my54<-c(my54,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[6,5])
	my64<-c(my64,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[7,5])
	my35<-c(my35,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[4,6])
	my45<-c(my45,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[5,6])
	my75<-c(my75,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[8,6])
	my26<-c(my26,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[3,7])
	my46<-c(my46,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[5,7])
	my76<-c(my76,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[8,7])
	my37<-c(my37,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[4,8])
	my57<-c(my57,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[6,8])
	my67<-c(my67,my.rate.mat.corr.cleanerer.fitz.no.diagn$solution[7,8])
}



my.means<-as.data.frame(matrix(nrow=8,ncol=8))
colnames(my.means)<-0:7
rownames(my.means)<-0:7
my.means[2,1]<-mean(my10)
my.means[3,1]<-mean(my20)
my.means[5,1]<-mean(my40)
my.means[1,2]<-mean(my01)
my.means[3,2]<-mean(my21)
my.means[5,2]<-mean(my41)
my.means[1,3]<-mean(my02)
my.means[2,3]<-mean(my12)
my.means[4,3]<-mean(my32)
my.means[5,3]<-mean(my42)
my.means[7,3]<-mean(my62)
my.means[3,4]<-mean(my23)
my.means[6,4]<-mean(my53)
my.means[8,4]<-mean(my73)
my.means[1,5]<-mean(my04)
my.means[2,5]<-mean(my14)
my.means[3,5]<-mean(my24)
my.means[6,5]<-mean(my54)
my.means[7,5]<-mean(my64)
my.means[4,6]<-mean(my35)
my.means[5,6]<-mean(my45)
my.means[8,6]<-mean(my75)
my.means[3,7]<-mean(my26)
my.means[5,7]<-mean(my46)
my.means[8,7]<-mean(my76)
my.means[4,8]<-mean(my37)
my.means[6,8]<-mean(my57)
my.means[7,8]<-mean(my67)
write.csv(my.means,file=paste(file.path,"_boot_my.rate.mat.corr.cleanerer.fitz.no.diagn_my.means.csv",sep=""))


my.mins<-as.data.frame(matrix(nrow=8,ncol=8))
colnames(my.mins)<-0:7
rownames(my.mins)<-0:7
my.mins[2,1]<-min(my10)
my.mins[3,1]<-min(my20)
my.mins[5,1]<-min(my40)
my.mins[1,2]<-min(my01)
my.mins[3,2]<-min(my21)
my.mins[5,2]<-min(my41)
my.mins[1,3]<-min(my02)
my.mins[2,3]<-min(my12)
my.mins[4,3]<-min(my32)
my.mins[5,3]<-min(my42)
my.mins[7,3]<-min(my62)
my.mins[3,4]<-min(my23)
my.mins[6,4]<-min(my53)
my.mins[8,4]<-min(my73)
my.mins[1,5]<-min(my04)
my.mins[2,5]<-min(my14)
my.mins[3,5]<-min(my24)
my.mins[6,5]<-min(my54)
my.mins[7,5]<-min(my64)
my.mins[4,6]<-min(my35)
my.mins[5,6]<-min(my45)
my.mins[8,6]<-min(my75)
my.mins[3,7]<-min(my26)
my.mins[5,7]<-min(my46)
my.mins[8,7]<-min(my76)
my.mins[4,8]<-min(my37)
my.mins[6,8]<-min(my57)
my.mins[7,8]<-min(my67)
write.csv(my.mins,file=paste(file.path,"_boot_my.rate.mat.corr.cleanerer.fitz.no.diagn_my.mins.csv",sep=""))


my.maxs<-as.data.frame(matrix(nrow=8,ncol=8))
colnames(my.maxs)<-0:7
rownames(my.maxs)<-0:7
my.maxs[2,1]<-max(my10)
my.maxs[3,1]<-max(my20)
my.maxs[5,1]<-max(my40)
my.maxs[1,2]<-max(my01)
my.maxs[3,2]<-max(my21)
my.maxs[5,2]<-max(my41)
my.maxs[1,3]<-max(my02)
my.maxs[2,3]<-max(my12)
my.maxs[4,3]<-max(my32)
my.maxs[5,3]<-max(my42)
my.maxs[7,3]<-max(my62)
my.maxs[3,4]<-max(my23)
my.maxs[6,4]<-max(my53)
my.maxs[8,4]<-max(my73)
my.maxs[1,5]<-max(my04)
my.maxs[2,5]<-max(my14)
my.maxs[3,5]<-max(my24)
my.maxs[6,5]<-max(my54)
my.maxs[7,5]<-max(my64)
my.maxs[4,6]<-max(my35)
my.maxs[5,6]<-max(my45)
my.maxs[8,6]<-max(my75)
my.maxs[3,7]<-max(my26)
my.maxs[5,7]<-max(my46)
my.maxs[8,7]<-max(my76)
my.maxs[4,8]<-max(my37)
my.maxs[6,8]<-max(my57)
my.maxs[7,8]<-max(my67)
write.csv(my.maxs,file=paste(file.path,"_boot_my.rate.mat.corr.cleanerer.fitz.no.diagn_my.maxs.csv",sep=""))


my.sds<-as.data.frame(matrix(nrow=8,ncol=8))
colnames(my.sds)<-0:7
rownames(my.sds)<-0:7
my.sds[2,1]<-sd(my10)
my.sds[3,1]<-sd(my20)
my.sds[5,1]<-sd(my40)
my.sds[1,2]<-sd(my01)
my.sds[3,2]<-sd(my21)
my.sds[5,2]<-sd(my41)
my.sds[1,3]<-sd(my02)
my.sds[2,3]<-sd(my12)
my.sds[4,3]<-sd(my32)
my.sds[5,3]<-sd(my42)
my.sds[7,3]<-sd(my62)
my.sds[3,4]<-sd(my23)
my.sds[6,4]<-sd(my53)
my.sds[8,4]<-sd(my73)
my.sds[1,5]<-sd(my04)
my.sds[2,5]<-sd(my14)
my.sds[3,5]<-sd(my24)
my.sds[6,5]<-sd(my54)
my.sds[7,5]<-sd(my64)
my.sds[4,6]<-sd(my35)
my.sds[5,6]<-sd(my45)
my.sds[8,6]<-sd(my75)
my.sds[3,7]<-sd(my26)
my.sds[5,7]<-sd(my46)
my.sds[8,7]<-sd(my76)
my.sds[4,8]<-sd(my37)
my.sds[6,8]<-sd(my57)
my.sds[7,8]<-sd(my67)
write.csv(my.sds,file=paste(file.path,"_boot_my.rate.mat.corr.cleanerer.fitz.no.diagn_my.sds.csv",sep=""))


my.meds<-as.data.frame(matrix(nrow=8,ncol=8))
colnames(my.meds)<-0:7
rownames(my.meds)<-0:7
my.meds[2,1]<-median(my10)
my.meds[3,1]<-median(my20)
my.meds[5,1]<-median(my40)
my.meds[1,2]<-median(my01)
my.meds[3,2]<-median(my21)
my.meds[5,2]<-median(my41)
my.meds[1,3]<-median(my02)
my.meds[2,3]<-median(my12)
my.meds[4,3]<-median(my32)
my.meds[5,3]<-median(my42)
my.meds[7,3]<-median(my62)
my.meds[3,4]<-median(my23)
my.meds[6,4]<-median(my53)
my.meds[8,4]<-median(my73)
my.meds[1,5]<-median(my04)
my.meds[2,5]<-median(my14)
my.meds[3,5]<-median(my24)
my.meds[6,5]<-median(my54)
my.meds[7,5]<-median(my64)
my.meds[4,6]<-median(my35)
my.meds[5,6]<-median(my45)
my.meds[8,6]<-median(my75)
my.meds[3,7]<-median(my26)
my.meds[5,7]<-median(my46)
my.meds[8,7]<-median(my76)
my.meds[4,8]<-median(my37)
my.meds[6,8]<-median(my57)
my.meds[7,8]<-median(my67)
write.csv(my.meds,file=paste(file.path,"_boot_my.rate.mat.corr.cleanerer.fitz.no.diagn_my.meds.csv",sep=""))






require(diversitree)
require(corHMM)
file.path<-"/home/mpnelsen/examl_bootstrapping/"
p2p<-diversitree:::plot2.phylo

for(q in 1:100){
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	tree<-my.rate.mat.corr.cleanerer.fitz.no.diagn$phy
	thtr<-tree
	png(paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn_NODES.png",sep=""),width=14,height=12.51,units="in",res=600)
	t <- max(branching.times(thtr))
	obj<-p2p(thtr,type="fan",label.offset = t * 1/4, show.tip.label = FALSE)
	cbPalette <- c("white", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999", "#D55E00", "#CC79A7","#9999CC","#0072B2")
	nodelabels(pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[,1:8],cex=0.2,piecol=cbPalette[1:8])
	title(main=paste("Tree ",q))
	dev.off()
}





require(diversitree)
require(corHMM)
file.path<-"/home/mpnelsen/examl_bootstrapping/"
p2p<-diversitree:::plot2.phylo
div.fa<-diversitree:::filled.arcs
dat<-read.csv(file="/home/mpnelsen/examl_bootstrapping/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
traits<-c("Form","Photobiont","Cephalodia")

pdf(paste(file.path,"BOOTS.my.rate.mat.corr.cleanerer.fitz.no.diagn_NODES.pdf",sep=""))
for(q in 1:100){
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	tree<-my.rate.mat.corr.cleanerer.fitz.no.diagn$phy
	thtr<-tree
	t <- max(branching.times(thtr))
	tiplabs<-thtr$tip.label
	dat.c<-dat[dat$edited_name %in% tiplabs,]
	dat.c<-dat.c[match(tiplabs,rownames(dat.c)),]
	dat.c<-dat.c[,c(1:11,35,30,21)]
	colnames(dat.c)[12:14]<-c("Form","Photobiont","Cephalodia")
	dat.c$Photobiont[dat.c$Photobiont %in% "0&2"]<-4
	dat.c$Photobiont[dat.c$Photobiont %in% "2&3"]<-5
	for(x in 12:14){
		dat.c[,x]<-as.numeric(dat.c[,x])
	}
	obj<-p2p(thtr,type="fan",label.offset = t * 1/4, show.tip.label = TRUE, edge.width=0.25, cex=0.045)
	xy<-obj$xy
	theta<-xy$theta[seq_along(thtr$tip.label)]
	dt<-diff(sort(theta))[1]/2
	cbPalette <- c("white", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999", "#D55E00", "#CC79A7","#9999CC","#0072B2")
	st.cols<-c("lemonchiffon1","lemonchiffon4","white","goldenrod2","forestgreen","blue", "light green","purple","darkolivegreen1","lavenderblush4")
	cols=list("Form"=c(st.cols[1],st.cols[2]),"Photobiont"=c(st.cols[3],st.cols[4],st.cols[5],st.cols[6],st.cols[7],st.cols[8]),"Cephalodia"=c(st.cols[9],st.cols[10]))
	dat.plot<-dat.c[,c(12:14)]
	w = 10
	for (i in seq_along(cols)){
		idx <- dat.plot[[names(dat.plot)[i]]]
		if (any(idx == 0, na.rm = TRUE)){
			idx <- idx + 1
		}	
		div.fa(theta-dt,theta+dt,rep(250+(i*0.02)+(i*w),length(theta)),w,cols[[i]][idx])
	}
	legend(x=-400,y=-250,title=NA,c(expression(bold("Growth Form")),"Micro","Micro","Micro","Macro","Micro","Macro","Micro","Macro"),col=c(NA,cbPalette),pch=19,bty="n",cex=0.25,pt.cex=0.5)
	legend(x=-350,y=-250,title=NA,c(expression(bold("Primary Photobiont")),"None","Trentepohliales","Trebouxiophyceae","Trebouxiophyceae","Cyanobacteria","Cyanobacteria","Trebouxiophyceae","Trebouxiophyceae"),col=NA,pch=19,bty="n",cex=0.25,pt.cex=0.5)
	legend(x=-250,y=-250,title=NA,c(expression(bold("Cephalodia")),"None","None","None","None","None","None","Present","Present"),col=NA,pch=19,bty="n",cex=0.25,pt.cex=0.5)
	segments(x0=-400,y0=-260,x1=-208,y1=-260)
	text(x=-300,y=-248,expression(bold("Node States")),cex=0.5)
	legend(x=300,y=-250,title=NA,c(expression(bold("Growth Form (inner)")),"Micro","Macro", NA,expression(bold("Primary Photobiont (middle)")),"None","Trentepohliales","Trebouxiophyceae","Cyanobacteria","None & Trebouxiophyceae","Trebouxiophyceae & Cyanobacteria",NA,expression(bold("Cephalodia (outer)")),"None","Present"),col=c(NA,"lemonchiffon1","lemonchiffon4",NA,NA,"white","goldenrod2","forestgreen","blue","light green","purple",NA,NA,"darkolivegreen1","lavenderblush4"),pch=19,bty="n",cex=0.25,pt.cex=0.5)
	segments(x0=300,y0=-260,x1=400,y1=-260)
	text(x=350,y=-248,expression(bold("Tip States")),cex=0.5)
	nodelabels(pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[,1:8],cex=0.2,piecol=cbPalette[1:8])
	title(main=paste("Bootstrap Topology ",q))
}
dev.off()








#tried re-plotting bootstraps on 18 Dec 2019


require(diversitree)
require(corHMM)
file.path<-"/home/mpnelsen/examl_bootstrapping/"
p2p<-diversitree:::plot2.phylo
div.fa<-diversitree:::filled.arcs
dat<-read.csv(file="/home/mpnelsen/examl_bootstrapping/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
traits<-c("Form","Photobiont","Cephalodia")

pdf(paste(file.path,"BOOTS.my.rate.mat.corr.cleanerer.fitz.no.diagn_NODES_18dec19.pdf",sep=""),width=11,height=11)
for(q in 1:100){
	load(file=paste(file.path,q,"_my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave",sep=""))
	tree<-my.rate.mat.corr.cleanerer.fitz.no.diagn$phy
	thtr<-tree
	t <- max(branching.times(thtr))
	tiplabs<-thtr$tip.label
	dat.c<-dat[dat$edited_name %in% tiplabs,]
	dat.c<-dat.c[match(tiplabs,rownames(dat.c)),]
	dat.c<-dat.c[,c(1:11,35,30,21)]
	colnames(dat.c)[12:14]<-c("Form","Photobiont","Cephalodia")
	dat.c$Photobiont[dat.c$Photobiont %in% "0&2"]<-4
	dat.c$Photobiont[dat.c$Photobiont %in% "2&3"]<-5
	for(x in 12:14){
		dat.c[,x]<-as.numeric(dat.c[,x])
	}
	obj<-p2p(thtr,type="fan",label.offset = t * 1/5, show.tip.label = TRUE, edge.width=0.25, cex=0.042)
	xy<-obj$xy
	theta<-xy$theta[seq_along(thtr$tip.label)]
	dt<-diff(sort(theta))[1]/2
	cbPalette <- c("white", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999", "#D55E00", "#CC79A7","#9999CC","#0072B2")
	st.cols<-c("lemonchiffon1","lemonchiffon4","white","goldenrod2","forestgreen","blue", "light green","purple","darkolivegreen1","lavenderblush4")
	cols=list("Form"=c(st.cols[1],st.cols[2]),"Photobiont"=c(st.cols[3],st.cols[4],st.cols[5],st.cols[6],st.cols[7],st.cols[8]),"Cephalodia"=c(st.cols[9],st.cols[10]))
	dat.plot<-dat.c[,c(12:14)]
	w = 10
	for (i in seq_along(cols)){
		idx <- dat.plot[[names(dat.plot)[i]]]
		if (any(idx == 0, na.rm = TRUE)){
			idx <- idx + 1
		}	
		div.fa(theta-dt,theta+dt,rep(250+(i*0.02)+(i*w),length(theta)),w,cols[[i]][idx])
	}
	legend(x=-350,y=-250,title=NA,c(expression(bold("Growth Form")),"Micro","Micro","Micro","Macro","Micro","Macro","Micro","Macro"),col=c(NA,cbPalette),pch=19,bty="n",cex=0.25,pt.cex=0.5)
	legend(x=-300,y=-250,title=NA,c(expression(bold("Primary Photobiont")),"None","Trentepohliales","Trebouxiophyceae","Trebouxiophyceae","Cyanobacteria","Cyanobacteria","Trebouxiophyceae","Trebouxiophyceae"),col=NA,pch=19,bty="n",cex=0.25,pt.cex=0.5)
	legend(x=-250,y=-250,title=NA,c(expression(bold("Cephalodia")),"None","None","None","None","None","None","Present","Present"),col=NA,pch=19,bty="n",cex=0.25,pt.cex=0.5)
	segments(x0=-350,y0=-260,x1=-208,y1=-260)
	text(x=-275,y=-248,expression(bold("Node States")),cex=0.5)
	legend(x=250,y=-260,title=NA,c(expression(bold("Growth Form (inner)")),"Micro","Macro", NA,expression(bold("Primary Photobiont (middle)")),"None","Trentepohliales","Trebouxiophyceae","Cyanobacteria","None & Trebouxiophyceae","Trebouxiophyceae & Cyanobacteria",NA,expression(bold("Cephalodia (outer)")),"None","Present"),col=c(NA,"lemonchiffon1","lemonchiffon4",NA,NA,"white","goldenrod2","forestgreen","blue","light green","purple",NA,NA,"darkolivegreen1","lavenderblush4"),pch=19,bty="n",cex=0.25,pt.cex=0.5)
	segments(x0=250,y0=-260,x1=300,y1=-260)
	text(x=275,y=-248,expression(bold("Tip States")),cex=0.5)
	nodelabels(pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[,1:8],cex=0.2,piecol=cbPalette[1:8])
	title(main=paste("Bootstrap Topology ",q))
}
dev.off()
