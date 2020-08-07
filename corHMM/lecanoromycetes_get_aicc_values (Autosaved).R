

#my.rate.mat.corr.cleanerer.fitz.no.diagn,
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave")

#my.rate.mat.corr.cleanerer.constr.fitz.no.diagn,
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.constr.fitz.no.diagn.Rsave")

#my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn,
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn.Rsave")

#my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn,
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn.Rsave")

#my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn,
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn.Rsave")


#my.rate.mat.corr.cleanerer.fitz.no.diagn
#my.rate.mat.corr.cleanerer.constr.fitz.no.diagn
#my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn
#my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn
#my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn


my.rate.mat.corr.cleanerer.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.constr.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.eq.notceph.more.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.constr.constr.fitz.no.diagn$AICc
my.rate.mat.corr.cleanerer.eq.cephnocephsame.more.fitz.no.diagn$AICc


attributes(my.rate.mat.corr.cleanerer.fitz.no.diagn)
traits<-as.matrix(my.rate.mat.corr.cleanerer.fitz.no.diagn$data$multistate)
rownames(traits)<-my.rate.mat.corr.cleanerer.fitz.no.diagn$data$edited_name
colnames(traits)<-"multistate"
traits[,"multistate"][traits[,"multistate"] %in% "0"]<-as.numeric(8)
traits[,"multistate"][traits[,"multistate"] %in% "1"]<-as.numeric(1)
traits[,"multistate"][traits[,"multistate"] %in% "2"]<-as.numeric(2)
traits[,"multistate"][traits[,"multistate"] %in% "3"]<-as.numeric(3)
traits[,"multistate"][traits[,"multistate"] %in% "4"]<-as.numeric(4)
traits[,"multistate"][traits[,"multistate"] %in% "5"]<-as.numeric(5)
traits[,"multistate"][traits[,"multistate"] %in% "6"]<-as.numeric(6)
traits[,"multistate"][traits[,"multistate"] %in% "7"]<-as.numeric(7)

#just do this for now
traits[,"multistate"][traits[,"multistate"] %in% "0&2"]<-as.numeric(8)
traits[,"multistate"][traits[,"multistate"] %in% "5&7"]<-as.numeric(5)
traits.num<-as.numeric(traits)
names(traits.num)<-row.names(traits)
traits.num<-as.matrix(traits.num)

my.dtt<-dtt(phy=tree,data=traits.num,index="num.states",nsim=10)



geo$dat[geo$dat[,"wingL"]>4.3,"wingL"]<-"long"
geo$dat[geo$dat[,"wingL"]<=4.3,"wingL"]<-"short"
geo$dat[,"wingL"]
disparity(phy=geo$phy, data=geo$dat[,"wingL"],index="num.states")