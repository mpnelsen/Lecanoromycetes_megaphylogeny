
require(diversitree)
require(corHMM)
p2p<-diversitree:::plot2.phylo
div.fa<-diversitree:::filled.arcs
dat<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lecanoromycetes_photobiont_growth_form/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
traits<-c("Form","Photobiont","Cephalodia")

png("/Users/matthewnelsen/Documents/papers_reviews/papers/lecanoromycetes_photobiont_growth_form/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ML.my.rate.mat.corr.cleanerer.fitz.no.diagn_NODES_19dec19.png",,width=11,height=11,units="in",res=1200)

load(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lecanoromycetes_photobiont_growth_form/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave")
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
#nodelabels(pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[,1:8],cex=0.2,piecol=cbPalette[1:8])
#title(main=paste("Bootstrap Topology ",q))

changes<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/lecanoromycetes_photobiont_growth_form/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/node.info.csv",stringsAsFactors=FALSE)

#mark nodes that differ from parent node
changes$change<-FALSE
for(x in 1:nrow(changes)){
	if(any(changes[x,c(7:15)]=="TRUE")){
		changes[x,16]<-"TRUE"
	}
}

#need to add in root
root<-changes[1,]
root$To<-3374
root$ToAge<-root$FromAge
#root state info is OK as is
#but say it changed, so it will plot the root state
root$change<-TRUE
#combine
changes<-rbind(root,changes)

#take TO nodes 3374 and up
changes.sb<-changes[changes$To>3373,]
nrow(changes.sb)
#good 3372 internal nodes

#change states at nodes that don't change
#changes.sb[changes.sb$change %in% "TRUE","To"]-3373

#skip plotting nodes
nodelabels(node=c(changes.sb[changes.sb$change %in% "TRUE","To"]), pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[changes.sb$change %in% "TRUE",],cex=0.3,piecol=cbPalette[1:8])

dev.off()
