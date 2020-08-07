require(BAMMtools)
require(diversitree)
require(ape)
library(circlize)
library(phytools)
mcmcout <- read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/bamm_analyses/bamm_medusa_med_families_65m/bamm_medusa_med_families_65m_mcmc_out.txt", header=T)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#tree<-read.tree("~/Desktop/other/lecanoromycetes_diversification_2017_analyses/bamm_analyses/bamm_medusa_med_families_65m/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/res.ind.eq.root.fitz.Rsave")
tree<-res.ind.eq.root.fitz$phy
edata <- getEventData(tree, eventdata = "~/Desktop/other/lecanoromycetes_diversification_2017_analyses/bamm_analyses/bamm_medusa_med_families_65m/bamm_medusa_med_families_65m_event_data.txt", burnin=0.1,nsamples=2000)
shift_probs <- summary(edata)
postfile<-mcmcout
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=50, burnin=0.1)

#plots single best shift configuration
msc.set <- maximumShiftCredibility(edata, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)

edata$tip.label->tiplabs
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% tiplabs,]
dat<-dat[match(tiplabs,rownames(dat)),]
dat<-dat[,c(1:11,35,30,21)]

colnames(dat)[12:14]<-c("Form","Photobiont","Cephalodia")
dat$Photobiont[dat$Photobiont %in% "0&2"]<-4
dat$Photobiont[dat$Photobiont %in% "2&3"]<-5

for(x in 12:14){
	dat[,x]<-as.numeric(dat[,x])
}

traits<-c("Form","Photobiont","Cephalodia")

thtr<-tree


png("~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/lecanoromycetes_tips_rates_w_nodestates_nodechanges.png",width=14,height=12.51,units="in",res=600)
#then do plotting


p2p<-diversitree:::plot2.phylo
t <- max(branching.times(thtr))
obj<-p2p(thtr,type="fan",label.offset = t * 1/4, show.tip.label = FALSE)

#get center by vz<-pb(...)
vz<-pb(msc.config,breaksmethod="jenks",pal="GnPu",spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.75,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE)
vz$coords[,1][row.names(vz$coords)==nrow(dat)]
max(vz$coords[,4])

ht<-max(nodeHeights(thtr))
#since bamm does not put root right at center, need to make height greater than root height
ht<-ht*max(vz$coords[,4])
circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1.5,1.5),canvas.ylim=c(-1.5,1.5))
circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 

timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
require(phytools)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<ht,]

for(x in 1:nrow(timey)){
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-timey$End[x])/ht,rou2=(ht-timey$Start[x])/ht,col=timey$RGB[x],border=NA)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-50)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-100)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-150)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-200)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-250)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
}
for(x in 3:7){
	#print(x)
	text(x=(ht-rev(timey$Midpoint[x]))/ht,y=-0.03,labels=timey$Abbrev[x],cex=0.75,col="gray50")
}
circos.clear()
par(new = T)

#get center by vz<-pb(...)
#bammplot<-pb(msc.config,breaksmethod="jenks",pal="terrain",spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.75,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE,legend=FALSE)
#bammplot<-pb(msc.config,breaksmethod="jenks",pal=c("#6B8E23", "#859F36", "#9FB14A", "#B9C25D", "#D3D471", "#EEE685", "#E7D26A", "#E0BF4F", "#DAAB35", "#D3981A", "#CD8500"),spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.75,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE,legend=FALSE)
require(RColorBrewer)
#myterrain.fun<-colorRampPalette(c("olivedrab","khaki3","orange3"))
#myterrain.fun<-colorRampPalette(c("olivedrab","darkgoldenrod1","plum4"))
#myterrain.fun<-colorRampPalette(c("plum4","olivedrab2","orange1"))
myterrain.fun<-colorRampPalette(c("plum4","yellowgreen","orange1"))
myterrain<-myterrain.fun(11)
#bammplot<-pb(msc.config,breaksmethod="jenks",pal=myterrain,spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.75,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE,legend=FALSE)
bammplot<-pb(msc.config,breaksmethod="jenks",pal="black",spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.75,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE,legend=FALSE)
#addBAMMshifts(msc.config,index=1,method="polar",cex=0.75,pch=21,col="slategray3",bg="slategray1")
xy<-obj$xy
theta<-xy$theta[seq_along(thtr$tip.label)]
dt<-diff(sort(theta))[1]/2


#st.cols<-c("white","plum4","white","chartreuse3","white","orange","white","blue","white","brown")
#cols=list("Macrolichen"=c(st.cols[1],st.cols[2]),"Trebouxiophyceae"=c(st.cols[3],st.cols[4]),"Trentepohliales"=c(st.cols[5],st.cols[6]),"Cyanobacteria"=c(st.cols[7],st.cols[8]),"Cephalodia"=c(st.cols[9],st.cols[10]))
#st.cols<-c("#CC79A7","#56B4E9", "gold", "#D55E00","#009E73", "#0072B2", "light green","purple")
st.cols<-c("lemonchiffon1","lemonchiffon4","white","goldenrod2","forestgreen","blue", "light green","purple","darkolivegreen1","lavenderblush4")
cols=list("Form"=c(st.cols[1],st.cols[2]),"Photobiont"=c(st.cols[3],st.cols[4],st.cols[5],st.cols[6],st.cols[7],st.cols[8]),"Cephalodia"=c(st.cols[9],st.cols[10]))

dat.plot<-dat[,c(12:14)]



dat$cladename<-dat$family
orders<-c("Pertusariales","Arctomiales","Trapeliales","Lecideales","Caliciales","Peltigerales","Teloschistales")
families<-c("Stictidaceae","Porinaceae","Graphidaceae","Pilocarpaceae","Psoraceae","Ramalinaceae","Lecanoraceae sl","Parmeliaceae")
subclass<-c("Acarosporomycetidae","Umbilicariomycetidae")
other<-"Cladoniineae"

keepers<-c(orders,families,subclass,other)
for(x in 1:length(orders)){
	dat$cladename[dat$order %in% orders[x]]<-orders[x]
}
for(x in 1:length(families)){
	dat$cladename[dat$family %in% families[x]]<-families[x]
}
for(x in 1:length(subclass)){
	dat$cladename[dat$subclass %in% subclass[x]]<-subclass[x]
}

dat$cladename[dat$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")]<-other

div.fa<-diversitree:::filled.arcs
w = 3/50
for (i in seq_along(cols)){
	idx <- dat.plot[[names(dat.plot)[i]]]
	if (any(idx == 0, na.rm = TRUE)){
		idx <- idx + 1
	}	
	div.fa(theta-dt,theta+dt,rep(1.19+(i*0.02)+(i*w),length(theta)),w,cols[[i]][idx])
}
div.arcs<-diversitree:::arcs
n.taxa <- obj$n.taxa
col.bar<-"black"
n <- obj$n.spp
offset.bar=w * (n + 3)
offset.lab=w * (n + 4)
dy <-1/6
dy <- dy/n * 2 * pi
yy <- obj$xy$theta[seq_len(obj$Ntip)]
#y0 <- tapply(yy - dy, tax2$Tribe, min)
#y1 <- tapply(yy + dy, tax2$Tribe, max)
y0 <- tapply(yy - dy, dat$cladename, min)
y1 <- tapply(yy + dy, dat$cladename, max)
x.bar <- rep(max(obj$xx) + offset.bar, length(y0))
x.lab <- rep(max(obj$xx) + offset.lab, length(y0))
#keepers<-names(table(dat$family))[table(dat$family)>19]
str[!str %in% keepers]<-NA
y0<-y0[names(y0) %in% keepers]
y1<-y1[names(y1) %in% keepers]
x.bar<-x.bar[1:length(keepers)]
x.lab<-x.lab[1:length(keepers)]
div.arcs(y0, y1, x.bar-max(x.bar)+1.51, col = col.bar, lwd = 2)
div.rt<-diversitree:::radial.text
font=3
cex=1
ym <- (y0 + y1)/2
str <- names(y0)
str<-gsub("_"," ",str)
col.lab="black"
div.rt(x.lab-x.lab+1.53, ym, str, col = col.lab, font = font,cex = cex)

require(corHMM)
load(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/my.rate.mat.corr.cleanerer.fitz.no.diagn.Rsave")
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9999CC","#000000")
cbPalette <- c("white", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999", "#D55E00", "#CC79A7","#9999CC","#0072B2")


#plot only nodes w changes
changes<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/node.info.csv",stringsAsFactors=FALSE)
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

nodelabels(node=c(changes.sb[changes.sb$change %in% "TRUE","To"]), pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[changes.sb$change %in% "TRUE",],cex=0.3,piecol=cbPalette[1:8])

#plot ALL nodes
#nodelabels(pie=my.rate.mat.corr.cleanerer.fitz.no.diagn$states[,1:8],cex=0.2,piecol=cbPalette[1:8])

#legend(x=-2.2,y=-1.3,title=NA,c(expression(bold("Growth Form")),"Micro","Micro","Micro","Macro","Micro","Macro","Micro","Macro"),col=c(NA,cbPalette),pch=19,bty="n",cex=0.9,pt.cex=2)
#legend(x=-1.9,y=-1.3,title=NA,c(expression(bold("Primary Photobiont")),"None","Trentepohliales","Trebouxiophyceae","Trebouxiophyceae","Cyanobacteria","Cyanobacteria","Trebouxiophyceae","Trebouxiophyceae"),col=NA,pch=19,bty="n",cex=0.9,pt.cex=2)
#legend(x=-0.9,y=-1.3,title=NA,c(expression(bold("Cephalodia")),"None","None","None","None","None","None","Present","Present"),col=NA,pch=19,bty="n",cex=0.9,pt.cex=2)

legend(x=0.8,y=2,title=NA,c(expression(bold("Growth Form")),"Micro","Micro","Micro","Macro","Micro","Macro","Micro","Macro"),col=c(NA,cbPalette),pch=19,bty="n",cex=0.9,pt.cex=2)
legend(x=1.2,y=2,title=NA,c(expression(bold("Primary Photobiont")),"None","Trentepohliales","Trebouxiophyceae","Trebouxiophyceae","Cyanobacteria","Cyanobacteria","Trebouxiophyceae","Trebouxiophyceae"),col=NA,pch=19,bty="n",cex=0.9,pt.cex=2)
legend(x=1.8,y=2,title=NA,c(expression(bold("Cephalodia")),"None","None","None","None","None","None","Present","Present"),col=NA,pch=19,bty="n",cex=0.9,pt.cex=2)

segments(x0=0.8,y0=1.915,x1=2.2,y1=1.915)
#segments(x0=1.17,y0=-1.51,x1=1.9,y1=-1.51)
text(x=1.5,y=1.97,expression(bold("Node States")),cex=1.5)




#legend(x=1.2,y=-1.2,title=NA,c(expression(bold("Growth Form (inner)")),"Microlichen","Macrolichen", NA,expression(bold("Primary Photobiont (outer)")),"None","Trentepohliales","Trebouxiophyceae","Cyanobacteria"),col=c(NA,"#CC79A7","#56B4E9",NA,NA,"gold", "#D55E00","#009E73", "#0072B2"),pch=19,bty="n",cex=1,pt.cex=2)
legend(x=1.35,y=-1,title=NA,c(expression(bold("Growth Form (inner)")),"Micro","Macro", NA,expression(bold("Primary Photobiont (middle)")),"None","Trentepohliales","Trebouxiophyceae","Cyanobacteria","None & Trebouxiophyceae","Trebouxiophyceae & Cyanobacteria",NA,expression(bold("Cephalodia (outer)")),"None","Present"),col=c(NA,"lemonchiffon1","lemonchiffon4",NA,NA,"white","goldenrod2","forestgreen","blue","light green","purple",NA,NA,"darkolivegreen1","lavenderblush4"),pch=19,bty="n",cex=0.9,pt.cex=2)
segments(x0=1.35,y0=-1.08,x1=2.35,y1=-1.08)
text(x=1.85,y=-1.04,expression(bold("Tip States")),cex=1.5)
#text(x=-1.62,y=-0.295,expression(bold("Tip States")),cex=1.5)
##text(x=1.7,y=-1.00,expression(bold("Tip States")),cex=1.4)
##segments(x0=1.17,y0=-1.04,x1=2.3,y1=-1.04)
#segments(x0=1.17,y0=-1.38,x1=2.3,y1=-1.38)
#segments(x0=1.17,y0=-1.69,x1=2.3,y1=-1.69)
#text(x=1.45,y=-0.92,"Diet (inner)",cex=1,pos=1)
#text(x=1.55,y=-1.3,"Foraging (middle)",cex=1)
#text(x=1.5,y=-1.55,"Nesting (outer)",cex=1)
#addBAMMlegend(bammplot,direction="horizontal",location="bottomleft")
#text(x=-1.45,y=-1.6,expression(bold("Speciation Rate")),cex=1.4)

dev.off()




















modified.addBAMMshifts<-function (ephy, index = 1, method = "phylogram", cex = 1, pch = 21, 
    col = 1, bg = 2, msp = NULL, shiftnodes = NULL, par.reset = TRUE) 
{
      as.phylo.bammdata<-BAMMtools:::as.phylo.bammdata
	if (!"bammdata" %in% class(ephy)) 
        stop("Object ephy must be of class bammdata")
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (par.reset) {
        op <- par(no.readonly = TRUE)
        par(lastPP$pp)
    }
    if (length(ephy$eventData) == 1) {
        index <- 1
    }
    times <- ephy$begin[ephy$edge[,1] %in% shiftnodes][c(TRUE,FALSE)]
    if (!is.null(msp)) {
        cex <- 0.75 + 5 * msp$edge.length[msp$edge[, 2] %in% 
            shiftnodes]
    }
    if (method == "phylogram") {
        XX <- times
        YY <- lastPP$yy[shiftnodes]
    }
    else if (method == "polar") {
        rb <- lastPP$rb
        XX <- (rb + times/max(branching.times(as.phylo.bammdata(ephy)))) * 
            cos(lastPP$theta[shiftnodes])
        YY <- (rb + times/max(branching.times(as.phylo.bammdata(ephy)))) * 
            sin(lastPP$theta[shiftnodes])
    }
    points(XX, YY, pch = pch, cex = cex, col = col, bg = bg)
    if (par.reset) {
        par(op)
    }
}


#took out ofs=0 and moved to options in plot.bammdata in BAMMtools 2.1.4
pb<-function (x, tau = 0.01, method = "phylogram", ofs=0, xlim = NULL, ylim = NULL, 
    vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, legend = FALSE, 
    spex = "s", lwd = 1, cex = 1, pal = "RdYlBu", mask = integer(0), 
    mask.color = gray(0.5), colorbreaks = NULL, logcolor = FALSE, 
    breaksmethod = "linear", color.interval = NULL, JenksSubset = 20000, 
    par.reset = FALSE, direction = "rightwards", ...) 
{
    as.phylo.bammdata<-BAMMtools:::as.phylo.bammdata
    colorMap<-BAMMtools:::colorMap
    setPolarTreeCoords<-BAMMtools:::setPolarTreeCoords
    mkdtsegsPolar<-BAMMtools::: mkdtsegsPolar
    div.arcs<-diversitree:::arcs
    arc<-BAMMtools:::arc
	if ("bammdata" %in% class(x)) {
        if (attributes(x)$order != "cladewise") {
            stop("Function requires tree in 'cladewise' order")
        }
        phy <- as.phylo.bammdata(x)
    }
    else stop("Object ephy must be of class bammdata")
    if (!spex %in% c("s", "e", "netdiv")) {
        stop("spex must be 's', 'e' or 'netdiv'.")
    }
    if (length(pal) == 1 && !pal %in% names(get("palettes", envir = .colorEnv)) && 
        pal != "temperature" && pal != "terrain") 
        pal <- rep(pal, 3)
    else if (length(pal) == 2) 
        pal <- c(pal, pal[2])
    if (breaksmethod == "linear" & !is.null(color.interval)) {
        if (length(color.interval) != 2) {
            stop("color.interval must be a vector of 2 numeric values.")
        }
    }
    if (!is.binary.tree(phy)) {
        stop("Function requires fully bifurcating tree")
    }
    if (any(phy$edge.length == 0)) {
        warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero")
    }
    if (!("dtrates" %in% names(x))) {
        x <- dtRates(x, tau)
    }
    if (is.null(colorbreaks)) {
        colorbreaks <- assignColorBreaks(x$dtrates$rates, 64, 
            spex, logcolor, breaksmethod, JenksSubset)
    }
    if (x$type == "trait") {
        colorobj <- colorMap(x$dtrates$rates, pal, colorbreaks, 
            logcolor, color.interval)
    }
    else if (x$type == "diversification") {
        if (tolower(spex) == "s") {
            colorobj <- colorMap(x$dtrates$rates[[1]], pal, colorbreaks, 
                logcolor, color.interval)
        }
        else if (tolower(spex) == "e") {
            colorobj <- colorMap(x$dtrates$rates[[2]], pal, colorbreaks, 
                logcolor, color.interval)
        }
        else if (tolower(spex) == "netdiv") {
            colorobj <- colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], 
                pal, colorbreaks, logcolor, color.interval)
        }
    }
    else {
        stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'")
    }
    edge.color <- colorobj$cols
    tH <- max(x$end)
    phy$begin <- x$begin
    phy$end <- x$end
    tau <- x$dtrates$tau
    if (method == "polar") {
        ret <- setPolarTreeCoords(phy, vtheta, rbf)
        rb <- tH * rbf
        p <- mkdtsegsPolar(ret$segs[-1, ], tau, x$edge)
    }
    else if (method == "phylogram") {
        ret <- setPhyloTreeCoords(phy)
        p <- mkdtsegsPhylo(ret$segs[-1, ], tau, x$edge)
    }
    else {
        stop("Unimplemented method")
    }
    x0 <- c(ret$segs[1, 1], p[, 1])
    x1 <- c(ret$segs[1, 3], p[, 2])
    y0 <- c(ret$segs[1, 2], p[, 3])
    y1 <- c(ret$segs[1, 4], p[, 4])
    offset <- table(p[, 5])[as.character(unique(p[, 5]))]
    if (length(mask)) {
        edge.color[p[, 5] %in% mask] <- mask.color
    }
    arc.color <- c(edge.color[1], edge.color[match(unique(p[, 
        5]), p[, 5]) + offset])
    edge.color <- c(edge.color[1], edge.color)
    if (show) {
        op <- par(no.readonly = TRUE)
        if (length(list(...))) {
            par(...)
        }
        if (legend) {
            par(mar = c(5, 4, 4, 5))
        }
        plot.new()
        #ofs <- 0     #####THIS
        if (labels) {
            if (method == "phylogram") 
                ofs <- max(nchar(phy$tip.label) * 0.03 * cex * 
                  tH)
            else ofs <- max(nchar(phy$tip.label) * 0.03 * cex)
        }
        if (method == "polar") {
            plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                ofs), asp = 1)
            segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, 
                lend = 2)
            arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + 
                phy$end/tH), border = arc.color, lwd = lwd)
            if (labels) {
                for (k in 1:length(phy$tip.label)) {
                  text(ret$segs[-1, ][phy$edge[, 2] == k, 3], 
                    ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k], 
                    cex = cex, srt = (180/pi) * ret$arcs[-1, 
                      ][phy$edge[, 2] == k, 1], adj = c(0, NA))
                }
            }
        }
        if (method == "phylogram") {
            direction <- match.arg(direction, c("rightwards", 
                "leftwards", "downwards", "upwards"))
            if (direction == "rightwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), 0)
                arcs <- redirect(ret$arcs, 0)
                bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
                arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
                ret$segs[-1, c(1, 3)] <- tH * ret$segs[-1, c(1, 
                  3)]
            }
            else if (direction == "leftwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), pi)
                bars[, c(2, 4)] <- abs(bars[, c(2, 4)])
                arcs <- redirect(ret$arcs, pi)
                arcs[, c(2, 4)] <- abs(arcs[, c(2, 4)])
                bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
                arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
                ret$segs[-1, c(1, 3)] <- -tH * ret$segs[-1, c(1, 
                  3)]
            }
            else if (direction == "downwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), -pi/2)
                arcs <- redirect(ret$arcs, -pi/2)
                bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
                arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
                ret$segs <- redirect(ret$segs, -pi/2)
                ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2, 4)]
            }
            else if (direction == "upwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), pi/2)
                bars[, c(1, 3)] <- abs(bars[, c(1, 3)])
                arcs <- redirect(ret$arcs, pi/2)
                arcs[, c(1, 3)] <- abs(arcs[, c(1, 3)])
                bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
                arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
                ret$segs <- redirect(ret$segs, pi/2)
                ret$segs[, c(1, 3)] <- abs(ret$segs[, c(1, 3)])
                ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2, 4)]
            }
            if (is.null(xlim) && direction == "rightwards") 
                xlim <- c(0, tH + ofs)
            if (is.null(xlim) && direction == "leftwards") 
                xlim <- c(-(tH + ofs), 0)
            if (is.null(ylim) && (direction == "rightwards" || 
                direction == "leftwards")) 
                ylim <- c(0, phy$Nnode)
            if (is.null(xlim) && (direction == "upwards" || direction == 
                "downwards")) 
                xlim <- c(0, phy$Nnode)
            if (is.null(ylim) && direction == "upwards") 
                ylim <- c(0, tH + ofs)
            if (is.null(ylim) && direction == "downwards") 
                ylim <- c(-(tH + ofs), 0)
            plot.window(xlim = xlim, ylim = ylim)
            segments(bars[-1, 1], bars[-1, 2], bars[-1, 3], bars[-1, 
                4], col = edge.color[-1], lwd = lwd, lend = 2)
            isTip <- phy$edge[, 2] <= phy$Nnode + 1
            isTip <- c(FALSE, isTip)
            segments(arcs[!isTip, 1], arcs[!isTip, 2], arcs[!isTip, 
                3], arcs[!isTip, 4], col = arc.color[!isTip], 
                lwd = lwd, lend = 2)
            if (labels) {
                if (direction == "rightwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 4, offset = 0.25)
                else if (direction == "leftwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 2, offset = 0.25)
                else if (direction == "upwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 4, srt = 90, offset = 0)
                else if (direction == "downwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 2, srt = 90, offset = 0)
            }
        }
    }
    index <- order(as.numeric(rownames(ret$segs)))
    if (show) {
        if (method == "phylogram") {
            assign("last_plot.phylo", list(type = "phylogram", 
                direction = direction, Ntip = phy$Nnode + 1, 
                Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                  3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)), 
                envir = .PlotPhyloEnv)
        }
        else if (method == "polar") {
            assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 
                1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                3], yy = ret$segs[index, 4], theta = ret$segs[index, 
                5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv)
        }
        if (legend) {
            addBAMMlegend(x = list(coords = ret$segs[-1, ], colorbreaks = colorbreaks, 
                palette = colorobj$colpalette, colordens = colorobj$colsdensity), 
                location = "right")
        }
    }
    if (par.reset) {
        par(op)
    }
    invisible(list(coords = ret$segs[-1, ], colorbreaks = colorbreaks, 
        palette = colorobj$colpalette, colordens = colorobj$colsdensity))
}










#from https://github.com/macroevolution/bammtools/blob/587117a2a79159cc9809e346d7dedf58a4aa2add/BAMMtools/R/bammcolors.R
palettes <- list(

BrBG =     rev(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3",

                 "#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e",

                 "#003c30")),

PiYG =     rev(c("#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef",

                 "#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221",

                "#276419")),

PRGn =     rev(c("#40004b","#762a83","#9970ab","#c2a5cf","#e7d4e8",

                 "#f7f7f7","#d9f0d3","#a6dba0","#5aae61","#1b7837",

                "#00441b")),

PuOr =     rev(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6",

                 "#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788",

                 "#2d004b")),

RdBu =     rev(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",

                 "#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac",

                 "#053061")),

RdYlBu =   rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",

                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",

                 "#313695")),

BuOr =     c("#002bff","#1a66ff","#3399ff","#66CCff","#99eeff",

            "#ccffff","#ffffcc","#ffee99","#ffee66","#ff9933",

            "#ff661a","#ff2b00"),

BuOrRd =   c("#085aff","#3377ff","#5991ff","#8cb2ff","#bfd4FF",

             "#e6eeff","#f7faff","#ffffcc","#ffff99","#ffff00",

             "#ffcc00","#ff9900","#ff6600","#ff0000"),

DkRdBu =   c("#2a0bd9","#264eff","#40a1ff","#73daff","#abf8ff",

             "#e0ffff","#ffffbf","#ffe099","#ffad73","#f76e5e",

             "#d92632","#a60021"),

BuDkOr =   c("#1f8f99","#52c4cc","#99faff","#b2fcff","#ccfeff",

             "#e6ffff","#ffe6cc","#ffca99","#ffad66","#ff8f33",

             "#cc5800","#994000"),

GnPu =     c("#005100","#008600","#00bc00","#00f100","#51ff51",

             "#86ff86","#bcffbc","#ffffff","#fff1ff","#ffbcff",

             "#ff86ff","#ff51ff","#f100f1","#bc00bc","#860086",

             "#510051"),

RdYlGn =   rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee08b",

                 "#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850",

                 "#006837")),

Spectral = rev(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b",

             "#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd",

             "#5e4fa2")),

grayscale = c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000"),

revgray = rev(c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000")),

greyscale = c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000"),

revgrey = rev(c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000"))

);

.colorEnv <- new.env();

assign("palettes", palettes, env = .colorEnv);


