library(sciClone)
library(clonevol)
library(fishplot)
OUT='/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/06Sciclone/BM36'
GROUP='BM36'
OUThead=paste(OUT, GROUP, sep = "/")
## read in the data - copy number
cn1 = read.table("/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/03CNV/BM36_L.cnv.txt",sep="\t",stringsAsFactors=F)[,c(1,2,3,5)]
cn2 = read.table("/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/03CNV/BM36_M.cnv.txt",sep="\t",stringsAsFactors=F)[,c(1,2,3,5)]
cn3 = read.table("/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/03CNV/BM36_P.cnv.txt",sep="\t",stringsAsFactors=F)[,c(1,2,3,5)]
head(cn1)
## read in vaf data
v1 = read.table("/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM36_L.vafs",sep="\t",stringsAsFactors=F,header=T)
v2 = read.table("/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM36_M.vafs",sep="\t",stringsAsFactors=F,header=T)
v3 = read.table("/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM36_P.vafs",sep="\t",stringsAsFactors=F,header=T)
samples = c('BM36_L', 'BM36_M', 'BM36_P')
## run sciclone to detect clusters
sc = sciClone(vafs=list(v1,v2,v3),
              copyNumberCalls=list(cn1,cn2,cn3),
              sampleNames=samples,
              minimumDepth=20,
              copyNumberMargins=0.4,
              useSexChrs=FALSE, doClusteringAlongMargins=FALSE) 
sc.plot1d(sc,paste(OUThead,".clusters.1d.pdf",sep=""))
writeClusterTable(sc, paste(OUThead,"clusters.xls",sep=""))
sc.plot2d(sc,paste(OUThead,".clusters.2d.pdf",sep=""))
summary(sc)
head(sc@vafs.merged$P.vaf)
## prepare clonevol input
vafs = data.frame(cluster=sc@vafs.merged$cluster,
                  BM36_L = sc@vafs.merged$BM36_L.vaf,
                  BM36_M = sc@vafs.merged$BM36_M.vaf,
                  BM36_P = sc@vafs.merged$BM36_P.vaf,
                  stringsAsFactors=F
                  )
vafs = vafs[!is.na(vafs$cluster) & vafs$cluster > 0,]
names(vafs)[2:length(samples)+1] = samples
## run clonevol
mycolors = c('#A473AE','#E7A72D',"#ED5401","#0072B2","#009E73", "#EE7777","#FF0000", "#FFCC00",'#777777',"#000000","#0000FF", '#6187C4','#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', "#000001",'#33a02c', '#fb9a99', '#d9d9d9','#999999', '#ff7f00', '#1f78b4', '#fca27e', '#ffffb3', '#fccde5', '#fb8072', '#b3de69', '#f0ecd7', '#e5f5f9')
clone.colors = mycolors[1:length(samples)]
res = infer.clonal.models(
                        variants=vafs, 
                        cluster.col.name="cluster", 
                        vaf.col.names=samples,
                        subclonal.test="bootstrap",
                        subclonal.test.model="non-parametric",
                        cluster.center="median", 
                        num.boots=2000, 
                        founding.cluster=1,
                        clone.colors=clone.colors,
                        min.cluster.vaf=0.01, 
                        sum.p=0.05, 
                        alpha=0.1, 
                        random.seed=334
                        )
res = convert.consensus.tree.clone.to.branch(res, branch.scale='sqrt')
# new clonevol
clone.colors = mycolors[1:7]
plot.clonal.models(res,
                 box.plot = TRUE,
                 fancy.boxplot = TRUE,
                 fancy.variant.boxplot.highlight = 'is.driver',
                 fancy.variant.boxplot.highlight.shape = 15,
                 fancy.variant.boxplot.highlight.fill.color = 'red',
                 fancy.variant.boxplot.highlight.color = 'black',
                 fancy.variant.boxplot.highlight.note.col.name = 'gene',
                 fancy.variant.boxplot.highlight.note.color = 'blue',
                 fancy.variant.boxplot.highlight.note.size = 2,
                 fancy.variant.boxplot.jitter.alpha = 1,
                 fancy.variant.boxplot.jitter.center.color = 'grey50',
                 fancy.variant.boxplot.base_size = 12,
                 fancy.variant.boxplot.plot.margin = 1,
                 fancy.variant.boxplot.vaf.suffix = '.VAF',
                 clone.shape = 'bell',
                 clone.time.step.scale = 0.8,
                 merged.tree.plot = TRUE,
                 tree.node.label.split.character = ',',
                 tree.node.shape = 'circle',
                 tree.node.size = 50,
                 tree.node.text.size = 0.5,
                 merged.tree.node.size.scale = 0.5,
                 merged.tree.node.text.size.scale = 1.3,
                 merged.tree.cell.frac.ci = FALSE,
                 merged.tree.clone.as.branch = TRUE,
                 mtcab.event.sep.char = ',\n',
                 mtcab.branch.text.size = 1,
                 mtcab.branch.width = 0.5,
                 mtcab.node.size = 2.3,
                 mtcab.node.label.size = 1.2,
                 mtcab.node.text.size = 0.8,
                 cell.plot = TRUE,
                 num.cells = 100,
                 cell.border.size = 0.25,
                 cell.border.color = 'black',
                 clone.grouping = 'horizontal',
                 scale.monoclonal.cell.frac = TRUE,
                 show.score = FALSE,
                 cell.frac.ci = TRUE,
                 disable.cell.frac = FALSE,
                 out.dir = './',
                 out.format = 'pdf',
                 overwrite.output = TRUE,
                 width = 12,
                 height = 9,
                 panel.widths = c(3,4,2,4,2)
                )
## create a list of fish objects - one for each model (in this case, there's only one)
f = generateFishplotInputs(results=res)
fishes = createFishPlotObjects(f)
## plot each of these with fishplot
pdf(paste(OUThead,'.fish.pdf',sep="") , width=10, height=5)
for (i in 1:length(fishes)){
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="spline", title.btm="PatientID", cex.title=0.7,
           vlines=seq(1, length(samples)), vlab=samples, pad.left=0.5)
}
 dev.off()