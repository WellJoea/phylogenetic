#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
***********************************************************
* Author : Zhou Wei                                       *
* Date   : Fri Jul 27 12:36:42 2018                       *
* E-mail : welljoea@gmail.com                             *
* You are using the program scripted by Zhou Wei.         *
* If you find some bugs, please                           *
* Please let me know and acknowledge in your publication. *
* Thank you!                                              *
* Best wishes!                                            *
***********************************************************
'''
import argparse
import re
import numpy as np
import pandas as pd
import time
import glob
from multiprocessing import Process, Pool
from functools import reduce
import multiprocessing.pool
import traceback
import sys
import os
import gzip

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)
# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.

class NoDaemonProcessPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

class Openf():
    def __init__(self, file, *array,compression='infer',comment=None,sep="\t",header=0):
        self.file  = file.strip()
        self.array = array
        self.compr = compression
        self.comme = comment
        self.sep   = sep
        self.header= header
    def opena(self):
        with open(self.file, 'r') as inp:
            outa = [i.strip('\n').split("\t") for i in inp.readlines() if not re.match(r'##',i)]
        return outa
    def openb(self):
        if re.match(r'.*xlsx$|.*xlss$', self.file, re.I):
            outp = pd.read_excel(self.file, header=self.header, encoding='utf-8', sheet_name=0, comment=self.comme, compression=self.compr).fillna(np.nan)
        else:
            outp = pd.read_table(self.file, header=self.header, encoding='utf-8', sep=self.sep, comment=self.comme,compression=self.compr).fillna(np.nan)
        return(pd.DataFrame(outp))
    def openw(self):
        os.system('mkdir -p '+ os.path.dirname(self.file))
        f= open(self.file, 'w',encoding='utf-8')
        f.write(self.array[0])
        f.close()
    def openv(self):
        os.system('mkdir -p '+ os.path.dirname(self.file))
        df = self.array[0]
        df.to_csv(self.file, sep=self.sep, index=False,header=True, line_terminator='\n')

class WorkShell():
    def __init__(self, *lis, **dic):
        self.lis = lis
        self.dic = dic
    def resub(self,strs):
        for i,k in self.dic.items():
            strs = re.sub(r'''@%s@'''%i, k, strs)
        return strs
    def annovar(self):
        worksh = r'''
FILE=@FILE@
Annovfi=@Annovfi@
mkdir -p `dirname $Annovfi`
perl /lustre/rde/user/zhouw/00software/annovar/table_annovar.pl  \
    $FILE \
    /lustre/rde/user/zhouw/00software/annovar/humandb/  \
    --outfile $Annovfi \
    --thread 10 \
    --buildver hg19 --otherinfo --remove --nastring . --thread 10 --vcfinput \
    --operation g,f,f,f,f,f,f,f,r \
    --argument ,,,,,,,, \
    --protocol refGene,1000g2015aug_all,exac03,gnomad_exome,snp150,cosmic80,HGMD,HGMD_mut,OMIM
'''.strip()
        return self.resub(worksh)
    def LicheeTex(self):
        worksh = r'''
INput=@Annovar@
OUTput=`realpath @OUTdir@`
Script=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/script
CLONE=$OUTput/03LICHEE
Textplot=$OUTput/04Texplot
mkdir -p $OUput $CLONE

perl $Script/change_format_mutect.pl $INput $CLONE

while read i
do
  java -jar /lustre/rde/user/zhouw/00software/lichee/LICHeE/release/lichee.jar \
      -build -maxVAFAbsent 0.02 \
      -minVAFPresent 0.08 \
      -n 0 \
      -e 0.10 \
      -minPrivateClusterSize 2 \
      -maxClusterDist 0.2 \
      -maxVAFValid  0.6 \
      -showTree 1 -color \
      -dot -i $i -o $i
  /usr/bin/dot -Tpdf ${i}.dot -O $CLONE
done <$CLONE/all_mut_evol.list.txt
python $Script/Lichee.site.plot.py -i $CLONE/all_mut_evol.list.txt -o $Textplot
'''.strip()
        return self.resub(worksh)
    def sciclone(self):
        worksh = r'''
library(sciClone)
library(clonevol)
library(fishplot)
OUT='@OUDDir@'
GROUP='@GROUP@'
OUThead=paste(OUT, GROUP, sep = "/")
## read in the data - copy number
@CNVRead@
head(cn1)
## read in vaf data
@SNVRead@
samples = c@SAMPLEs@
## run sciclone to detect clusters
sc = sciClone(vafs=list(@SNVList@),
              copyNumberCalls=list(@CNVList@),
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
                  @SNVexpr@
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
'''.strip()
        return self.resub(worksh)

    def pyclone(self):
        worksh = r'''
OUTput=@OUTdir@
mkdir -p $OUTput
/lustre/rde/user/zhouw/00software/anaconda2/bin/PyClone run_analysis_pipeline \
    --num_iters 10000 \
    --burnin 1000 \
    --plot_file_format pdf \
    --in_files \
        @TSVs@ \
   --working_dir $OUTput
'''.strip()
        return self.resub(worksh)

class Pipline():
    def __init__(self, arg, *array):
        self.arg = arg
        self.array = array
    def Get(self,efile):
        if 'annovar' in self.arg.analysis:
            efile = efile.reset_index(drop=True)
            (FILE,SAMPLE,GROUP)  = efile.loc[0,['FILE','SAMPLE','GROUP']]
            Annovsh = self.arg.outdir + '/00Worksh/S01_%s.annovar.sh'%(GROUP + '-' + SAMPLE)
            Annovfi = self.arg.outdir + '/01ANNOVAR/' + GROUP + '-' + SAMPLE + '-annovar'
            Annov = WorkShell( FILE = FILE, Annovfi =Annovfi).annovar()
            Openf(Annovsh ,(Annov)).openw()
            if not os.path.exists(Annovfi + '.hg19_multianno.txt'):
                os.system('qsub -cwd %s'%Annovsh)
            return Annovfi + '.hg19_multianno.txt'

        if 'sciclone' in self.arg.analysis:
            Group,Groupdf = efile
            Sampall = Groupdf['SAMPLE'].tolist()
            SAMPLElist = list(enumerate(Sampall, start=1))
            SNVpath = self.arg.outdir + '/02SNVCNV/02SNV/'
            CNVpath = self.arg.outdir + '/02SNVCNV/03CNV/'
            OUTpath = self.arg.outdir + '/06Sciclone/' + Group
            OUTrscr = OUTpath + '/' + Group + '.sci.R'
            OUTsh   = self.arg.outdir + '/00Worksh/S06_%s.sciclone.sh'%(Group)
            cnvread = [r'''cn%s = read.table("%s.cnv.txt",sep="\\t",stringsAsFactors=F)[,c(1,2,3,5)]'''%(i[0], CNVpath + i[1]) for i in SAMPLElist ]
            cnvread = r'\n'.join(cnvread)
            snvread = [r'''v%s = read.table("%s.vafs",sep="\\t",stringsAsFactors=F,header=T)'''%(i[0], SNVpath + i[1]) for i in SAMPLElist ]
            snvread = r'\n'.join(snvread)
            cnvlist = [r'''cn%s'''%(i[0]) for i in SAMPLElist ]
            cnvlist = r','.join(cnvlist)
            snvlist = [r'''v%s'''%(i[0]) for i in SAMPLElist ]
            snvlist = r','.join(snvlist)
            snvexpr = [r'''%s = sc@vafs.merged$%s.vaf,'''%(i[1],i[1]) for i in SAMPLElist ]
            snvexpr = (r'''\n''' + ' '*18).join(snvexpr)
            samples = str(tuple(Sampall))
            SciR    = WorkShell(OUDDir=OUTpath, GROUP=Group, CNVRead=cnvread, SNVRead=snvread,SAMPLEs=samples,SNVList=snvlist,CNVList=cnvlist, SNVexpr=snvexpr ).sciclone()
            Openf(OUTrscr ,(SciR)).openw()
            Openf(OUTsh ,('Rscript %s'%OUTrscr)).openw()
            if self.arg.do:
                os.system('nohup sh %s &'%OUTsh)

        if 'pyclone' in self.arg.analysis:
            GROUP,Groupdf = efile
            SAMPLElist = Groupdf['SAMPLE'].tolist()
            TSVS = r''' \\\n        '''.join([ self.arg.outdir + '/02SNVCNV/02SNV/%s.tsv'%i for i in SAMPLElist])
            OUTsh= self.arg.outdir + '/00Worksh/S05_%s.pyclone.sh'%(GROUP)
            OUT  = self.arg.outdir + '/05Pyclone/%s/'%(GROUP)
            Pysh = WorkShell( TSVs = TSVS, OUTdir =OUT).pyclone()
            Openf(OUTsh ,(Pysh)).openw()
            if self.arg.do:
                os.system('nohup sh %s &'%OUTsh)
                #os.system('qsub -cwd %s'%OUTsh)

    def SUB(self):
        self.arg.outdir = os.path.realpath(self.arg.outdir)
        files = Openf(self.arg.input,comment="#").openb()
        colum = [i.upper() for i in files.columns]
        files.columns = colum
        self.arg.analysis = re.split(r'[,;_]',self.arg.analysis)

        if 'annovar' in self.arg.analysis:
            print('''doing annovar file...''')
            file_split = np.array_split(files.copy(), len(files))
            pool  = NoDaemonProcessPool(processes=self.arg.pool)
            result= pool.map(self.Get, file_split)
            pool.close()
            pool.join()
            Openf(self.arg.outdir + '/01ANNOVAR/All.annovar.files.txt', ('\n'.join(result))).openw()

        if 'Lichee' in self.arg.analysis:
            print('''doing Lichee...''')
            LicheeTexsh = WorkShell(Annovar=self.arg.outdir + '/01ANNOVAR/All.annovar.files.txt',OUTdir=self.arg.outdir).LicheeTex()
            Openf(self.arg.outdir + '/00Worksh/S02get.lichee_tex.sh' ,(LicheeTexsh)).openw()
            if self.arg.do:
                print('''doing licheee and Texplot...''')
                #os.system('qsub -cwd %s'%(self.arg.outdir + '/00Worksh/S02get.lichee_tex.sh'))
                os.system('nohup sh  %s & '%(self.arg.outdir + '/00Worksh/S02get.lichee_tex.sh'))

        if ('sciclone' in self.arg.analysis) or ('pyclone' in self.arg.analysis):
            print('''doing sciclone or/and pyclone...''')
            file_split = files.copy().groupby(['GROUP'])
            pool  = NoDaemonProcessPool(processes=self.arg.pool)
            pool.map(self.Get, file_split)
            pool.close()
            pool.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="\nClone phylogeny analysis :\n")
    parser.add_argument("-i", "--input",     type=str,  required=True, help="the input file including Id,snv.vcf.gz and bam columns")
    parser.add_argument("-o", "--outdir",    type=str,  default=os.getcwd(), help="output file dir, default=current dir")
    parser.add_argument("-p", "--pool",      type=int,  default=23,  help="the CPU numbers that can be used")
    parser.add_argument("-a", "--analysis",  type=str,  default="annovar,Lichee,sciclone,pyclone", help="the analysis type, default=Lichee,choices=annovar,Lichee,sciclone,pyclone")
    parser.add_argument("-d", "--do",        action="store_true",default=False, help='do the workshell')

    global Arags
    args  = parser.parse_args()
    Arags = args.__dict__
    print('start time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print("***The global argument you have set as follows:*****")
    i = 0
    for arg in vars(args):
        i +=1
        print ('**' + str(i).zfill(2) + '|' + arg + ': ' + str(getattr(args, arg)))
    print("****************************************************")

    try:
        Pipline(args).SUB()
        print('end time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        sys.exit(0)
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)
