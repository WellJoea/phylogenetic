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
import time
import glob
import pandas as pd
from multiprocessing import Process, Pool
import multiprocessing.pool
import traceback
import sys
import os

class Openf():
    def __init__(self, file, *array):
        self.file  = file
        self.array = array
    def opena(self):
        with open(self.file, 'r', encoding='utf-8') as inp:
            outa = [i.strip('\n') for i in inp.readlines() if not re.match(r'#|^$',i)]
            if self.array : outa[0] += self.array[0]
        return outa
    def openw(self):
        f= open(self.file, 'w',encoding='utf-8')
        f.write(self.array[0])
        f.close()

class WorkShell():
    def __init__(self, arg, *array):
        self.arg = arg
        self.array = array
    def texplot(self,Plotlist):
        tex = r'''
%-*- coding: UTF-8 -*-
%%%%%%%%%% package %%%%%%%%%% 
\documentclass[UTF8,twoside]{ctexrep}
\usepackage{blkarray}
\usepackage{fancyvrb}
\usepackage{ltxtable}
\usepackage{tabu}
\usepackage{multirow,makecell,multicol}
\usepackage[cmyk,table]{xcolor}
\usepackage{colortbl}
\usepackage{pagecolor,afterpage}
\usepackage{makeidx}
\usepackage{threeparttable}
\usepackage{array}
%\usepackage[a4paper,total={200mm,380mm},left=1cm,right=1cm,top=1cm,bottom=1cm]{geometry}
%\usepackage[paperheight=9in,paperwidth=9in,left=0.3cm,right=1cm,top=3cm,bottom=0.0cm]{geometry}
\usepackage[paperheight=10in,paperwidth=9in,left=0.3cm,right=1cm,top=3cm,bottom=0.0cm]{geometry}
\usepackage{fancyhdr}
\usepackage{graphicx}
\usepackage{float}
\usepackage{balance}
\usepackage{tabu,tikz}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage[justification=RaggedRight,labelformat=empty]{caption}
\usepackage{tocloft}
\usepackage{lastpage}
\usepackage{dashrule}
\usepackage{xhfill}
\usepackage{paralist} 
\usepackage{wallpaper}
\usepackage{ifthen,changepage}
\usepackage{pdfpages}
\usepackage{ragged2e}
\usepackage{pifont}
\usepackage{adjustbox}
\usepackage{pstricks}
\usepackage{hhline}
\usepackage{titletoc}
\usepackage{rotating}
\usepackage{titlesec}
\usepackage{stackengine}
\usepackage{forest}

%\usepackage[colorlinks,linkcolor=black,anchorcolor=blue,citecolor=green]{hyperref} 
\usepackage{hyperref}
%%%%%%%%%% table format command  %%%%%%%%%
\renewcommand{\raggedright}{\leftskip=0pt \rightskip=0pt plus 0cm}
%\justifying\let\raggedright\justifying
\newcolumntype{C}[1]{>{\centering\arraybackslash}m{#1}}
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}m{#1}}
\newcolumntype{Y}{>{\centering\arraybackslash}X}
\newcolumntype{V}{!{\color{white}\vline width 1.8pt}}
\newcolumntype{Y}{!{\color{white}\vline width 3pt}}
\newcolumntype{?}{!{\vrule width 1.5pt}}
\setlength{\aboverulesep}{0.0pt}
\setlength{\belowrulesep}{0.0pt}

%\tikzset{dot/.style={circle,draw,inner sep=1.2,fill=black},}

\begin{document}
\setlength{\aboverulesep}{0pt}
\setlength{\belowrulesep}{0pt}
\setlength\tabcolsep{0 pt}

\definecolor{greyi}{RGB}{230,230,230}
\definecolor{grayi}{RGB}{220,220,220}

\begin{center}
%\begin{flushleft}
%\begin{sidewaystable}
  \color{black}
  \fontsize{14}{14}\selectfont
  \renewcommand\arraystretch{0.4}
\begin{adjustbox}{angle=90}

%[CONTENT]%

\end{tabularx}

\end{adjustbox}
%\end{sidewaystable}
%\end{flushleft}
\end{center}
\end{document}
'''
        tex = re.sub(r'%\[CONTENT\]%',Plotlist, tex)
        return tex

    def HISTO(self,DataPd):
        NumVal = DataPd.columns.tolist().index('Values')
        DataPda = DataPd.iloc[:,1:NumVal+1].copy()
        for i in range(DataPd.shape[0]):
            if i==0 :
                DataPda.iloc[i] = DataPd.iloc[i].replace(['0','1'],['\\addstackgap[20pt]{\\cellcolor{greyi}}','\\addstackgap[20pt]{\\cellcolor{' + DataPd.loc[i,'Colors'] + '}}' ])
            else:
                DataPda.iloc[i] = DataPd.iloc[i].replace(['0','1'],['\\cellcolor{greyi}','\\cellcolor{' + DataPd.loc[i,'Colors'] + '}' ])
            DataPda.loc[i,'Values'] = '\\cellcolor{' + DataPd.loc[i,'Colors'] + '}'

        Tabularx = '@{}V'.join(['C{' + str(i) + '\\linewidth}' for i in DataPd['Values']] )
        Tabularx = '\\begin{tabularx}{\\linewidth}{' + Tabularx + '@{}Y>{\\color{black}\\bfseries}b{0.1\\linewidth}@{}}'
        Plotlist = []
        for i in DataPda.columns:
            if i != 'Values':
                #Plotlist.append(' & '.join(DataPda[i]) + ' & \\rotatebox[origin=l]{200}{\\multirow{1}{*}{' + i + '}} \\\\ \n\\specialrule{0em}{2.6pt}{2.6pt}')
                Plotlist.append(' & '.join(DataPda[i]) + ' & {\\multirow{1}{*}{' + i + '}} \\\\ \n\\specialrule{0em}{2.6pt}{2.6pt}')
            else :
                Plotlist.append('\\specialrule{0em}{3.0pt}{3.0pt}\n' + ' & '.join(DataPda[i]) + ' &  \\\\ \n\\specialrule{0em}{4.5pt}{4.5pt}\n\\addstackgap[20pt]{\\cellcolor{white} } %s \\\\'%('&'*(DataPd.shape[0]-1)))
        Plotlist = '\n'.join([Tabularx] + Plotlist)
        return Plotlist

    def Change(self, DataPd):
        b = DataPd[['o','t','Length']].copy()
        b.sort_values(by=['o','t','Length'], ascending=[True, True, False], axis=0, inplace=True)
        b= b.astype('str')
        st = [['0']]
        sd = ['0']
        for i in range(len(set(b['o']))):
            tv = b[b['o'].isin(st[-1])]['t'].tolist()
            ta = list(set(tv) & set(b['o']))
            if len(ta)>0:
                st.append(ta)
                sd += ta
        strall ='[0]'
        for i in sd:
            i = str(i) 
            endi = b[b['o']==i]['t']
            endp = '\n' +'\n'.join(['[' + i + ']' for i in endi]) + '\n'
            strall = re.sub(r'(?P<HE>[\n\[\]])(?P<STA>%s)(?P<TA>[\n\[\]])'%i, r'\g<HE>\g<STA>%s\g<TA>'%endp, strall)
        #print(r'%r'%strall)
        sign = set(b['t']) - set(b['o'])
        for k, v in DataPd.iterrows():
            label  = r''',label={right: \\\\bf{\\\\makecell[l]{%s}}}'''%(r''' \\\\\\\\ '''.join(v['Names'])) if (v['t'] in sign) else ''
            strty  = r''' ,name=%s,l*=%s, edge={%s,line width=0.3cm}%s'''%(v['t'],v['Length'],v['Colors'],label)
            strall = re.sub(r'''(?P<HE>[\n\[\]])(?P<STA>%s)(?P<TA>[\n\[\]])''' %v['t'], r'''\g<HE>%s\g<TA>'''%strty, strall)
        strall = re.sub(r'\[0\n', r'''[ ,name=root, label={left: \\\\bf{Nomal}}\n''', strall)
        #print(strall)
        return strall

    def ETree(self,DataPd):
        CONTENT = self.Change(DataPd)
        LENT    = len(DataPd) + 1
        Forest  = r'''
\multicolumn{%[LENT]%}{@{}l}{
\begin{forest}
    for tree={l sep=1cm,grow=0,circle,outer sep=0pt,inner sep=4pt,line width=3pt,fill=grayi,draw=black,calign=fixed edge angles,calign angle=30},
%[CONTENT]%
\end{forest}
} \\
'''
        Forest = re.sub(r'%\[LENT\]%' , str(LENT), Forest)
        Forest = re.sub(r'%\[CONTENT\]%', CONTENT, Forest)
        return Forest

    def Get(self,file):
        print(file)
        Dottxt = Openf(file).opena()
        SampNm = os.path.basename(file).split('_')[0]

        colors = ['red','blue','green','yellow','magenta','cyan','purple','black','orange','violet','pink','brown','gray','lime','teal']
        try:
            Nodes  = Dottxt[1:Dottxt.index('****Tree 0****')]
            print(file)
            Name   = [ re.sub(r'^.*Sample lineage decomposition: ','',i ) for i in Dottxt if re.search(r'Sample lineage decomposition',i)]
            Name   = [ re.sub(r'_',r'\_',i ) for i in Name ]
            TREE   = [ re.sub(r'\s','',i) for i in Dottxt if re.search(r'[0-9]+.*\->.*[0-9]+',i) ]
            TREE   = [i.split('->') for i in TREE]
            TREE   = pd.DataFrame(TREE, columns=['o', 't'], dtype='str')
            matrix = [ list(j.split('\t')[1]) + [len(j.split('\t'))-3, j.split('\t')[0]] for j in Nodes ]
            DataPd = pd.DataFrame(matrix,columns=Name + ['Values','t'])
            DataPd['Values'] = DataPd['Values']/sum( DataPd['Values'] ) * 0.9
            #DataPd.sort_values(by=Name[1:],ascending=[False]*(len(Name)-1),axis=0,inplace=True)
            #DataPd = DataPd.reset_index(drop=True)
            DataPd = pd.merge(DataPd, TREE, how ='left', on=['t'])

            Rns=[]
            for i in DataPd.index:
                Rn = [ n for n in Name if int(DataPd.loc[i,n]) ]
                Rns.append(Rn)
            DataPd['Names']=Rns
            DataPd['Length'] = DataPd['Values']/sum( DataPd['Values'] ) * 15
            DataPd['Colors'] = colors[:len(DataPd)]

            Tabular = self.HISTO(DataPd.copy())
            Forest  = self.ETree(DataPd.copy())

            Plotlist = Tabular + '\n' + Forest
            Plotlist = re.sub(r'\\',r'\\\\',Plotlist)
            tex = self.texplot(r'''%s''' %Plotlist)
            out = self.arg.outdir + '/' + SampNm + '_'  + self.arg.prefix + ".plot.tex"
            Openf(out, (tex)).openw()
            os.system('xelatex -output-directory=%s %s'%(self.arg.outdir,out))
        except ValueError:
            print('ERROR: the fomat of the file %s is wrong beacuse of the none tree in lichee!!'%file)

    def SUB(self):
        self.arg.outdir = os.path.realpath(self.arg.outdir)
        os.system('mkdir -p ' + self.arg.outdir)
        files = Openf(self.arg.input).opena()
        pool  = Pool(processes=self.arg.pool)
        pool.map(self.Get, files)
        pool.close()
        pool.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="\nTMB analysis:\n")
    parser.add_argument("-i", "--input",  type=str,  required=True, help="the input file including Id,snv.vcf.gz and bam columns")
    parser.add_argument("-o", "--outdir", type=str,  default=os.getcwd(), help="output file dir, default=current dir")
    parser.add_argument("-f", "--prefix", type=str,  default='', help="output file prefix")
    parser.add_argument("-p", "--pool",   type=int,  default=20,  help="the CPU numbers that can be used")

    global Arags
    args  = parser.parse_args()
    Arags = args.__dict__
    print('start time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print("*****The argument you have set as follows:*****")
    i = 0
    for arg in vars(args):
        i +=1
        print ('**' + str(i).zfill(2) + '|' + arg + ': ' + str(getattr(args, arg)))
    print("***********************************************")

    try:
        WorkShell(args).SUB()
        print('end time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        sys.exit(0)
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)
