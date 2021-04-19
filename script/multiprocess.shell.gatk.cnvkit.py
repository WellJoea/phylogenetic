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
import os
import numpy as np
import pandas as pd
import time
import glob
from multiprocessing import Process, Pool
import multiprocessing.pool

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
    def __init__(self, file, *array):
        self.file  = file
        self.array = array
    def opena(self):
        with open(self.file, 'r', encoding='utf-8') as inp:
            #outa = inp.read().strip().split('\n')
            outa = [i.strip('\n') for i in inp.readlines() if not re.match(r'#',i)]
        return outa
    def openb(self):
        if re.match(r'.*xlsx$|.*xlss$', self.file, re.I):
            outp = pd.read_excel(self.file, header=0, encoding='utf-8', sheet_name=0).fillna('NaN')
        else:
            outp = pd.read_table(self.file, header=0, encoding='utf-8', sep='\t', comment=None).fillna('NaN')
        return(pd.DataFrame(outp))
    def openw(self):
        f= open(self.file, 'w',encoding='utf-8')
        f.write(self.array[0])
        f.close()

class WorkShell():
    def __init__(self, arg, *array):
        self.arg = arg
        self.array = array

    def gatk(self, file):
        Bam    = file.strip()
        Output = self.arg.output
        Sample = os.path.basename(file.strip()).split('_')[0]

        worksh =  \
r'''
Bam=%[Bam]%
Sample=%[Sample]%
Output=%[Output]%

REF=/work/public/gatk_bundle/2.8/hg19/ucsc.hg19.fasta
GatkBD=/lustre/work/public/gatk_bundle/2.8/hg19/
BED=/lustre/rde/user/zhouw/02database/03BED/medexome.bed
ANnovDB=/lustre/rde/user/zhouw/00software/annovar/humandb/

OutID=$Output/$Sample
SampleID=$Output/$Sample/$Sample

mkdir -p $OutID

/lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk BaseRecalibrator \\
        -R $REF \\
        -I $Bam \\
        -O ${SampleID}.sort_markdup.grp \\
        --known-sites $GatkBD/dbsnp_138.hg19.vcf \\
        --known-sites $GatkBD/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \\
        --known-sites $GatkBD/1000G_phase1.indels.hg19.sites.vcf \\
        --known-sites $GatkBD/1000G_phase1.snps.high_confidence.hg19.sites.vcf

    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk  ApplyBQSR \\
        -R $REF \\
        -I $Bam  \\
        -O ${SampleID}.sort_markdup.grp.bam \\
        -bqsr ${SampleID}.sort_markdup.grp

    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk  AnalyzeCovariates \\
        -bqsr ${SampleID}.sort_markdup.grp
        -plots ${SampleID}.sort_markdup.grp.compa.pdf

    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk  PrintReads \\
        -R $REF \\
        -I ${SampleID}.sort_markdup.grp.bam \\
        -O ${SampleID}.sort_markdup.grp.call.bam
    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk  HaplotypeCaller \\
        -R $REF \\
        -L $BED \\
        -I ${SampleID}.sort_markdup.grp.call.bam \\
        --genotyping-mode DISCOVERY -stand-call-conf 10 \\
        -O ${SampleID}.gatk.call.vcf

perl /lustre/rde/user/zhouw/00software/annovar/table_annovar.pl \\
        ${SampleID}.gatk.call.vcf \\
        $ANnovDB \\
        --outfile ${SampleID}.gatk.annovar \\
        --buildver hg19 --otherinfo --remove --nastring . --thread 10 --vcfinput \\
        --operation g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r  \\
        --argument '--splicing_threshold 10 --exonicsplicing','--splicing_threshold 10 --exonicsplicing',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, \\
        --protocol refGene,ensGene,knownGene,cytoBand,snp138,snp150,avsnp150add,gnomad_exome,gnomad_genome,exac03,esp6500siv2_all,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,1000g2015aug_amr,1000g2015aug_chb,1000g2015aug_chs,popfreq_max_20150413,gwasCatalog,wgRna,targetScanS,tfbsConsSites,simpleRepeat,genomicSuperDups,dbscsnv11,dbnsfp33a,intervar_20180118,kaviar_20150923,hrcr1,mcap,revel,clinvar_20180603,icgc21,cosmic80,HGMD,HGMD_mut,rmsk,BIC,OMIM,reliability,Pro_CancerRepeat  #mitimpact24

python /lustre/rde/user/zhouw/21Project/05gatkcall/germline.somatic.filter.site.py \\
    -i ${SampleID}.gatk.annovar.hg19_multianno.txt \\
    -o $SampleID \\
    -g /lustre/rde/user/zhouw/02database/03BED/47genelist.txt

<<'COMMENT'
    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk VariantRecalibrator \\
        -R $REF \\
        -V ${SampleID}.gatk.call.vcf \\
        -O ${SampleID}.gatk.recalibrate_SNP.recal \\
        --resource hapmap,known=false,training=true,truth=true,prior=15.0:$GatkBD/hapmap_3.3.hg19.sites.vcf \\
        --resource omni,known=false,training=true,truth=true,prior=12.0:$GatkBD/1000G_omni2.5.hg19.sites.vcf \\
        --resource 1000G,known=false,training=true,truth=false,prior=10.0:$GatkBD/1000G_phase1.snps.high_confidence.hg19.sites.vcf \\
        --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$GatkBD/dbsnp_138.hg19.excluding_sites_after_129.vcf \\
        -an DP \\
        -an QD \\
        -an FS \\
        -an SOR \\
        -an MQ \\
        -an MQRankSum \\
        -an ReadPosRankSum \\
        -mode SNP \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        --tranches-file ${SampleID}.gatk.recalibrate_SNP.tranches \\
        --rscript-file  ${SampleID}.gatk.recalibrate_SNP_plots.R #-AS 
        #-an InbreedingCoeff \\
        #--max-gaussians = 4 \\
    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk ApplyVQSR \\
        -R $REF \\
        -V ${SampleID}.gatk.call.vcf \\
        -O ${SampleID}.gatk.snp.recalibrated.vcf \\
        --recal-file ${SampleID}.gatk.recalibrate_SNP.recal \\
        --tranches-file ${SampleID}.gatk.recalibrate_SNP.tranches \\
        --truth-sensitivity-filter-level 99.0 \\
        -mode SNP
    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk VariantRecalibrator \\
        -R $REF \\
        -V ${SampleID}.gatk.call.vcf \\
        -O ${SampleID}.gatk.recalibrate_INDEL.recal \\
        --resource mills,known=false,training=true,truth=true,prior=12.0:$GatkBD/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf  \\
        --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$GatkBD/dbsnp_138.hg19.excluding_sites_after_129.vcf \\
        -an QD \\
        -an DP \\
        -an FS \\
        -an SOR \\
        -an MQ \\
        -an MQRankSum \\
        -an ReadPosRankSum \\
        -mode INDEL \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        --tranches-file ${SampleID}.gatk.recalibrate_INDEL.tranches \\
        --rscript-file ${SampleID}.gatk.recalibrate_INDEL_plots.R \\
        --max-gaussians = 4 
        #-an InbreedingCoeff \\
    /lustre/rde/user/zhouw/00software/gatk-4.0.5.1/gatk ApplyVQSR \\
        -R $REF \\
        -V ${SampleID}.gatk.call.vcf \\
        -O ${SampleID}.gatk.indel.recalibrated.vcf \\
        --recal-file ${SampleID}.gatk.recalibrate_INDEL.recal \\
        --tranches-file ${SampleID}.gatk.recalibrate_INDEL.tranches \\
        --truth-sensitivity-filter-level 99.0 \\
        -mode INDEL
COMMENT
'''
        worksh = re.sub(r'%\[Bam\]%'   , Bam   , worksh)
        worksh = re.sub(r'%\[Sample\]%', Sample, worksh)
        worksh = re.sub(r'%\[Output\]%', Output, worksh)

        os.system('mkdir -p ' + Output + '/01WorkSH/')
        f = open(Output + '/01WorkSH/' + Sample + '.gatk.work.sh', 'w', encoding='utf-8')
        f.write(worksh)
        f.close()
        os.system('sh ' + Output + '/01WorkSH/' + '.gatk.work.sh')

    def cnvkit(self, file):
        (Sample,TYpe,BED,Bam,NBam) = file.split("\t")
        BamID  = '.'.join(os.path.basename(Bam).split('.')[:-2])
        Output = self.arg.output + '/01CNVkit/'
        worksh =  \
r'''
#python2=/lustre/rde/user/zhouw/00software/anaconda2/bin/python2
python2=python2
cnvkit=/lustre/rde/user/zhouw/00software/cnvkit/cnvkit.py
REF=/work/public/gatk_bundle/2.8/hg19/ucsc.hg19.fasta
CNVGEN=/lustre/rde/user/zhouw/02database/02UCSC/refFlat.txt
CNVacces=/lustre/rde/user/zhouw/00software/cnvkit/data/

Bam=%[Bam]%
OUTput=%[Output]%
SampleID=%[Sample]%
NBam=%[NBam]%
BED=%[BED]%
BamID=%[BamID]%
Ckitout=$OUTput/$SampleID/
mkdir -p $Ckitout

if [ ! -f "$Ckitout/${BamID}.cns" ]
then
if [[ "$NBam" ]]
then
$python2 ${cnvkit} batch \
    $Bam -n $NBam -c --drop-low-coverage \
    -t ${BED} -f ${REF} \
    --access ${CNVacces}/access-5k-mappable.hg19.bed \
    --output-reference $Ckitout/${SampleID}_reference.cnn \
    --output-dir $Ckitout/ --scatter --diagram  -p 5
    #--annotate ${CNVGEN} \
else
$python2 ${cnvkit} batch \
    $Bam -n -c --drop-low-coverage \
    -t ${BED} -f ${REF} \
    --access ${CNVacces}/access-5k-mappable.hg19.bed \
    --output-reference $Ckitout/${SampleID}_reference.cnn \
    --output-dir $Ckitout/ --scatter --diagram  -p 5
    #--annotate ${CNVGEN} \
fi
fi

$python2 ${cnvkit} call \
    $Ckitout/${BamID}.cns \
    -o $Ckitout/${SampleID}.call.cns \

if [[ "$NBam" ]]
then
$python2 ${cnvkit} genemetrics \
    -m 5 -t 0 \
    -o $Ckitout/${SampleID}.genemetrics.txt  \
    -s $Ckitout/${SampleID}.call.cns \
    $Ckitout/${BamID}.cnr
else
$python2 ${cnvkit} genemetrics \
    -m 5 -t 0.2  \
    -o $Ckitout/${SampleID}.genemetrics.txt  \
    -s $Ckitout/${SampleID}.call.cns \
    $Ckitout/${BamID}.cnr
fi

$python2 ${cnvkit} breaks \
    $Ckitout/${BamID}.cnr \
    $Ckitout/${SampleID}.call.cns \
    -o $Ckitout/${SampleID}.breaks

$python2 ${cnvkit} diagram \
    -s $Ckitout/${SampleID}.call.cns \
    $Ckitout/${BamID}.cnr \
    -o $Ckitout/${SampleID}.diagram.pdf

$python2 ${cnvkit} sex \
    $Ckitout/${SampleID}_reference.cnn \
    $Ckitout/${BamID}.cnr \
    $Ckitout/${SampleID}.call.cns \
    -o $Ckitout/${SampleID}.sex
'''
        worksh = re.sub(r'%\[Bam\]%'   , Bam   , worksh)
        worksh = re.sub(r'%\[Sample\]%', Sample, worksh)
        worksh = re.sub(r'%\[Output\]%', Output, worksh)
        worksh = re.sub(r'%\[NBam\]%'  , NBam  , worksh)
        worksh = re.sub(r'%\[BED\]%'   , BED   , worksh)
        worksh = re.sub(r'%\[BamID\]%' , BamID , worksh)

        os.system('mkdir -p ' + self.arg.output + '/00WorkSH/')
        Openf(self.arg.output + '/00WorkSH/' + Sample + '.cnvkit.work.sh', (worksh)).openw()
        os.system('qsub -cwd ' + self.arg.output + '/00WorkSH/' + Sample + '.cnvkit.work.sh')

    def contra(self, file):
        (Sample,TYpe,BED,Bam,NBam) = file.split("\t")
        Output = self.arg.output + '/02CONTRA/'
        worksh =  \
r'''
python2=/lustre/rde/user/zhouw/00software/Python-2.7.13/bin/python2
contra=/lustre/rde/user/zhouw/00software/CONTRA.v2.0.8/contra.py
REF=/work/public/gatk_bundle/2.8/hg19/ucsc.hg19.fasta

Bam=%[Bam]%
OUTput=%[Output]%
SampleID=%[Sample]%
NBam=%[NBam]%
BED=%[BED]%
Ctraout=$OUTput/$SampleID/
mkdir -p $Ctraout

if [[ "$NBam" ]]
then
$python2 ${contra} \
    -t $BED \
    -s $Bam \
    -c $NBam \
    -o $Ctraout \
    -f $REF
fi
'''
        worksh = re.sub(r'%\[Bam\]%'   , Bam   , worksh)
        worksh = re.sub(r'%\[Sample\]%', Sample, worksh)
        worksh = re.sub(r'%\[Output\]%', Output, worksh)
        worksh = re.sub(r'%\[NBam\]%'  , NBam  , worksh)
        worksh = re.sub(r'%\[BED\]%'   , BED   , worksh)

        if TYpe =='somatic':
            os.system('mkdir -p ' + self.arg.output + '/00WorkSH/')
            Openf(self.arg.output + '/00WorkSH/' + Sample + '.contra.work.sh', (worksh)).openw()
            os.system('sh ' + self.arg.output + '/00WorkSH/' + Sample + '.contra.work.sh')

    def SUB(self):
        file = Openf(self.arg.input).opena()
        self.arg.output=os.path.realpath(self.arg.output)
        for i in re.split(r'[,;]', self.arg.type):
            pool = Pool(processes=self.arg.pool)
            pool.map(eval('self.'+i), file)
            pool.close()
            pool.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input"  , "-i", type=str, help="the input file used for filter")
    parser.add_argument("--output" , "-o", type=str, default=os.getcwd(),help="the output file name header")
    parser.add_argument("--pool"   , "-p", type=int, default=10, help="the CPU numbers that can be used")
    parser.add_argument("--type"   , "-t", type=str, default='cnvkit', help="the software you need to ananlyse,such as: cnvkit,contra,gatk")

    global Arags
    args  = parser.parse_args()
    Arags = args.__dict__
    print(args)
    print('start time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    WorkShell(args).SUB()
    print('end time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
