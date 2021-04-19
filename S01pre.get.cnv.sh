IN=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic
OU=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/cnv
CNV=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/03CNV
#/lustre/rde/user/zhouw/00software/anaconda3/bin/python  $IN/script/multiprocess.shell.gatk.cnvkit.py -i $IN/sample.bam.infor.txt -o $OU
for i in `cut -f 1 $IN/sample.bam.infor.txt`;
do
    grep -v "chromosome" $OU/01CNVkit/$i/${i}.call.cns | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$7,$6}' > $CNV/${i}.cnv.txt
done
