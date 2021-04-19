INput=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/01ANNOVAR/All.annovar.files.txt
OUTput=`realpath /lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL`
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