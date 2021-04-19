while read a b
do
  ln -s /lustre/rde/user/guoxc/project/180824_RuiJ_pancreatic/CNVKIT/$b/${b}.merge.cns ./
  egrep -v "CHROM|ENSG" /lustre/rde/user/guoxc/project/180824_RuiJ_pancreatic/CNVKIT/$b/${b}.merge.cns |awk 'BEGIN{FS="\t";OFS="\t"}{if($4 !="-") print $1,$2,$3,$8,$6}' > $a.cnv.txt
done < zxm.id.txt
