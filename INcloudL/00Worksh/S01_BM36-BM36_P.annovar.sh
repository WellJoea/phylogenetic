FILE=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/mutect2/BM36_P.filter.vcf.gz
Annovfi=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/01ANNOVAR/BM36-BM36_P-annovar
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