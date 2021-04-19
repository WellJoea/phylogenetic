OUTput=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/05Pyclone/BM32/
mkdir -p $OUTput
/lustre/rde/user/zhouw/00software/anaconda2/bin/PyClone run_analysis_pipeline \
    --num_iters 10000 \
    --burnin 1000 \
    --plot_file_format pdf \
    --in_files \
        /lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM32_L.tsv \
        /lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM32_M.tsv \
        /lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM32_P.tsv \
   --working_dir $OUTput