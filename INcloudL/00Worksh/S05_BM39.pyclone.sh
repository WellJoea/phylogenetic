OUTput=/lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/05Pyclone/BM39/
mkdir -p $OUTput
/lustre/rde/user/zhouw/00software/anaconda2/bin/PyClone run_analysis_pipeline \
    --num_iters 10000 \
    --burnin 1000 \
    --plot_file_format pdf \
    --in_files \
        /lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM39_M.tsv \
        /lustre/rde/user/guoxc/project/181126_12Sample_BrM/phylogenetic/INcloudL/02SNVCNV/02SNV/BM39_P.tsv \
   --working_dir $OUTput