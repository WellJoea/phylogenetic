#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

unless (@ARGV==4) {
	print "perl\t$0\t<sampleInfo>\t<outdir>\t<threads>\t<bedfile>\n";
	exit 1;
}

my ($info,$outdir,$thread,$bed) = @ARGV;
my $lancet = "/lustre/rde/user/guoxc/software/app/lancet/lancet";
my $ref = " /work/public/gatk_bundle/2.8/hg19/ucsc.hg19.fasta";
#my $bed = "/lustre/rde/user/guoxc/project/BED/WESG46M_only_gene.bed";

my $pm=new Parallel::ForkManager($thread);
open I,$info or die $!;
while (<I>) {
	chomp;next if(/^$/||/^\#/);
	my @item = split/\t/;
	$pm->start and  next;
	#my $sa_dir = $outdir."/"."$item[0]";
	#`mkdir -p $sa_dir`;
	#chdir($sa_dir);
	#my $vcf= $sa_dir."/".$item[0].".vcf";
#	lancet($item[1],$item[2],$vcf,$item[0]);
	mipleup($item[1],$item[0],$outdir);
	$pm->finish;
}

$pm->wait_all_children;


sub lancet{
	my ($tumor,$normal,$outfile,$prefix) = @_;
	my $cmd = "$lancet -t $tumor -n $normal -B $bed -r $ref \\
		--num-threads 10 \\
		--min-map-qual 0 \\
		--min-vaf-tumor 0.02 \\
		--max-vaf-normal 0.01 \\
		-m 2 1>$outfile 2>$prefix.log  
	";
	system($cmd);
}


sub mipleup{
	my ($bam, $id, $outdir) = @_;
	my $cmd = "samtools mpileup $bam -f $ref -L 1000000 -d 1000000 -l $bed -q 30 -Q 20 |gzip > $outdir/$id.mpileup.gz";
	system($cmd);
}

=cut 

nohup  /lustre/rde/user/guoxc/software/app/lancet/lancet -t ../BAM/B1700.sorted.mkdup.realign.bam -n 
../BAM/B17NC.sorted.mkdup.realign.bam -B ../../BioInfo_EQA_171225/CNVkit/results/Illumina_pt2.target.bed 
-r /work/public/gatk_bundle/2.8/hg19/ucsc.hg19.fasta --num-threads 20 --min-map-qual 0 --min-coverage-tumor 500 --min-vaf-tumor 0.01  1>B1700.vcf 2>lancet.log &

=cut
