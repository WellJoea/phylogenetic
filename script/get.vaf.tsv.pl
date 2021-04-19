#!/usr/bin/perl -w
use strict ;
use File::Basename  qw(dirname basename);
use Data::Dumper;


unless (@ARGV==4) {
	print "perl\t$0\t<ID>\t<anno_dir>\t<mpileip_dir>\t<outdir>\n";
	exit 1;
}
my ($info,$anno_dir,$mpile,$outdir) = @ARGV;
system("mkdir -p $outdir");
my $groupa ='';
open I,$info or die $!;
while (<I>) {
	chomp;next if(/^$/ || /^\#/ || /^Group/i);
    my $group = (split/\t/,$_)[0];
    next if ($group eq $groupa);
    $groupa = $group;
	my @anno = sort glob("$anno_dir/${group}*annovar.hg19_multianno.txt");
	print Dumper @anno;
	my (%A,%M,%C,@sam);
	for my $file(@anno) {
		my $id = (split/\-/,basename($file))[1];
		if($id !~ /BC/){
            push @sam,$id;
		    my $mp = $mpile."/".$id.".mpileup.gz";
		    anno($file,\%A,$id);
		    mpileup($mp,\%M,$id);
        }else{
		    anno($file,\%C,$id);
        }
	}

	for my $samID(@sam) {
		open O,">$outdir/$samID.vafs" or die $!;
        open T,">$outdir/$samID.tsv" or die $!;
		print O "chr\tstart\trefCount\tvarCount\tVAF\n";
		print T "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n";
        foreach my $chr (sort{my ($c)= $a =~/chr(\S+)/;my ($d) =$b=~/chr(\S+)/; $c<=>$d} keys %A) {
        #foreach my $chr ( keys %A) {
			foreach my $pos (sort keys %{$A{$chr}}){ 
                if (exists $C{$chr}{$pos}){
                    print "$chr\t$pos: is in BC file"
                }elsif(exists $A{$chr}{$pos}{$samID}){
					print O "$chr\t$pos\t$A{$chr}{$pos}{$samID}{'R'}\t$A{$chr}{$pos}{$samID}{'V'}\t$A{$chr}{$pos}{$samID}{'Ra'}\n";
                    print T "$chr:$pos\t$A{$chr}{$pos}{$samID}{'R'}\t$A{$chr}{$pos}{$samID}{'V'}\t2\t0\t2\n";
				}elsif(exists $M{$chr}{$pos}{$samID}){
					print O "$chr\t$pos\t$M{$chr}{$pos}{$samID}\t0\t0\n";
				}else{
					print O "$chr\t$pos\t0\t0\t0\n";
				}
			}
		}
        close O;
        close T;
        system("/lustre/rde/user/zhouw/00software/anaconda2/bin/PyClone  build_mutations_file --in_file $outdir/$samID.tsv --out_file $outdir/$samID.yaml");
	}
}
close I;

###################sub_read_anno#########
sub anno {
	my ($file,$ref,$samID) = @_;
	open A,$file or die $!;
	while (<A>) {
		chomp;next if(/^$/ || /^\#/ || /^Chr/);
		my @item = split/\t/;
		if($item[5] eq 'exonic'){
            my $flag = 1;
            for my $index(11..18) {
                next if($item[$index] eq '.');
                $flag++ if($item[$index] >0.01);
            }
            for my $index(19..27) {
               next if($item[$index] eq '.');
               $flag++ if($item[$index] >0.01);
            }
            next if($flag!=1);

            my ($Ref, $var, $ra);
            if($item[-2]!~/0\/0/){
                ($Ref,$var) = (split/,/,(split/:/,$item[-2])[1]);
            }elsif($item[-1]!~/0\/0/){
                ($Ref,$var) = (split/,/,(split/:/,$item[-1])[1]);
            }else{
                print "The sample $samID caonnot find the site :$item[0]\t$item[1] type!!\n";
            }
            $ra = sprintf("%.4f",$var/($Ref+$var)*100);
            #if($ra>0.005 && $var>=4 && ($Ref+$var)>=20){
            if($ra>0.05 && $var>=5 && ($Ref+$var)>=40){
			    $ref->{$item[0]}{$item[1]}{$samID}{'R'}=$Ref;
			    $ref->{$item[0]}{$item[1]}{$samID}{'V'}=$var;
			    $ref->{$item[0]}{$item[1]}{$samID}{'Ra'}=$ra;
            }
		}
	}close A;
}

####################parse_mpileup#########
sub mpileup{
	my ($file,$ref,$samID) = @_;
	open M,"gzip -dc $file|" or die $!;
	while (<M>) {
		my @s = split/\t/;
		$ref->{$s[0]}{$s[1]}{$samID}=$s[3];
	}
	close M;
}
