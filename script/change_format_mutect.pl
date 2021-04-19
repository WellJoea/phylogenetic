#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Data::Dumper;
use feature qw/say/;

unless (@ARGV>=2) {
	print "perl\t$0\t<mutation.list.txt>\t<outpath>\t<change>\n";
	exit 1;
}

my %G  = ();
my %n  = ();
my %k  = ();
if ($ARGV[2]){
    open IF,$ARGV[2] or die $!;
        %k=map{chomp;my @cut=split/\t/,$_ ; $cut[0] =>$cut[1] ; }<IF>;
    close IF;
}

open IN,$ARGV[0] or die $!;
while(<IN>){
	chomp;next if /^$|^#/;
	my ($grop,$filename) = (split/[-]/,basename($_))[0,1];
#my ($filename) =join('_', (split/[_.]/,basename($_))[1,2]);
#	my ($filename) =(split/[.]/,basename($_))[0];
    say $filename;
    ($ARGV[2]) && ($filename = $k{$filename}) ;
	push @{$n{$grop}}, $filename;
	open I,$_ or die $!;
	while(<I>){
		chomp;
		next if /^$|^Chr|intronic/;
		my @s = split/\t/;
		my $flag=1;
		if(/PASS/ && /exonic/i){
			for my $n (11..27) {
				next if($s[$n] eq '.');
				$s[$n]>=0.002 && $flag++ ;
			}
#			next if($s[8] eq 'synonymous SNV');
			next if($flag!=1);

            my $aa = ($s[-1] !~/0\/0/ ) ? $s[-1] :
                     ($s[-2] !~/0\/0/ ) ? $s[-2] : '???' ;
            my @d = split/:/,$aa;

			my ($ref,$alt) = split/,/,$d[1];
			my $ra = sprintf("%.4f",$alt/($ref+$alt));
			if($ra>0.001 && $alt>=1 && ($ref+$alt)>=10){
				$G{$grop}{join "\t", @s[0,1,6,9]}{$filename} = "$ra\t$aa" ;
                #@{$G{$grop}{join "\t", @s[0,1,6]}{$filename}} = ($ra, @s[0..6,9,43]);
			}
		}
	}
	close I;
}
close IN;

open  J,">$ARGV[1]/all_mut_evol.list.txt" or die $!;
foreach my $pid (sort keys %G){
    my %count;
    my @uniq = sort(grep { ++$count{ $_ } < 2; } @{$n{$pid}});
    my @titl = map{ "$_\t${_}_info" }@uniq;
	open O, ">$ARGV[1]/${pid}_mutect_evol.txt" or die $!;
	open P, ">$ARGV[1]/${pid}_mutect_evol_massage.xls" or die $!;
	say  O "#chrom\tpos\tDESC\tnormal\t".join("\t",@uniq);
	say  P "#chrom\tpos\tDESC\tAAchange\tnormal\t".join("\t",@titl);
    say  J "$ARGV[1]/${pid}_mutect_evol.txt";
	foreach my $info (sort keys %{$G{$pid}}){
		my @ramatrix = (0);
		my @Pamatrix = (0);
		foreach my $samp ( @uniq ){
			my $ra = (exists $G{$pid}{$info}{$samp}) ? $G{$pid}{$info}{$samp} : "0\t-" ;
			push @ramatrix, (split/\t/,$ra)[0];
            push @Pamatrix , $ra;
		}
		say O join"\t",((split/\t/,$info)[0..2],@ramatrix);
		say P join"\t",($info,@Pamatrix);
	}
	close O;
	close P;
}
close J;
