#!/usr/bin/perl
use strict; 
use warnings;
my $F1_read=200;
my $cis_trans_prop=1; #cis-only effect simulation 1; trans-only effect simulation 0;
my $r_F0=0.7;

my $outfile="interval_simu_F1_exp_ance_prop_cis".$cis_trans_prop."_".$F1_read."_rF0_".$r_F0;
open(OUT, ">$outfile");
print OUT "Read_num\t","r_F0\t","cis_num\t","noncis_num\t","trans_num\t","nontrans_num\n";


my $cis_num=0;my $trans_num=0; 
my $non_cis=0;my $non_trans=0;


for (my $m=1; $m<=1000; $m+=1){
my %boot_prop;
my $infile1="replicate".$m."/bootstrap_ap_gene_resample30allele_".$F1_read."_cis".$cis_trans_prop."_ExpProp".$r_F0."_for_cis";
open(FILE,"<", "$infile1")||die"$!"; 
while(my $sig = <FILE>){
	chomp($sig);
	my @a = split("\t", $sig);
my $anc_pro=$a[1];
$boot_prop{$a[0]}=$anc_pro;
}
my $test=0;my @boot_ap;
foreach my $key (keys %boot_prop) {
push @boot_ap, $boot_prop{$key};
$test++;
}
if($test>=50){
my @sorted_boot = sort { $a <=> $b } @boot_ap;
my $len_b=scalar(@sorted_boot);
my $low_b = $len_b*0.025; my $high_b=$len_b*0.975;

if (($sorted_boot[$low_b]>0.5)||($sorted_boot[$high_b]<0.5)){
$cis_num++;
}else{
$non_cis++;
}}}

for (my $m=1; $m<=1000; $m+=1){
my %resam_parent;

my $infile2="replicate".$m."/parentalRatio".$cis_trans_prop."_".$F1_read."read_exp_parental".$r_F0;#parentalRatio0_200read_exp_prop0.35
open(FILE2,"<", "$infile2")||die"$!"; 
while(my $count = <FILE2>){
chomp($count);
	my @c = split("\t", $count);
$resam_parent{$c[0]}=($c[1]+$c[2])/($c[1]+$c[2]+$c[3]);
}
my %boot_prop;

my $infile1="replicate".$m."/bootstrap_ap_gene_resample_30allele_".$F1_read."_cis".$cis_trans_prop."_ExpProp".$r_F0."_for_trans";
open(FILE,"<", "$infile1")||die"$!"; 
while(my $sig = <FILE>){
	chomp($sig);
	my @a = split("\t", $sig);
my $anc_pro=$a[1];
$boot_prop{$a[0]}=$anc_pro;
}

my $test=0; my @dif;
foreach my $key (keys %boot_prop) {
my $dif_ac_pp=$boot_prop{$key}-$resam_parent{$key};
push @dif, $dif_ac_pp;
$test++;
}
if($test>=50){

my @sorted_dif = sort { $a <=> $b } @dif;
my $len=scalar(@sorted_dif);
my $low = $len*0.025; my $high=$len*0.975;
if (($sorted_dif[$low]>0)||($sorted_dif[$high]<0)){
$trans_num++;
}else{
my $non_trans++;
}
}
}

print OUT "$F1_read\t","$r_F0\t","$cis_num\t","$non_cis\t","$trans_num\t","$non_trans\n";
