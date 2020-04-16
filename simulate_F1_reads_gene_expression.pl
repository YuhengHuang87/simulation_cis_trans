#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
my $cis_trans_prop=1; #simulate cis-only effect set it to 1; simulate trans-only effect set it to 0;
my $F1_read=200; #number of reads for F1 sample
my $F1_af=resample((1,1)); #the cold allele frequency in 30 F1 offspring given the cold adapted parent is heterozygous
my $r_F0=0.7;#the ratio of the cold parental strain expression to the total expression of both parental strains
my $a_ef=(1-$r_F0)/2; #effect size on gene expression for the cold adapted allele
my $A_ef=$r_F0-$a_ef;  #effect size on gene expression for the alternative warm adapted allele

for (my $m=1; $m<=1000; $m+=1){
my $outfile1="replicate".$m."/gene_ExpProp_cis".$cis_trans_prop."_".$F1_read."_".$r_F0;
open(OUT, ">$outfile1");
my @snp_posi;
for (my $i=0; $i < 5; $i++) {
my $random_number = int(rand(100));
my $posi=$i*100+$random_number;
push @snp_posi, $posi;
}
my @cis=($F1_read*$F1_af*$A_ef,$F1_read*(1-$F1_af)*$a_ef);
my @trans=($F1_read*$F1_af,$F1_read*(1-$F1_af));
my $read_num=0;
for (my $j=0; $j < $F1_read; $j++) {
my $sel_al;my @weight;

if ($j<$F1_read*$cis_trans_prop){
@weight=(("A")x$cis[0],("a")x$cis[1]);
}else{
@weight=(("A")x$trans[0],("a")x$trans[1]);
}
$sel_al=$weight[ rand(@weight)];

my $read1_str=int(rand(500)); my $read1_end=$read1_str+75;
my $read2_str=int(rand(500-$read1_end))+$read1_end+1; my $read2_end=$read2_str+75;

my $info=0;
foreach my $posi (@snp_posi){
if ((($posi >=$read1_str)&&($posi <= $read1_end))||(($posi >=$read2_str)&&($posi <= $read2_end))){
if ($info == 0){
$read_num++;
print OUT "$read_num\t";
}
$info++;
print OUT "$posi\t";
}
}
if ($info >0){
print OUT "$sel_al\n";
}}
}

close OUT;

sub resample {
my @al=@_;
my @list=(("A")x$al[0],("a")x$al[1]);
my @resam;
for (my $i=1; $i<=30; $i+=1){
push @resam, $list[ rand(@list)];
}
my $major=0;
foreach my $resam (@resam){
if ($resam eq "A"){
$major++;
}}
my $af=$major/60;
return $af;
}
