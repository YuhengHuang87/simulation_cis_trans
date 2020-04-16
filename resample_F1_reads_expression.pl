#!/usr/bin/perl
use strict; 
use warnings;
use List::Util qw(sum);
my $cis_trans_prop=1; #cis-only effect eq 1; trans-only effect eq 0;
my $F1_read=200;
my $r_F0=0.7;


for (my $m=1; $m<=1000; $m+=1){
    for (my $i=1; $i<=1000; $i+=1){
my $outfile="replicate".$m."/bootstrap_30allele_cis".$cis_trans_prop."_".$F1_read."_".$r_F0."_".$i;
open(OUT, ">$outfile");
my $read_num=0; my @read=();
my $infile1="replicate".$m."/gene_ExpProp_cis".$cis_trans_prop."_".$F1_read."_".$r_F0;
open(FILE1,"<", "$infile1")||die"$!";
while(my $sig = <FILE1>){
	chomp($sig);
	push @read, $sig;
$read_num++;
}
	
my $new_len=scalar(@read);
my @boot;
    for (my $j=1; $j<=$new_len; $j+=1){
     my $ind= int(rand($new_len));
push @boot, $read[$ind];
}

foreach my $boot_read (@boot){
print OUT "$boot_read\n";
}
}
}
