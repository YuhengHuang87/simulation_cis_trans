#!/usr/bin/perl
use strict; 
use warnings;
my $F1_read=200;
my $cis_trans_prop=1; #cis-only effect eq 1; trans-only effect eq 0;
my $r_F0=0.7;
my $al_num=60;
my $ini_h=0.5;my $ini_l=0; #if both parents are homozygous, set my $ini_h=1;my $ini_l=0; 

###analysis for cis- effect
for (my $m=1; $m<=1000; $m+=1){
my $obs_p_h=resample($ini_h*$al_num,$al_num-$ini_h*$al_num);
my $obs_p_l=resample($ini_l*$al_num,$al_num-$ini_l*$al_num);
my $al_h=$obs_p_h*$al_num;my $al_l=$obs_p_l*$al_num;

my $outfile="replicate".$m."/bootstrap_ap_gene_resample30allele_".$F1_read."_cis".$cis_trans_prop."_ExpProp".$r_F0."_for_cis";
open(OUT, ">$outfile");

for (my $i=1; $i<=1000; $i+=1){
my %alle_h; my %alle_l;
my $infile1="replicate".$m."/bootstrap_30allele_cis".$cis_trans_prop."_".$F1_read."_".$r_F0."_".$i;
open(FILE1,"<", "$infile1")||die"$!"; 
while(my $sig = <FILE1>){
	chomp($sig);
	my @a = split("\t", $sig);
    for (my $k=1; $k<@a-1; $k+=1){
my $s=$a[$k];
if ($a[-1] eq "A"){
if (exists $alle_h{$s}){
$alle_h{$s}=$alle_h{$s}+1;
}else{
$alle_h{$s}=1;
}
}elsif($a[-1] eq "a"){
if (exists $alle_l{$s}){
$alle_l{$s}=$alle_l{$s}+1;
}else{
$alle_l{$s}=1;
}
}}}
my %counted_loci;my @ac_ap=();
my $sum_ap=0; my $num_site=0;my $exp_all_h=0; my $exp_all_l=0;
foreach my $key (keys %alle_h) {
$exp_all_h=$alle_h{$key};
$counted_loci{$key}=1;
if (exists $alle_l{$key}){
$exp_all_l=$alle_l{$key};
}else{
$exp_all_l=0;
}
if(($exp_all_h+$exp_all_l)>=10){
my $af_f1=$exp_all_h/($exp_all_h+$exp_all_l);

my $sam_h=resample($al_h/2,($al_num-$al_h)/2);
my $sam_l=resample($al_l/2,($al_num-$al_l)/2);

my $paf_h=resample($al_num*$sam_h/2,$al_num*(1-$sam_h)/2);
my $paf_l=resample($al_num*$sam_l/2,$al_num*(1-$sam_l)/2);
if($paf_h!=$paf_l){

my $anc_prop=($af_f1-$paf_l)/($paf_h-$paf_l);
push @ac_ap, $anc_prop;
$sum_ap+=$anc_prop;
$num_site++;
}}
}

foreach my $key (keys %alle_l) {
$exp_all_l=$alle_l{$key};
if (exists $counted_loci{$key}){
}else{
$exp_all_h=0;
if(($exp_all_h+$exp_all_l)>=10){
my $af_f1=$exp_all_h/($exp_all_h+$exp_all_l);
my $sam_h=resample($al_h/2,($al_num-$al_h)/2);
my $sam_l=resample($al_l/2,($al_num-$al_l)/2);
my $paf_h=resample($al_num*$sam_h/2,$al_num*(1-$sam_h)/2);
my $paf_l=resample($al_num*$sam_l/2,$al_num*(1-$sam_l)/2);

if($paf_h!=$paf_l){
my $anc_prop=($af_f1-$paf_l)/($paf_h-$paf_l);
push @ac_ap, $anc_prop;
$sum_ap+=$anc_prop;
$num_site++;
}}}
}
if ($num_site>0){
my $ap_mean=$sum_ap/$num_site;my $ap_median=median(@ac_ap);
print OUT "$i\t","$ap_mean\t","$ap_median\t","$num_site\n";
}}
close OUT;
}

###analysis for trans- effect
 for (my $m=1; $m<=1000; $m+=1){
my $al_h=$ini_h*$al_num;my $al_l=$ini_l*$al_num;
my $al_num_h; my $al_num_l;
my $a_ef;my $A_ef;
my $file=$r_F0."_parental_ExpCis".$cis_trans_prop."_".$F1_read."read_parental.txt";
open(FILE,"<", "$file")||die"$!";
while(my $gene = <FILE>){
		chomp($gene);
		my @a = split("\t", $gene);
if ($a[0] == $m){
$al_num_h=$a[1]+$a[2];$al_num_l=$a[3];
$A_ef=$a[1];$a_ef=$a[2];
}}
my $outfile="replicate".$m."/bootstrap_ap_gene_resample_30allele_".$F1_read."_cis".$cis_trans_prop."_ExpProp".$r_F0."_for_trans";
open(OUT, ">$outfile");

    for (my $i=1; $i<=1000; $i+=1){
my %alle_h; my %alle_l;
my $infile1="replicate".$m."/bootstrap_30allele_cis".$cis_trans_prop."_".$F1_read."_".$r_F0."_".$i;
open(FILE1,"<", "$infile1")||die"$!"; 
while(my $sig = <FILE1>){
	chomp($sig);
	my @a = split("\t", $sig);
    for (my $k=1; $k<@a-1; $k+=1){
my $s=$a[$k];
if ($a[-1] eq "A"){
if (exists $alle_h{$s}){
$alle_h{$s}=$alle_h{$s}+1;
}else{
$alle_h{$s}=1;
}
}elsif($a[-1] eq "a"){
if (exists $alle_l{$s}){
$alle_l{$s}=$alle_l{$s}+1;
}else{
$alle_l{$s}=1;
}
}
}}
my %counted_loci;my @ac_ap=();
my $sum_ap=0; my $num_site=0;my $exp_all_h=0; my $exp_all_l=0;
foreach my $key (keys %alle_h) {
$exp_all_h=$alle_h{$key};
$counted_loci{$key}=1;
if (exists $alle_l{$key}){
$exp_all_l=$alle_l{$key};
}else{
$exp_all_l=0;
}
if(($exp_all_h+$exp_all_l)>=10){
my $af_f1=$exp_all_h/($exp_all_h+$exp_all_l);
my $sam_h=resample($al_h,$al_num-$al_h);
my $sam_l=resample(($al_l,$al_num-$al_l));
my $p_h=($sam_h*$A_ef)/($sam_h*$A_ef+(1-$sam_h)*$a_ef);my $p_l=($sam_l*$A_ef)/($sam_l*$A_ef+(1-$sam_l)*$a_ef); 
my $sam_al_h=$p_h*$al_num_h;my $sam_al_l=$p_l*$al_num_l;
my $paf_h=resample($sam_al_h/2,($al_num_h-$sam_al_h)/2);
my $paf_l=resample($sam_al_l/2,($al_num_l-$sam_al_l)/2);

if($paf_h!=$paf_l){
my $anc_prop=($af_f1-$paf_l)/($paf_h-$paf_l);
push @ac_ap, $anc_prop;
$sum_ap+=$anc_prop;
$num_site++;
}
}}

foreach my $key (keys %alle_l) {
$exp_all_l=$alle_l{$key};
if (exists $counted_loci{$key}){
}else{
$exp_all_h=0;
if(($exp_all_h+$exp_all_l)>=10){
my $af_f1=$exp_all_h/($exp_all_h+$exp_all_l);
my $sam_h=resample($al_h,$al_num-$al_h);
my $sam_l=resample(($al_l,$al_num-$al_l));
my $p_h=($sam_h*$A_ef)/($sam_h*$A_ef+(1-$sam_h)*$a_ef);my $p_l=($sam_l*$A_ef)/($sam_l*$A_ef+(1-$sam_l)*$a_ef); 
my $sam_al_h=$p_h*$al_num_h;my $sam_al_l=$p_l*$al_num_l;
my $paf_h=resample($sam_al_h/2,($al_num_h-$sam_al_h)/2);
my $paf_l=resample($sam_al_l/2,($al_num_l-$sam_al_l)/2);
if($paf_h!=$paf_l){
my $anc_prop=($af_f1-$paf_l)/($paf_h-$paf_l);
push @ac_ap, $anc_prop;
$sum_ap+=$anc_prop;
$num_site++;
}
}}}
if ($num_site>0){
my $ap_mean=$sum_ap/$num_site; my $ap_median=median(@ac_ap);
print OUT "$i\t","$ap_mean\t","$ap_median\t","$num_site\n";
}}
}


sub resample {
my @al=@_;
my @list=(("A")x (int $al[0]),("a")x(int $al[1]));
my $num = scalar(@list);
my @offspring_resam;
for (my $i=1; $i<=$num; $i+=1){
push @offspring_resam, $list[ rand(@list)];
}
my $major=0;
foreach my $offspring_resam (@offspring_resam){
if ($offspring_resam eq "A"){
$major++;
}}
my $af=$major/$num;
return $af;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}