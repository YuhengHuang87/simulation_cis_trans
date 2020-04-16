Set up simulation folders. enter the commands in the working directary to generate 1000 folders:

n=1
max=1000
set -- # this sets $@ [the argv array] to an empty list.

while [ "$n" -le "$max" ]; do
    set -- "$@" "replicate$n" # this adds s$n to the end of $@
    n=$(( $n + 1 ));
done 
mkdir "$@"

simulation phrase
1. run Rscript simulate_F0_parental_gene_expression.R 
2. run perl simulate_F1_reads_gene_expression.pl

analysis phrase
1. run Rscript resample_F0_expression.R to resample reads for F0 expression
2. run perl resample_F1_reads_expression.pl to resample reads for F1 expression
3. run perl cis_trans_effect.pl to compute cis- and trans- effect for each replicate
4. run perl interval_significant.pl to determine the numbers of cis- and trans- effects
