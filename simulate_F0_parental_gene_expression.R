r_F0=c(0.7); 
num=c(200);
cis=1; #1 for cis-only simulation; 0 for trans-only simulation
for (read_num in num){
  for (r1 in r_F0){
    comb<-data.frame(1:1000)
    r2=1-r1;
if (cis == 1){
    pA= r1-r2*0.5;
    A=pA*read_num;
    B=r2*read_num;
    C=r2*read_num*0.5;
  }
if (cis == 0){
    A=r1*read_num*0.5; 
    B=r2*read_num;
    C=r1*read_num*0.5;
}
  
    N=1000000;
    p=A/N;q1=C/N; q=B/N; 
    readA = rbinom(1000, N, p); readB = rbinom(1000, N, q)
    readC = rbinom(1000, N, q1);
    comb<-data.frame(comb,readA,readC,readB)
    write.table(comb, file=paste(r1,"_parental_ExpCis",cis,"_",read_num,"read_parental.txt",sep=""),sep="\t",eol="\n",row.names = F,col.names = F)
  }}

