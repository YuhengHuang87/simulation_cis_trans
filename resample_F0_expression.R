cis=1; r_F0=0.7;read_num=200;
par<-read.table("0.7_parental_ExpCis1_200read_parental.txt",header=F) #need to change file name for different simulations

for (m in 1:1000){
  comb<-data.frame(1:1000)
    N=1000000;
    p=par[m,2]/N;q1=par[m,3]/N; q=par[m,4]/N; 
    readA = rbinom(1000, N, p); readB = rbinom(1000, N, q)
    readC = rbinom(1000, N, q1);
    comb<-data.frame(comb,readA,readC,readB)
    write.table(comb, file=paste("replicate",m,"/parentalRatio",cis,"_",read_num,"read_exp_parental",r_F0,sep=""),sep="\t",eol="\n",row.names = F,col.names = F)
}
    
  