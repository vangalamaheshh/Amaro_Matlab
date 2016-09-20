function plot_allelic_fraction_bpdfs( maf_dir )
% Takes in list of mafs produces graphs
mafs=dir(maf_dir);
mafs(1)=[];mafs(1)=[];mafs(1)=[]; %lazy
afb=0:0.005:1;

for i=1:length(mafs)
    
    
M=load_table(sprintf('%s/%s',maf_dir,mafs(i).name)); %load each samples maf
M=rmfield(M,'header');
M=rmfield(M,'headline');
M=reorder_struct(M,~(cellfun(@isempty,M.dbSNP_RS)));
M=reorder_struct(M,(M.t_alt_count+M.t_ref_count)>100);
X.afb=repmat(afb,length(mafs),1);
X.sample=mafs(i).name;
X.af=zeros(size(X.afb));
X.afbeta=zeros(size(X.afb));
x1=zeros(size(slength(M)),201);
x2=zeros(size(slength(M)),201);
for j=1:slength(M)
     %binomial pdf for each mutation
    x1(j,:)=binopdf(M.t_alt_count(j),(M.t_alt_count(j)+M.t_ref_count(j)),afb);
    x1(j,:)=x1(j,:)/sum(x1(j,:));
    X.af(i,:)=X.af(i,:)+x1(j,:);
    x2(j,:)=betapdf(afb,M.t_alt_count(j)+1,M.t_ref_count(j)+1);
    X.afbeta(i,:)=X.afbeta(i,:)+x2(j,:);
end
subplot(4,1,1)
xlabel('allele fraction summed across mutations');
ylabel('binopdf');
plot(afb,X.af(i,:));

subplot(4,1,2)
xlabel('allele fraction')
ylabel('beta pdf summed across mutations');
plot(afb,X.afbeta(i,:));
subplot(4,1,3)
xlabel('allele fraction')
ylabel('bino pdf by mutation');
plot(afb,x1);
subplot(4,1,4)
plot(afb,x2)
ylim([0 50])
xlabel('allele fraction');
ylabel('beta pdf by mutation');

end


end



function test

plot_allelic_fraction_bpdfs('/Users/amaro/Downloads/No_Normal_Mafs/');


end