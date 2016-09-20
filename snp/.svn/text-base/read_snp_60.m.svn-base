function [H1,H2]=read_snp_60(datadir,snp_fname,ext,genome_info_file)

%datadir='/xchip/cancergenome/1chip/data/Nsp-Sty_Mixes/';
%snp_fname='6.0_nspsty_mixes.snp';
%ext='_070630';

fname=[datadir snp_fname];
% unix(['wc -l ' datadir snp_fname]);
if exist([fname ext '.dat.mat'],'file')
  disp('Reading existing mat files ...')
  load([fname ext '.dat.mat']); H.dat=x;
  load([fname ext '.sdesc.mat']); H.sdesc=x;
  load([fname ext '.marker.mat']); H.marker=cellstr(x);
else
  H=read_modelled_data_file(fname,-1,-1,0,0,0,0,0,5e5); 
  disp('Writing mat files ...');
  x=H.dat; save([fname ext '.dat.mat'],'x','-v7.3');
  x=H.sdesc; save([fname ext '.sdesc.mat'],'x','-v7.3');
  x=strvcat(H.marker); save([fname ext '.marker.mat'],'x','-v7.3');
end

H1=H;

m=strvcat(H.marker);
msz=sum(m~=' ',2);
last_char_idx=sub2ind(size(m),(1:size(m,1))',msz);
pre_lst_idx=sub2ind(size(m),(1:size(m,1))',msz-1);
m_no_AB=m;
m_no_AB(last_char_idx)=' ';
m_no_AB(pre_lst_idx)=' ';

Hsnp_idx=find(m(:,1)=='S');
Hcn_idx=find(m(:,1)=='C');
Hother_idx=setdiff(setdiff(1:size(m,1),Hsnp_idx),Hcn_idx);

% break the other probes
% AFFX_???-5Q
%
% AFR_xxx
iafr=grep('AFR_...-',cellstr(m(Hother_idx,:)),1);
iafrsb=grep('AFR_..._SB',cellstr(m(Hother_idx,:)),1);
iafrnp=grep('AFR_..._NP',cellstr(m(Hother_idx,:)),1);

% AFFX_SNP
isnp=grep('AFFX-SNP',cellstr(m(Hother_idx,:)),1);

% RandomGC3 .. RandomGC25
ir=grep('Random',cellstr(m(Hother_idx,:)),1);

HsnpA_idx=Hsnp_idx(m(last_char_idx(Hsnp_idx))=='A');
HsnpB_idx=Hsnp_idx(m(last_char_idx(Hsnp_idx))=='B');
HcnA_idx=Hcn_idx(m(last_char_idx(Hcn_idx))=='A');
CN_no_AB=0;
if isempty(HcnA_idx)
  HcnA_idx=Hcn_idx;
  CN_no_AB=1;
end

if nnz(m_no_AB(HsnpA_idx,:)-m_no_AB(HsnpB_idx,:))>0
  error('SNP names of A and B do not match');
end


H1.adat=cat(3,H.dat(HsnpA_idx,:), H.dat(HsnpB_idx,:));
tmp=cat(3,H.dat(HcnA_idx,:),zeros(length(HcnA_idx),size(H.dat,2)));
H1.adat=cat(1,H1.adat,tmp);
v=zeros(size(H1.adat,1),1);
v(1:length(HsnpA_idx))=1;
H1=add_D_sup(H1,'SNPs','SNPs',v','rows');
v=zeros(size(H1.adat,1),1);
v(length(HsnpA_idx)+(1:length(HcnA_idx)))=1;
H1=add_D_sup(H1,'CNs','CNs',v','rows');
if CN_no_AB
  H1.marker=cellstr(strvcat(m_no_AB(HsnpA_idx,:),m(HcnA_idx,:)));
else
  H1.marker=cellstr(m_no_AB([ HsnpA_idx; HcnA_idx],:));
end  
H1.dat=sum(H1.adat,3);
H1.marker=strvcat(H1.marker);

v=str2int_matrix(H1.marker);
H1=add_D_sup(H1,'ID','ID',v','rows');


if ~exist('genome_info_file','var')
  H2=[];
  return
end

%% ----- load genome info file

%genome_info_file='/xchip/cancergenome/1chip/annotation/genome.info.6.0_build35';
f=fopen(genome_info_file);
GI=textscan(f,'%s%s%f','bufSize',50000000,'headerlines',0,'treatAsEmpty',{'NA'});
fclose(f);

% match SNP id
g=strvcat(GI{1});
Gsnp_idx=find(g(:,1)=='S' & ~isnan(GI{3}));
Gcn_idx=find(g(:,1)=='C' & ~isnan(GI{3}));

G_snp_num=str2int_matrix(g(Gsnp_idx,:));
[uu,ui,uj]=unique(G_snp_num);
qq=setdiff(1:size(G_snp_num,1),ui);

done=[];
for i=1:length(qq)
  a=find(G_snp_num==G_snp_num(qq(i)));
  a=setdiff(a,done);
  done=[done; a];
end

if (0)
  fid=fopen('/xchip/cancergenome/1chip/annotation/non-uniqe.txt','w');
  done=[];
  tmp=H1.gsupdat(3,find(H1.gsupdat(1,:)));
  ri=1;
  for i=1:length(qq)
    a=find(G_snp_num==G_snp_num(qq(i)));
    a=setdiff(a,done);
    done=[ done; a];
    if ~isempty(a)
      fprintf(fid,'-- %d ---\n', ri);
      ri=ri+1;
      t1=find(tmp==G_snp_num(qq(i)));
      if ~isempty(t1)
        fprintf(fid,'In 906,600\n');
      else
        fprintf(fid,'NOT in 906,600\n');
      end
      for j=1:length(a)
        fprintf(fid,'%s\n',[ GI{1}{Gsnp_idx(a(j))} ' ' GI{2}{Gsnp_idx(a(j))} ' ' num2str(GI{3}(Gsnp_idx(a(j))))]);
      end
    end
  end
  fclose(fid);
end

good_idx=setdiff(1:length(Gsnp_idx),done);
Ggood_snp_idx=Gsnp_idx(good_idx);
Ggood_snp_num=G_snp_num(good_idx);
% tmp=unique(Ggood_snp_num);


H1_snp_idx=find(H1.gsupdat(1,:));
[Mt,m1,m2]=match_num_sets(Ggood_snp_num,H1.gsupdat(3,H1_snp_idx));

%% match CN probes
G_cn_num=str2int_matrix(g(Gcn_idx,:));
[uc,uci,ucj]=unique(G_cn_num);
% copy number probes are unique!

if isunique(H1.gsupdat(3,find(H1.gsupdat(2,:))))~=1
  error('CN probes not unique');
end

H1_cn_idx=find(H1.gsupdat(2,:));
[Ct,c1,c2]=match_num_sets(G_cn_num,H1.gsupdat(3,H1_cn_idx));

H2=reorder_D_rows(H1,[ H1_snp_idx(m2) H1_cn_idx(c2)]);
H2=rmfield(H2,{'orig','history'});

gchr=strvcat(GI{2});

H2.chr=gchr([ Ggood_snp_idx(m1); Gcn_idx(c1)],:);
H2.pos=GI{3}([ Ggood_snp_idx(m1); Gcn_idx(c1)]);
H2=add_chrn(H2);

H2=order_by_pos(H2);
H2=rmfield(H2,{'orig','history'});
