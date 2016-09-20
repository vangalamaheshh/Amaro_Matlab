function prepare_breakpointer_temp_files(drangerfile,startfrom,endat,P)
% Yotam Drier 2010

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'quiet',1);
P=impose_default_value(P,'refdir','/xchip/tcga/gbm/analysis/lawrence/genome/hg18');
P=impose_default_value(P,'blacklist','none');
P=impose_default_value(P,'readlen',101);
P=impose_default_value(P,'expand_pairs_extraction',600);
P=impose_default_value(P,'expand_seq_around_bkpt',400);

if P.expand_seq_around_bkpt>P.expand_pairs_extraction
  error('P.expand_seq_around_bkpt>P.expand_pairs_extraction');
end

X = load_struct(drangerfile);
if slength(X)==0, error('dRanger results file empty'); end
if strncmpi(X.chr1{1},'chr',3), X.chr1 = convert_chr(X.chr1); X.chr2 = convert_chr(X.chr2); end
if strncmpi(X.str1{1},'(',1), ss={'(+)','(-)'}; X.str1=listmap(X.str1,ss)-1; X.str2=listmap(X.str2,ss)-1; end
X = make_numeric(X,{'tumreads','normreads','chr1','chr2','pos1','pos2','str1','str2','num'});
if ~exist('startfrom','var'), startfrom = 1; end
if ~exist('endat','var'), endat = slength(X); end
lf = P.expand_pairs_extraction+P.readlen;
rt = P.expand_pairs_extraction;
chrname=cell(1,2);
for i=startfrom:endat
    dname = num2str(i); 
    if ~exist(dname,'dir'), mkdir(dname); end
    fid=fopen([dname '/splitreads.helper'],'w');   
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',X.num(i),X.chr1(i),X.pos1(i),X.str1(i),...
        X.chr2(i),X.pos2(i),X.str2(i),X.tumreads(i),P.expand_seq_around_bkpt+P.readlen,lf,rt);
    pos=[X.pos1(i) X.pos2(i)]-P.expand_seq_around_bkpt-P.readlen; 
    chrname{1}=index2name(X.chr1(i));
    chrname{2}=index2name(X.chr2(i));
    for k=1:2
        seqf=fopen([P.refdir '/chr' chrname{k} '.txt']);
        fseek(seqf,pos(k),'bof');
        seq=fread(seqf,2*P.expand_seq_around_bkpt+P.readlen,'*char');
        fclose(seqf);
        fprintf(fid,'%s\n',upper(seq));    
    end    
    fclose(fid);    
end, fprintf('\n');