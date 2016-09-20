function prepare_breakpointer_temp_files2(drangerfile,startfrom,endat,P)
% Yotam Drier 2010

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'quiet',1);
P=impose_default_value(P,'refdir','/xchip/tcga/gbm/analysis/lawrence/genome/hg18');
P=impose_default_value(P,'blacklist','none');
P=impose_default_value(P,'readlen',101);
P=impose_default_value(P,'expand_pairs_extraction',180);
P=impose_default_value(P,'max_insertion_size',450);
P=impose_default_value(P,'expand_seq_around_bkpt',400);

if P.expand_seq_around_bkpt>P.expand_pairs_extraction+P.max_insertion_size
  error('P.expand_seq_around_bkpt > P.expand_pairs_extraction + P.max_insertion_size');
end

X = load_struct(drangerfile);
if slength(X)==0, error('dRanger results file empty'); end
if strncmpi(X.chr1{1},'chr',3), X.chr1 = convert_chr(X.chr1); X.chr2 = convert_chr(X.chr2); end
if strncmpi(X.str1{1},'(',1), ss={'(+)','(-)'}; X.str1=listmap(X.str1,ss)-1; X.str2=listmap(X.str2,ss)-1; end
if exist('X.num','var')
    X = make_numeric(X,{'num'});
else
    X.num=1:length(X.chr1);
end
X = make_numeric(X,{'tumreads','normreads','chr1','chr2','pos1','pos2','str1','str2'});
if ~exist('startfrom','var'), startfrom = 1; end
if ~exist('endat','var'), endat = slength(X); end
chrname=cell(1,2);
for i=startfrom:endat
    dname = num2str(i); 
    if ~exist(dname,'dir'), mkdir(dname); end
    fid=fopen([dname '/splitreads.helper2'],'w');   
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',X.num(i),X.chr1(i),X.pos1(i),X.str1(i),...
        X.chr2(i),X.pos2(i),X.str2(i),X.tumreads(i),P.expand_seq_around_bkpt+P.readlen,P.expand_pairs_extraction,P.max_insertion_size);
    pos=[X.pos1(i) X.pos2(i)]-P.expand_seq_around_bkpt-P.readlen; 
    chrname(1)=index2name(X.chr1(i));
    chrname(2)=index2name(X.chr2(i));
    for k=1:2
        seqf=fopen([P.refdir '/chr' char(chrname{k}) '.txt']);
        fseek(seqf,pos(k),'bof');
        seq=fread(seqf,2*P.expand_seq_around_bkpt+P.readlen,'*char');
        fclose(seqf);
        fprintf(fid,'%s\n',upper(seq));    
    end    
    fclose(fid);    
end, fprintf('\n');