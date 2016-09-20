function endat=prepare_breakpointer_temp_files3XXX(drangerfile,bamfile,startfrom,endat,P)
% Yotam Drier, yotamd@gmail.com

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'quiet',1);
P=impose_default_value(P,'refdir','/xchip/tcga/gbm/analysis/lawrence/genome/hg18');
P=impose_default_value(P,'blacklist','none');
P=impose_default_value(P,'readlen',101); 
P=impose_default_value(P,'expand_pairs_extraction',500); 
P=impose_default_value(P,'max_insertion_size',450);
P=impose_default_value(P,'min_insertion_size',P.max_insertion_size);
P=impose_default_value(P,'expand_seq_around_bkpt',P.expand_pairs_extraction);
P=impose_default_value(P,'gsr_jar','/xchip/cga1/ydrier/CancerGenomeAnalysis/trunk/analysis_pipeline/tools/dist/GrabSplitReads.jar');
P=impose_default_value(P,'fastq_name','splitreads.fastq');
P=impose_default_value(P,'maxreads',1000000);
P=impose_default_value(P,'fish_allowedmm',80); %precentage
P=impose_default_value(P,'fish_tipsize',15);
P=impose_default_value(P,'fish_minmmintip',5);
P=impose_default_value(P,'fish_maxN',10);      %precentage

fprintf('min_insertion_size: %d, max_insertion_size: %d \n',P.min_insertion_size,P.max_insertion_size);
%if P.expand_seq_around_bkpt>P.expand_pairs_extraction+(P.max_insertion_size-P.min_insertion_size)
%    error('P.expand_seq_around_bkpt > P.expand_pairs_extraction + P.max_insertion_size');
%end

X = load_struct(drangerfile);
if slength(X)==0, error('dRanger results file empty'); end
if strncmpi(X.chr1{1},'chr',3), X.chr1 = convert_chrXXX(X.chr1,P); X.chr2 = convert_chrXXX(X.chr2,P); end
if strncmpi(X.str1{1},'(',1), ss={'(+)','(-)'}; X.str1=listmap(X.str1,ss)-1; X.str2=listmap(X.str2,ss)-1; end
if isfield(X,'num')
    X = make_numeric(X,{'num'});
else
    X.num=1:length(X.chr1);
end
X = make_numeric(X,{'tumreads','normreads','chr1','chr2','pos1','pos2','str1','str2'});
if ~exist('startfrom','var'), startfrom = 1; end
if ~exist('endat','var'), endat = slength(X); end
chrname=cell(1,2);

javaaddpath(P.gsr_jar);
%javaaddpath('/home/unix/stewart/CancerGenomeAnalysis/trunk/analysis_pipeline/tools/src/java');
import org.broadinstitute.cga.tools.*;
import java.lang.*;

%gsr = GrabSplitReads(String(bamfile),String(P.blacklist),P.maxreads,String(P.refdir));

try gsr = org.broadinstitute.cga.tools.seq.GrabSplitReads(String(bamfile),P.blacklist,P.maxreads,String(P.refdir));
catch gsr = org.broadinstitute.cga.tools.seq.GrabSplitReads(String(bamfile),P.blacklist,P.maxreads,String(P.refdir));
end

warnState = warning; %Save the current warning state
warning('off','Bioinfo:fastqwrite:AppendToFile');
for i=startfrom:endat
    dname = num2str(i);
    if ~exist(dname,'dir'), mkdir(dname); end
    fid=fopen([dname '/splitreads.helper2'],'w');
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',X.num(i),X.chr1(i),X.pos1(i),X.str1(i),...
        X.chr2(i),X.pos2(i),X.str2(i),X.tumreads(i),P.expand_seq_around_bkpt+P.readlen,P.expand_pairs_extraction,P.max_insertion_size);
    pos=[X.pos1(i) X.pos2(i)]-P.expand_seq_around_bkpt-P.readlen;
%    chrname(1)=index2name(X.chr1(i));
%    chrname(2)=index2name(X.chr2(i));
%    for k=1:2
%        seqf=fopen([P.refdir '/chr' char(chrname{k}) '.txt']);
%        fseek(seqf,pos(k),'bof');
%        seq=fread(seqf,2*P.expand_seq_around_bkpt+P.readlen,'*char');
%        fclose(seqf);
%        fprintf(fid,'%s\n',upper(seq));
%    end
fprintf(fid,'%s\n',upper(genome_regionXXX(X.chr1(i),X.pos1(i)-P.expand_seq_around_bkpt-P.readlen+1,X.pos1(i)+P.expand_seq_around_bkpt+P.readlen,P)));
fprintf(fid,'%s\n',upper(genome_regionXXX(X.chr2(i),X.pos2(i)-P.expand_seq_around_bkpt-P.readlen+1,X.pos2(i)+P.expand_seq_around_bkpt+P.readlen,P)));
    fclose(fid);
    if (X.str1(i)==0)
      shift_by_ins1=-1;
    else
      shift_by_ins1=1;
    end
    if (X.str2(i)==0)
      shift_by_ins2=-1;
    else
      shift_by_ins2=1;
    end
%    gsr.WriteSplitReads(X.chr1(i),X.pos1(i)+shift_by_ins1*P.max_insertion_size-P.expand_pairs_extraction,X.pos1(i)+shift_by_ins1*P.max_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.1']),P.fish_allowedmm,P.fish_tipsize,P.fish_minmmintip,P.fish_maxN);
%    fprintf('extract1: %d:%d-%d\n',X.chr1(i),X.pos1(i)+shift_by_ins1*P.max_insertion_size-P.expand_pairs_extraction,X.pos1(i)+shift_by_ins1*P.max_insertion_size+P.expand_pairs_extraction);
%    gsr.WriteSplitReads(X.chr2(i),X.pos2(i)+shift_by_ins2*P.max_insertion_size-P.expand_pairs_extraction,X.pos2(i)+shift_by_ins2*P.max_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.2']),P.fish_allowedmm,P.fish_tipsize,P.fish_minmmintip,P.fish_maxN);
%    fprintf('extract2: %d:%d-%d\n',X.chr2(i),X.pos2(i)+shift_by_ins2*P.max_insertion_size-P.expand_pairs_extraction,X.pos2(i)+shift_by_ins2*P.max_insertion_size+P.expand_pairs_extraction);
    if (shift_by_ins1==1)
        gsr.WriteUnmapped(X.chr1(i),X.pos1(i)+P.min_insertion_size-P.expand_pairs_extraction,X.pos1(i)+P.max_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.1']));
    else
        gsr.WriteUnmapped(X.chr1(i),X.pos1(i)-P.max_insertion_size-P.expand_pairs_extraction,X.pos1(i)-P.min_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.1']));
    end
    if (shift_by_ins2==1)
        gsr.WriteUnmapped(X.chr2(i),X.pos2(i)+P.min_insertion_size-P.expand_pairs_extraction,X.pos2(i)+P.max_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.2']));    
    else
        gsr.WriteUnmapped(X.chr2(i),X.pos2(i)-P.max_insertion_size-P.expand_pairs_extraction,X.pos2(i)-P.min_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.2']));
    end
    %gsr.WriteUnmapped(X.chr1(i),X.pos1(i)+shift_by_ins1*P.max_insertion_size-P.expand_pairs_extraction,X.pos1(i)+shift_by_ins1*P.max_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.1']));
    %gsr.WriteUnmapped(X.chr2(i),X.pos2(i)+shift_by_ins2*P.max_insertion_size-P.expand_pairs_extraction,X.pos2(i)+shift_by_ins2*P.max_insertion_size+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.2']));
    gsr.WriteColoredTips(X.chr1(i),X.pos1(i)-P.expand_pairs_extraction,X.pos1(i)+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.3']),P.fish_allowedmm,P.fish_tipsize,P.fish_minmmintip,P.fish_maxN);
    %fprintf('colored tips 1: %d:%d-%d\n',X.chr1(i),X.pos1(i)-P.expand_pairs_extraction,X.pos1(i)+P.expand_pairs_extraction);
    gsr.WriteColoredTips(X.chr2(i),X.pos2(i)-P.expand_pairs_extraction,X.pos2(i)+P.expand_pairs_extraction,String([dname '/' P.fastq_name '.4']),P.fish_allowedmm,P.fish_tipsize,P.fish_minmmintip,P.fish_maxN);
    %fprintf('colored tips 2: %d:%d-%d\n',X.chr2(i),X.pos2(i)-P.expand_pairs_extraction,X.pos2(i)+P.expand_pairs_extraction);
    FASTQStruct1 = fastqread([dname '/' P.fastq_name '.1']);
    FASTQStruct2 = fastqread([dname '/' P.fastq_name '.2']);
    FASTQStruct3 = fastqread([dname '/' P.fastq_name '.3']);
    FASTQStruct4 = fastqread([dname '/' P.fastq_name '.4']);
    t=[FASTQStruct1 FASTQStruct2 FASTQStruct3 FASTQStruct4];
    if (~isempty(t))
    [u,m]=unique(strcat({t.Header},{t.Sequence}));
    if exist([dname '/' P.fastq_name],'file')
        delete([dname '/' P.fastq_name]);
    end
    fastqwrite([dname '/' P.fastq_name], t(m));    
    end
end
warning(warnState) %Reset warning state to previous settings
gsr.closeFile();
fprintf('\n');
