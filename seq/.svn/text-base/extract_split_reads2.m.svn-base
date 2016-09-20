function extract_split_reads2(bamname,drangerfile,startfrom,endat,P)
% Yotam Drier 2010

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'quiet',1);
P=impose_default_value(P,'refdir','*required*');
P=impose_default_value(P,'blacklist','none');
P=impose_default_value(P,'readlen',101);
P=impose_default_value(P,'expand_pairs_extraction',600);
P=impose_default_value(P,'expand_seq_around_bkpt',400);
P=impose_default_value(P,'bamgrasp_jar','*required*');
P=impose_default_value(P,'maxreads',30000);

if P.expand_seq_around_bkpt>P.expand_pairs_extraction
  error('P.expand_seq_around_bkpt>P.expand_pairs_extraction');
end

X = load_struct(drangerfile);
if slength(X)==0, error('dRanger results file empty'); end
if strncmpi(X.chr1{1},'chr',3), X.chr1 = convert_chr(X.chr1); X.chr2 = convert_chr(X.chr2); end
if strncmpi(X.str1{1},'(',1), ss={'(+)','(-)'}; X.str1=listmap(X.str1,ss)-1; X.str2=listmap(X.str2,ss)-1; end
X = make_numeric(X,{'tumreads','normreads','chr1','chr2','pos1','pos2','str1','str2','num'});

javaclasspath(P.bamgrasp_jar);
import org.broadinstitute.cga.tools.bamgrasp.*;
import java.lang.*;
try x = org.broadinstitute.cga.tools.bamgrasp.BamGrasp();
catch x = BamGrasp();
end
if P.quiet, x.setQuietModeOn(); end
x.set_maxReads(P.maxreads) 
x.openFile(String(bamname),String(P.blacklist),String(P.refdir));

if ~exist('startfrom','var'), startfrom = 1; end
if ~exist('endat','var'), endat = slength(X); end

lf = P.expand_pairs_extraction+P.readlen;
rt = P.expand_pairs_extraction;
trim = P.expand_pairs_extraction-P.expand_seq_around_bkpt;
for i=startfrom:endat
    if ~mod(i,10), fprintf('%d ',i); end
    [R1,B1,S1,tmp1,tmp2,tmp3,rname1] = BamGrasp_load_region(x,X.chr1(i),max(X.pos1(i)-lf,1),X.pos1(i)+rt,P);
    [R2,B2,S2,tmp1,tmp2,tmp3,rname2] = BamGrasp_load_region(x,X.chr2(i),max(X.pos2(i)-lf,1),X.pos2(i)+rt,P);

    unmapped1=find(R1(:,4)==-200);
    unmapped2=find(R2(:,4)==-200);
    ind1=getcoloredtips(R1,B1);
    ind2=getcoloredtips(R2,B2);

    dname = num2str(i);
    if ~exist(dname,'dir'), mkdir(dname); end
    fastq=fopen([dname '/splitreads.fastq'],'w');
    dump2fastq(fastq,R1(unmapped1,:),B1,rname1(unmapped1));
    dump2fastq(fastq,R2(unmapped2,:),B2,rname2(unmapped2));
    dump2fastq(fastq,R1(ind1,:),B1,rname1(ind1));
    dump2fastq(fastq,R2(ind2,:),B2,rname2(ind2));
    fclose(fastq);

    fid=fopen([dname '/splitreads.helper'],'w');   
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',X.num(i),X.chr1(i),X.pos1(i),X.str1(i),...
        X.chr2(i),X.pos2(i),X.str2(i),X.tumreads(i),P.expand_seq_around_bkpt+P.readlen);
%    fprintf(fid,'%s\n%s\n',S1(trim+2:end-trim),S2(trim+2:end-trim));  (gave coordinates one base too small)
  fprintf(fid,'%s\n%s\n',S1(trim+1:end-trim),S2(trim+1:end-trim));

    fclose(fid);
end, fprintf('\n');

x.closeFile();

