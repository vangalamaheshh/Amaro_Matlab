function tab=load_table(fname,dlm,headlines, getcols, maxlines, varargin)
%function tab=load_table(fname,dlm,headlines, getcols, maxlines, varargin)
% 
% returns matlab structure based on delimited tex file
% input:    fname:  file name
%           dlm: optional delimiter (tab = char(9))
%           headlines: default 1 header line with field names
%           getcols: specify indices of fields to return
%           maxlines (ignored): maximum number of lines to load 
%
% output:  structure with fields
%
set=exist('dlm','var');
if (set), set=~isempty(dlm); end
if (~set)
    dlm=char(9);
end
set=exist('headlines','var');
if (set), set=~isempty(headlines); end
if (~set)
    headlines=1;
    headlines=autodetectHeader(fname);
    if (headlines>1)
        fprintf('autodetect %d header lines\n',headlines);
    end
end
set=exist('getcols','var');
if (set), set=~isempty(getcols); end
if (~set)
    getcols=[];
end
set=exist('maxlines','var');
if (set), set=~isempty(maxlines); end
if (~set)
    maxlines=100000000;
else
   warning('maxlines is ignored, since it can cause out of memory errors.  Number of lines present is what will be loaded.') 
end

q=dir(fname);
if (numel(q)<1)
    error([' no file ' fname]);
end

tab={};
fid=fopen(fname,'r');
if (fid<1)
    error([' failed to open ' fname]);
end
% first line
nline=0;
if (headlines<1)
    if (iscell(getcols) )
        % nop header use column names from getcols - must be right
        tab.header=getcols;
        %
        tline = fgetl(fid);
        frewind(fid);
        if ~ischar(tline), return; end
        s=regexp(tline,dlm,'split');
        if (length(s)~=length(getcols))
            fprintf('%s\n',join('\t',getcols));
            fprintf('%s\n',tline);
            error(['getcols do not match columns in file ' fname '\n']);
            return;
        end
                
    else ( (headlines<1)&(~iscell(getcols)) )
        % no header, use v# column names 
        tline = fgetl(fid);
        frewind(fid);
        if ~ischar(tline), return; end
        v=regexp(tline,dlm,'split'); 
        NV=length(v);
        v=cellstr(sprintf('\tv%d',(1:NV))); v=regexp(v,'\t','split'); v=v{1}; v(1)=[];
        tab.header=v;            
    end
    headlines=0;
else    
    while (nline<headlines)
        tline = fgetl(fid);
        nline=1+nline;
        tab.headline{nline,1}=tline;
    end
    if ~ischar(tline), return; end
    tab.header=regexp(tline,dlm,'split');
    if (~isempty(getcols)&(iscell(getcols)) )
        %tline = fgetl(fid);
        
        [keep,km]=ismember(upper(tab.header),upper(getcols));
        %tab.header=regexp(tline,dlm,'split');
        if (sum(keep)<length(getcols))
            fprintf('%s\n',join('\t',getcols));
            fprintf('%s\n',char(strcat(tab.header))');
            error(['getcols mismatch file: ' fname '\n']);
            return;
        end    
        %tab.header=getcols;
       
    end
end
% don't allow blank in last header column
sl=cellfun(@length,tab.header);

if (sl(end)<1)
    tab.header=tab.header(1:(end-1));
    fprintf('%s', ['remove tab from end of header line: ' fname '\n']);
end

NV=length(tab.header);
format=[repmat('%s',1,NV)] ; % \n

fseek(fid,0,'bof');

tab.dat=textscan(fid,format,'HeaderLines',headlines,'Delimiter',dlm,'BufSize',1024*16-1); %,varargin{:});
fclose(fid);
NL=min(cellfun(@length,tab.dat((2:end)-1)));
% problem with last line skipping if no end-of-line
% C = textscan(fid, '%s', 'delimiter', sprintf('\n'));C=C{1};C=regexp(C,dlm,'split');
% fix missing CR on last line...
NLend=min(cellfun(@length,tab.dat(end)));
if (NLend==(NL-1))
    tab.dat{end}=[tab.dat{end}; {''}];
end
    
if (~isempty(getcols)&(~iscell(getcols)) )
    tab.header=tab.header(getcols);
    tab.dat=tab.dat(:,getcols);
    NV=length(tab.header);
end

if (~isempty(getcols)&(iscell(getcols)) )    
    [keep,km]=ismember(upper(tab.header),upper(getcols));            
    tab.header=tab.header(keep);
    tab.dat=tab.dat(:,keep);
    NV=length(tab.header);
end


for i=1:NV
    v1=tab.dat{i};
    tab.dat{i}=NaN;
    v1=v1(1:NL);
    s1=tab.header{i};
    s1=strtrim(s1);
    s1=regexprep(s1,' ','_');
    s1=regexprep(s1,'(','_');
    s1=regexprep(s1,')','_');
    s1=regexprep(s1,'/','_');
    s1=regexprep(s1,'#','');
    s1=regexprep(s1,'%','');
    s1=regexprep(s1,'&','');
    s1=regexprep(s1,'>','_');
    s1=regexprep(s1,'-','_');
    s1=regexprep(s1,'\.','_');
    s1=regexprep(s1,'__','_');
    if (length(s1)<1), continue; end
    if (~isletter(s1(1)))
        s1=['v_' s1];
    end
    
    
    abc=any(cellfun(@length,regexp(regexprep(v1,'NaN','0'),'[b-c|f-m|o-z|B-C|F-M|O-Z]','start'))>0);
     
    
    if (~abc)
        v1(strmatch('NA',v1,'exact'))={'NaN'};
        v1(strmatch('---',v1,'exact'))={'NaN'};
        v1(strmatch('.',v1,'exact'))={'0'};

        k=find(cellfun(@length,v1)<1);
        if (~isempty(k)) 
            v1(k)={'NaN'};
        end
        
        if (NL<100000)
            x1=str2num(char(v1));
        else
            x1=zeros(size(v1));
            for i=1:NL
                x1(i)=str2double(v1{i});
            end
        end
        if (length(x1)==NL)
            v1=x1;
            x1=[];
        end
    end
    tab.(s1)=v1;
end
tab=rmfield(tab,'dat');


function headlines=autodetectHeader(fname)
commentchars={'#','%'};
fid=fopen(fname,'r');
if (fid<1)
    error([' failed to open ' fname]);
end
headlines=1;
tline = fgetl(fid);
while ischar(tline)
    if ~strcmp(tline(1),commentchars)
        fclose(fid);
        return
    end
    headlines=headlines+1;
    tline = fgetl(fid);
end
fclose(fid);




function test
unix('grep -100 Hugo /local/cga-fh/cga/An_ESO/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/mut/aggregated/PR_Esophageal_CIP_WGS.maf.annotated > ~/test.maf')
X=load_table('~/test.maf',[],[], [1:17 64:65 81:86])



f='/local/cga-fh/cga/An_ESO/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/mut/aggregated/PR_Esophageal_CIP_WGS.maf.annotated'
f='/local/cga-fh/cga/An_ESO/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/mut/aggregated/PR_Esophageal_CIP_WGS.maf.annotated'


clear
cd ~/Projects/Cancer/DoubleNormal/FH
f1='~/Projects/Cancer/DoubleNormal/FH/LUSC_Samples.tsv'
f1='~/Projects/Cancer/DoubleNormal/FH/test.tsv'
S1=load_table(f1);

f1='/Users/stewart/Projects/Cancer/ExonCapture/refGene.hg19.20100825.sorted.header.txt'
G=load_table(f1,char(9),1);

f2='/Users/stewart/Projects/Cancer/ExonCapture/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list.notabtab.txt'
T=load_table(f2,char(9),0,{'chr','p1','p2','strand','targ'});
    
f3='/Users/stewart/Projects/Cancer/ExonCapture/UCSC.HG19.wgEncodeCrgMapabilityAlign40mer.bwig'
M40=load_table(f3,char(9),1,{'chr','p1','p2','map'});
    M40=rmfield(M40,{'dat'});    
    M40.a=chrom2num(char(M40.chr));

    

f='/Users/stewart/Projects/Cancer/Pediatric/EwingsSarcoma/covbed/EwngSRC-SJDES001.bed'
t=load_table(f,char(9),-1,{'chr','p1','p2','str','targ','nt','nn'});

f='/Users/stewart/Projects/Cancer/Pediatric/Rhabdoid/PR_SIGMA_Rhabdoid.review.clinicalinfo.absolute.maf.txt'
t=load_table(f,char(9));


fl{ln}=dlmsep(tline,dlm);
if ischar(fname)
    fid=fopen(fname,'r');
else
    fid=fname;
end

ln=1;
if ~exist('nlines','var')
    nlines=Inf;
    do_close=1;
else
    do_close=0;
end

had_output=0;
fl={};
while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fl{ln}=dlmsep(tline,dlm);
    ln=ln+1;
    if mod(ln,1000)==0
        verbose(['...' num2str(ln)],30);
        had_output=1;
    end
end
if do_close
    fclose(fid);
    fid=-1;
end
ln=ln-1;
if had_output
    verbose([newline],30);
end



fpos=0;
if headerlines~=0
    if ischar(fname)
        f=fopen(fname,'r');
    else
        f=fname;
        fpos=ftell(f);
    end
    if length(dlm)~=1
        tab.dlm=find_dlm(fname,dlm);
    else
        tab.dlm=dlm;
    end
    if headerlines>0
        tab.headers=read_dlm_file(f,dlm,headerlines);
    elseif headerlines==-1 % R convention
        headerlines=1;
        tab.headers=read_dlm_file(f,dlm,headerlines);
        tab.headers{1}=['EMPTY' tab.headers{1,:}];
    end
    %   fclose(f);
else
    if ischar(fname)
        f=fopen(fname,'r');
    else
        f=fname;
        fpos=ftell(f);
    end
    tab.headers={};
    tab.dlm=dlm;
end

verbose(['Reading file using format:' format],10);
fseek(f,fpos,'bof');
tab.dat=textscan(f,format,'headerLines',headerlines,'delimiter',tab.dlm,varargin{:});
fclose(f);