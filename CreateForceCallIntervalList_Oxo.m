function CreateForceCallIntervalList_Oxo(infile1,infile2,outfile,mouse)
% a function to take 2 callstats, maflite, or MAF files and find the union of all
% calls in order to output a list of intervals for re-running the mutation
% caller in --force-output mode at those sites.
% can mix&match file types, as long as the contigs are consistent (i.e., 
% X vs X, not X vs 23...
% for Callstats, currently does not discriminate for "KEEP" rows.
% By Peleg Horowitz, 2012

if ~exist(infile1, 'file')
    error([infile1 ' File not found!'])
end

if ~exist(infile2, 'file')
    error([infile2 ' File not found!'])
end

%[DIR, FNAME, EXT] = fileparts(outfile);
%DIR = add_slash_if_needed(DIR);

%if ~exist(DIR, 'dir')
%    error([DIR ' Directory not found!'])
%end
DIR=pwd;

[s1, w1] = unix(['head -n 1 ' infile1]);
[s2, w2] = unix(['head -n 3 ' infile2]);

if s1==0 && s2==0
is1callstat = strcmp('## muTector',w1(1:11)); %use columns 1, 2
is2callstat = strcmp('## muTector',w2(1:11));
is1maflite = strcmp('build',w1(1:5)); %no header row, use columns 2, 3
is2maflite = strcmp('build',w2(1:5));
is1annotatedmaf = regexp(w1,'Hugo_Symbol'); %use columns 5, 6
is2annotatedmaf = ~isempty(regexp(w2,'## Oncotator'));
else
    error('Something went wrong with opening either infile')
end
%%%%      is1annotatedmaf=1; is2annotatedmaf=2;
if ~is1callstat && ~is1maflite && ~is1annotatedmaf
    error(['Format not recognized for file ' infile1])
end

if ~is2callstat && ~is2maflite && ~is2annotatedmaf
    error(['Format not recognized for file ' infile2])
end

if is1callstat
    % pipe1pre = ' | grep KEEP '
    chrcol1 = '1';
    poscol1 = '2';
    pipe1 = ' | grep -v t';
elseif is1maflite
    chrcol1 = '2';
    poscol1 = '3';
    pipe1 = ' | grep -v chr';
elseif is1annotatedmaf
    chrcol1 = '5';
    poscol1 = '6';
    pipe1 = ' | grep -v o | grep -v GL';
end

if is2callstat
    % pipe2pre = ' | grep KEEP '
    chrcol2 = '1';
    poscol2 = '2';
    pipe2 = ' | grep -v t';
elseif is2maflite
    chrcol2 = '2';
    poscol2 = '3';
    pipe2 = ' | grep -v chr';
elseif is2annotatedmaf
    chrcol2 = '5';
    poscol2 = '6';
    pipe2 = ' | grep -v o | grep -v \# | grep -v GL';
end
if ~mouse
unix(['cat ' infile1 ' | cut -f' chrcol1 ',' poscol1 pipe1 ' > ' DIR '/tmp.txt']);
unix(['cat ' infile2 ' | cut -f' chrcol2 ',' poscol2 pipe2 ' >> ' DIR '/tmp.txt']);
end
if mouse
    unix(['cat ' infile1 ' | cut -f' chrcol1 ',' poscol1 pipe1 ' > ' DIR '/tmp.txt']);
    unix(['cat ' infile2 ' | cut -f' chrcol2 ',' poscol2 pipe2 ' >> ' DIR '/tmp.txt']);
    %unix(['sed -i s/^/chr/ ' DIR '/tmp.txt']);
end
s=unix(['sed s/\t/\:/ ' DIR '/tmp.txt > ' outfile]);
if s==0
    unix(['rm -f ' DIR '/tmp.txt']);
end

%% more complex stuff: finding unique sites

unix(['cat ' outfile ' | grep -v X | grep -v Y | grep -v MT | grep -v M > ' DIR '/tmp.txt']);

[sX, w] = unix(['grep -c X ' outfile]);
if sX==0
    unix(['cat ' outfile ' | grep X | cut -f2 > ' DIR '/tmpX.txt']);
end

%[sY, w] = unix(['grep -c Y ' outfile]);
%if sY==0
%    unix(['cat ' outfile ' | grep Y | cut -f2 > ' DIR '/tmpY.txt']);
%end

%[sMT, w] = unix(['grep -c MT ' outfile]);
%if sMT==0
%    unix(['cat ' outfile ' | grep MT | cut -f2 > ' DIR '/tmpMT.txt']);
%end

unix(['rm -f ' outfile]);

A = load([DIR '/tmp.txt']);
An = A(:,1)*1000000000+A(:,2);
An = unique(An);
A = nan(length(An),2);
A(:,1) = floor(An/1000000000);
A(:,2) = mod(An,1000000000);

FID = fopen(outfile,'w');
for k=1:length(A)
    fprintf(FID,'%s%s%s\n',num2str(A(k,1)),':',num2str(A(k,2)));
end
unix(['rm -f ' DIR '/tmp.txt']);

if sX==0;
    B = load([DIR '/tmpX.txt']);
    B = unique(B);
    for k=1:length(B)
        fprintf(FID,'%s%s\n','X:',num2str(B(k)));
    end
    unix(['rm -f ' DIR '/tmpX.txt']);
end

%if sY==0
%    C = load([DIR '/tmpY.txt']);
%    C = unique(C);
%    for k=1:length(C)
%        fprintf(FID,'%s%s\n','Y:',num2str(A(k,2)));
%    end
%    unix(['rm -f ' DIR '/tmpY.txt']);
%end
%if sMT==0
%    D = load([DIR '/tmpMT.txt']);
%    D = unique(D);
%    for k=1:length(D)
%        fprintf(FID,'%s%s\n','MT:',num2str(A(k,2)));
%    end
%    unix(['rm -f ' DIR '/tmpMT.txt']);
%end
fclose(FID);
if mouse
   unix(['sed -i s/^/chr/ ' DIR '/' outfile]);
end

