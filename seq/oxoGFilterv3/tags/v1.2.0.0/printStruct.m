function [X1]=printStruct(X,k,file,preheader)

f=fieldnames(X);

if ~isfield(X,'N')   
    X.N=length(X.(f{end}));
end

if (nargin<2)
    k=1:X.N;
else
    if length(k)>0
        X=trimStruct(X,k);
    end
end

if (nargin<3)
    fid=1;
else
    fid=fopen(file,'wt');
end



if (nargin>3)
    if iscell(preheader)
        preheader=char(preheader);
    end
    fprintf(fid,'%s',preheader);
    if (preheader(end)~=char(10))
        fprintf(fid,'\n');
    end
else
    preheader='';
end    

kmax=X.N;
fh={};
fx={};


for n=1:length(f)
    f1=char(f(n));
    x=X.(f1);
    m=size(x);
    if (m(2)>1), continue; end;
    % if ( (m(1)>0) && ~(strcmp(f1,'vlab')))
    %if ( (m(1)>1) && ~(strcmp(f1,'vlab')) && ~(strcmp(f1,'N')) && ~(strcmp(f1,'head'))&& ~(strcmp(f1,'hd')) )
    % if ( (m(1)==kmax) && ~(strcmp(f1,'vlab')) && ~(strcmp(f1,'head'))&& ~(strcmp(f1,'hd')) )
    if ( (m(1)==kmax) && ~(strcmp(f1,'vlab')) && ~(strcmp(f1,'N')) && ~(strcmp(f1,'head'))&& ~(strcmp(f1,'headline')) )
        if (m(1)<kmax)
            display('ooops')
            return;
        end
        fx{length(fx)+1}=f1;
    else
        fh{length(fh)+1}=f1;
    end
end

if strfind('columns',preheader)
    printStructColumns(X,fx);
    X1=X;
    return
end
%disp(X)
fmt={};
Nf=length(fx);
Nk=length(k);
Nx=[];
for n=1:Nf
    f1=char(fx(n));
    
    x=X.(f1);
    if ischar(x)
        x=cellstr(x);
        X.(f1)=x;
    end
    
    Nx=[Nx, size(x,2)];
    c=class(x);
    switch c(1)
        %case {'d','s','u','i'}
        case {'d','u','i'}
            
            fmt{n}='%d';
            x1=x; x1(isnan(x1))=[];
            if any(round(x1)~=x1)
                % if %f = %d then use %f
                if all(abs(x1)>1e-4)
                     % LTL: Changed from %f, so that we get the decimal points, instead of e.
                    fmt{n}='%.6f';
                end
            end            
        case {'c'}
            fmt{n}='%s';
        otherwise
            fmt{n}='%g';
     
            
    end
    for m=1:Nx(end)
        fprintf(fid,'%s',f1);
        if m<Nx(n)
            fprintf(fid,'\t');
        end
    end
    if n<Nf
        fprintf(fid,'\t');

    end
    
end
fprintf(fid,'\n');

for k1=1:X.N
    for n=1:Nf
        f1=char(fx(n));
        x=X.(f1);
        
        for m=1:Nx(n)
            if (isstruct(x))
                fprintf(fid,'*');
            elseif (iscell(x))
                c=char(x{k1});
                if (strcmpi(c,'NAN'))
                    c='';
                end
                fprintf(fid,fmt{n},c)    ;
            elseif (~iscellstr(x))
                if (~isnan(x(k1,m)))
                    fprintf(fid,fmt{n},x(k1,m));
                else
                    fprintf(fid,'%s','');
                end
            else
                fprintf(fid,fmt{n},x{k1});
            end
            if m<Nx(n)
                fprintf(fid,'\t');
            end
        end
        if n<Nf
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end

if (fid~=1)
    fclose(fid);
end

if (nargout>0)
    X1=X;
end


function [X1]=printStructColumns(X,fx)

%disp(X)
fmt={};
Nf=length(fx);
Nx=[];
for n=1:Nf
    f1=char(fx(n));
    
    x=X.(f1);
    if ischar(x)
        x=cellstr(x);
        X.(f1)=x;
    end
    
    Nx=[Nx, size(x,2)];
    c=class(x);
    switch c(1)
        %case {'d','s','u','i'}
        case {'d','u','i'}
            fmt='%d';
            x1=x; x1(isnan(x1))=[];
            if any(round(x1)~=x1)
                % if %f = %d then use %f
                if all(abs(x1)>1e-4)
                    % LTL: Changed from %f, so that we get the decimal points, instead of e.
                    fmt='%.6f'; 
                end
            end   
            for i=1:length(x)
                fprintf(fid,['\t' fmt],x(i))
            end
        case {'c'}
            for i=1:length(x)
                fprintf(fid,['\t%s'],x{i})
            end
        otherwise
            for i=1:length(x)
                fprintf(fid,['\t%g'],x(i))
            end
        fprintf(fid,'\n');
    end
end
         
function test()

clear
cd /Users/stewardg/Projects/Spanner/src/matlab
f='/Volumes/STEWARDG/share0/data/1000G/Pilot2/Spanner/build/NA12878/illumina/ERR001700_aligned.1.multi.span'
X13=loadMultiSpan(f)
k=1:10
printStruct(X13,k)

% load all retro fragments 
A=loadRetroSpan(f);
k=find(bitget(A.element,1)>0);
A1=trimSPanStruct(A,k);