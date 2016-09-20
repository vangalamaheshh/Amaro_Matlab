function [X1]=printStruct(X,k,file,preheader)
% function [X1]=printStruct(X,k,file,preheader)
% print Structure in tab-delimited format
% inputs X: matlab structure with fields. fieldnames are printed in header.
%          field names 'N', 'header', 'headerline', and 'vlab' have special
%          meanings (field counters and header info) and will not appear in 
%          output file
%          field values of type numeric (vector) or cell (string vector)
%          are printed. Matrices, character arrays, or other more complex
%          structures do not print. 
%          fields only with value vector length same as the vector length 
%          of the last field are printed.  
%        k: indices of rows to print. k=[] prints all rows. 
%        file: output file. prints to screen if left empty
%        preheader: text to prepend to outputfile. 
% output: trimmed structure that was printed. 
%
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
%         preheader=char(preheader)';
        for i = 1:length(preheader)
           tmpLine = preheader{i};
           fprintf(fid,'%s\n',tmpLine);
        end
    else
        fprintf(fid,'%s',preheader);
        if (preheader(end)~=char(10))
            fprintf(fid,'\n');
        end
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

% if strfind(preheader,'columns')
%     printStructColumns(X,fx);
%     X1=X;
%     return
% end
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
    if n == 96
       disp(' ') 
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
                if c(1) == 'd'
                    fmt{n} = '%g';
                end
                if all(abs(x1)>1e-4)
                    fmt{n}='%0.6f';
                end
            end            
        case {'c'}
            fmt{n}='%s';
        otherwise
            fmt{n}='%.6g';
     
            
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


% function [X1]=printStructColumns(X,fx)
% 
% %disp(X)
% fmt={};
% Nf=length(fx);
% Nx=[];
% for n=1:Nf
%     f1=char(fx(n));
%     
%     x=X.(f1);
%     if ischar(x)
%         x=cellstr(x);
%         X.(f1)=x;
%     end
%     
%     Nx=[Nx, size(x,2)];
%     c=class(x);
%     switch c(1)
%         %case {'d','s','u','i'}
%         case {'d','u','i'}
%             fmt='%d';
%             x1=x; x1(isnan(x1))=[];
%             if any(round(x1)~=x1)
%                 % if %f = %d then use %f
%                 if all(abs(x1)>1e-4)                    
%                     fmt='%f'; 
%                 end
%             end   
%             for i=1:length(x)
%                 fprintf(fid,['\t' fmt],x(i))
%             end
%         case {'c'}
%             for i=1:length(x)
%                 fprintf(fid,['\t%s'],x{i})
%             end
%         otherwise
%             for i=1:length(x)
%                 fprintf(fid,['\t%.6g'],x(i))
%             end
%         fprintf(fid,'\n');
%     end
% end
%          
% function test()
% 
% clear
% cd /Users/stewardg/Projects/Spanner/src/matlab
% f='/Volumes/STEWARDG/share0/data/1000G/Pilot2/Spanner/build/NA12878/illumina/ERR001700_aligned.1.multi.span'
% X13=loadMultiSpan(f)
% k=1:10
% printStruct(X13,k)

