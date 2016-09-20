function write_mit_cls_file(fname,D,tps,use_multiclass,start_at_zero)
% WRITE_MIT_CLS_FILE writes an cls file
%    WRITE_MIT_CLS_FILE(FNAME,D,TPS,USE_MULTICLASS,START_AT_ZERO)


if ~exist('use_multiclass','var')
    use_multiclass=0;
end

if ~exist('start_at_zero','var')
    start_at_zero=0;
end

if use_multiclass
    if length(tps)>1
        error('too many supids');
    else
        n_cls=D_n_cls(D,tps);
        [typeacc,typedesc,D1,range,non_empty]=decollapse_supdat(D,tps,1:n_cls);
        write_mit_cls_file(fname,D1,range,0,start_at_zero);        
    end
else        
    f=fopen(fname,'w');
    fprintf(f,'%d %d 1\r\n',size(D.dat,2),max(length(tps),2));
    for i=1:length(tps)
        st=deblank(D.supacc(tps(i),:));
        st(st==' ')='_';
        sacc{i}=st;
    end
    if length(tps)==1
%        sacc={ ['not-' sacc{1}], sacc{:} };
        sacc={ ['not_' sacc{1}], sacc{:} };
    end
    
    s=D.supdat(tps,:);
    
    if any(isnan(s(:)))
        disp('Outputting NaNs as null');
    end
    %s(isnan(s))=0;
    if (max(nansum(s,1))>1) %nan
        error('A sample belongs to more than one group');
    end
    if size(s,1)>1
        [m,ind]=nanmax(s,[],1); %nan
    else
        ind=s+1;
    end
    
    % newlabel=zeros(max(ind),1);
    % curl=1;
    % for i=1:length(ind)
    %     if newlabel(ind(i))==0
    %         newlabel(ind(i))=curl;
    %         curl=curl+1;
    %     end
    % end
    
    
    fprintf(f,'#');
    for i=1:length(sacc)
        fprintf(f,' %s',sacc{i});
    end
    fprintf(f,'\r\n');
    for i=1:length(ind)
        if isnan(ind(i))
            fprintf(f,'null ');
        else
            if length(tps)>1
              if start_at_zero
                fprintf(f,'%d ',ind(i)-1);
              else
                fprintf(f,'%d ',ind(i));
              end                
            else
                fprintf(f,'%d ',ind(i)-1);
            end
        end
    end
    fclose(f);
end

