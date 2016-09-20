function C2=run_cbs(prefix,C,ch,should_log,should_smooth,window_size,rewrite_data_if_exists,write_header_for_output_file, localLSF)

nargs = nargin
if (nargs < 9)
    localLSF= 'LSF'
end

if should_log
  if(round(nanmean(reshape(C.dat,[],1)))==0)
    warning('data is already logratio! (mean=0)')
  end
  disp('taking log2(C)-1');
  C.dat(C.dat<0.1)=0.1;
  C.dat=log2(C.dat)-1;
end
if ~exist('rewrite_data_if_exists','var')
    rewrite_data_if_exists=1;
end

if isempty(ch)
    all_chrom=1;
    ch=1;
else
    all_chrom=0;
end

% write R script
orig_prefix=prefix;

if ~iscell(C.chr) %% make sure C.chr is a cell array
  C.chr = cellstr(C.chr);
end

for i=1:size(C.dat,2)
    %%%%%%% WRITE DATA to prefix.dat
    %  C.chrn=chromosome2num(C.chr);
    short_prefix=[ orig_prefix '.' num2str(i, ['%0' num2str(length(num2str(size(C.dat,2)))) 'd'])  ];
    if ~exist([short_prefix '.dat'],'file') || ...
            ( exist([short_prefix '.dat'],'file') && rewrite_data_if_exists )
        f=fopen([short_prefix '.dat'],'w');
        Cdat=C.dat(:,i);
        fprintf(f,'Marker\tChromosome\tPhysicalPosition\t%s\n', ...
                C.sdesc{i});
        %% Old behavior
        % for j=1:size(C.dat,1)
        %     fprintf(f,['%s %d %f' newline],C.chr{j},C.pos(j),Cdat(j));
        % end
        %%
        %% Write in chunks
        chunk=50000;
        for j=1:chunk:size(Cdat,1)
          cur_idx=j:min(j+chunk-1,size(Cdat,1));
          c_marker=as_row(C.marker(cur_idx));
          c_chr=as_row(C.chr(cur_idx));  
          c_pos=mat2cell(as_row(C.pos(cur_idx)),1,ones(length(cur_idx),1));
          c=[c_marker; c_chr; c_pos];
          c_dat=mat2cell(Cdat(cur_idx,1)',1,ones(length(cur_idx),1));
          cc=[c; c_dat];
          fprintf(f,'%s\t%s\t%d\t%0.3f\n',cc{:});
        end
        fclose(f);
    end
    for chi=ch
        prefix=[ orig_prefix '.' num2str(i, ['%0' num2str(length(num2str(size(C.dat,2)))) 'd']) '.' num2str(chi,'%02d') ];
        
        %%%%%%% WRITE R SCRIPT
        f=fopen([prefix '.R'],'w');
        fprintf(f,['library(DNAcopy)' newline]);
        fprintf(f,['genomdat <- read.table("' short_prefix ...
                   '.dat",header=T, row.names=1, check.names = FALSE, sep="\\t")' newline]);
        fprintf(f,['dim(genomdat)' newline]);
        if should_log
            log_st='data.type="logratio",';
        else
            log_st='';
        end
        fprintf(f,['CNA.object <- CNA(as.numeric(genomdat[,3]), ' ...
                   'genomdat[,1], as.numeric(genomdat[,2]),' log_st ...
                ' sampleid=colnames(genomdat)[3])' newline]);
        if ~all_chrom
            fprintf(f,['CNA.object <- subset.CNA(CNA.object, chromlist=' num2str(chi) ...
                    ')' newline]);
        end
        if should_smooth
            fprintf(f,['smoothed.CNA.object <- smooth.CNA(CNA.object)' newline]);
        else
            fprintf(f,['smoothed.CNA.object <- CNA.object' newline]);
        end    
        
        if ~isinf(window_size)
        else
          warning(['windowing is deprecated in DNAcopy >= v1.11. ' ...
                   'IGNORING value.']);
        end  
        fprintf(f,['segment.smoothed.CNA.object <- ' ...
                   'segment(smoothed.CNA.object' ...
                   ', min.width=5, verbose=2)' newline]);
        
        if write_header_for_output_file
          header='T';
        else
          header='F';
        end
        
        fprintf(f,['write.table(segment.smoothed.CNA.object$output,"' prefix ...
                     '.seg.dat",sep="\\t",row.names=F,col.names=' header ', quote=F)' ...
                     newline]);
        fprintf(f,['q()' newline]);
        fclose(f);
        if (strcmp(localLSF, 'LSF'))
             unix_str=['bsub ',...
                  '-E ~/CancerGenomeAnalysis/trunk/shell/chk_lsf_cga ',...
                  '-R "rusage[mem=4]" ',... 
                  '-mig 5 ',...
                  '-R "select[cpuf>100]" ',... 
                  '-Q "EXCLUDE(127)" ',...
                  '-q hour ',...
                  '-P cbs_cancerfolk ',...
                  '-o ' prefix '.out.txt -e ' prefix '.err.txt ',...
                  '-r ',... 
                  'R CMD BATCH --no-save --no-restore ' prefix '.R']; 
             %% The R version that is used above is the one that was defined
             %% by the 'use' command (dotkit) to setup the right $PATH
             disp(unix_str);
             [r1,r2]=unix(unix_str); 
             disp([strtrim(r2) ' Exit code: ' num2str(r1)])
        end
        if (strcmp(localLSF, 'local'))
             unix_str =['R CMD BATCH --no-save --no-restore ' prefix '.R'];
             disp(unix_str);
             [r1, r2] = unix(unix_str);
             disp([strtrim(r2) ' Exit code: ' num2str(r1)])
	     C2 = [prefix '.seg.dat'];
             return;
        end
    end 
end
C2=[];

% run script
% read output


