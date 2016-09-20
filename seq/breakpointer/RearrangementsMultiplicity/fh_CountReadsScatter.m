function fh_CountReadsScatter(sample,dRangerFile,BreakPointerFile,bam,blacklist,refdir,insertionsize,fish_expand_pairs_extraction,fish_maxreads,readlen,align_min_qual,align_minscorediff,libdir,firstidx,lastidx,firstsample,lastsample)
% Yotam Drier, yotamd@gmail.com


if nargin~=17, error('Usage: fh_CountReadsScatter(sample,dRangerFile,BreakPointerFile,BAM file(s),blacklist,refdir,insertionsize,fish_expand_pairs_extraction,fish_maxreads,readlen,align_min_qual,align_minscorediff,libdir,firstidx,lastidx,firstsample,lastsample)'); end

if strcmpi(sample,'NO_DATA')
  fprintf('fh_CountReadsScatter exiting because received "NO_DATA" signal\n');
  return;
end

demand_file(dRangerFile); demand_file(BreakPointerFile); demand_file(bam);
if ~strcmp(blacklist,'none'), demand_file(blacklist); end
if ~isnumeric(readlen), readlen=str2double(readlen); end
if ~isnumeric(fish_expand_pairs_extraction), fish_expand_pairs_extraction=str2double(fish_expand_pairs_extraction); end
if ~isnumeric(fish_maxreads), fish_maxreads=str2double(fish_maxreads); end
if ~isnumeric(align_min_qual), align_min_qual=str2double(align_min_qual); end
if ~isnumeric(align_minscorediff), align_minscorediff=str2double(align_minscorediff); end
if ~isnumeric(insertionsize)
    if isempty(strfind(insertionsize,'.isz'))
        insertionsize = str2double(insertionsize);
    else
        demand_file(insertionsize);
        inshist = dlmread(insertionsize, '\t');
        [m, i]=max(sum(inshist(:,2:end),2));
        insertionsize=inshist(i,1);
    end
end
if ~isdir(refdir), error('Not found: reference directory %s',refdir); end
if ~isnumeric(firstidx), firstidx=str2double(firstidx); end
if ~isnumeric(lastidx), lastidx=str2double(lastidx); end
if ~isnumeric(firstsample), firstsample=str2double(firstsample); end
if ~isnumeric(lastsample), lastsample=str2double(lastsample); end

if firstidx<1 || ((lastidx~=-1) && (lastidx<1 || firstidx>lastidx)), error('problem with firstidx,lastidx'); end
if firstsample<1 || ((lastsample~=-1) && (lastsample<1 || firstsample>lastsample)), error('problem with firstsample,lastsample'); end

%javaaddpath('/xchip/cga1/ydrier/CancerGenomeAnalysis/trunk/analysis_pipeline/tools/dist/ReferenceInfo.jar')
%ReferenceInfoObj.init(refdir)
P=[];
P.joindump_cmd = fullfile(libdir,'joindump.pl'); demand_file(P.joindump_cmd);
P.gsr_jar = fullfile(libdir,'GrabSplitReads.jar'); demand_file(P.gsr_jar);
P.refdir = refdir;
P.blacklist = blacklist;
P.maxreads = fish_maxreads;
P.readlen = readlen;
P.pairs_extraction=insertionsize+fish_expand_pairs_extraction;
P.align_min_qual=align_min_qual;
P.minscorediff=align_minscorediff;
outfn=[sample '.support.txt'];
if ((firstidx > lastidx)&&(lastidx > -1))||((firstsample > lastsample)&&(lastsample > -1))
    fid = fopen(outfn,'w');
    fprintf(fid,'NO_DATA\n');
    fclose(fid);
else
    if strcmp(bam(end-2:end),'bam')
        if exist(outfn,'file')
            disp([outfn ' already exists, skipping.']);
        else
            support=countreads(dRangerFile,BreakPointerFile,bam,firstidx,lastidx,P); 
            save_struct(support, outfn);
        end
    else if (strcmp(bam(end-2:end),'txt')||strcmp(bam(end-3:end),'list'))
            bamfiles=load_struct(bam);
            if isfield(bamfiles,'cluster')
              bamfiles=make_numeric(bamfiles,'cluster');
              cluster=bamfiles.cluster(strcmp(sample,bamfiles.individual));
              inc=bamfiles.cluster==cluster;
              if ~any(inc)
                  error('No other normal samples in the same cluster were found');
              else
                  bamfiles=structfun(@(x)x(inc),bamfiles,'UniformOutput',false);
              end
            end
%             mkdir(bamfiles.sample{1});
%             cd(bamfiles.sample{1});
%             support=countreads(dRangerFile,BreakPointerFile,bamfiles.bam_file{1},firstidx,lastidx,P); 
%             save_struct(support, 'support.onenormal.txt');    
%             flds=setdiff(fieldnames(support),{'num','drnum','bpnum'});
%             l=slength(bamfiles);
%             for i=2:l
%                 mkdir(['../bam_' num2str(i)]);
%                 cd(['../bam_' num2str(i)]);            
%                 support1=countreads(dRangerFile,BreakPointerFile,bamfiles.fn{i},firstidx,lastidx,P);                
%                 save_struct(support1, 'support.onenormal.txt');    
%                 for j=1:length(flds)
%                     support.(flds{j})=support.(flds{j})+support1.(flds{j});
%                 end
%             end 
%             for j=1:length(flds)
%                 support.(flds{j})=support.(flds{j})/l;
%             end
%            cd ..    
            l=slength(bamfiles);
            if (lastsample > l) || (lastsample == -1)
                lastsample=l;
            end
            for i=firstsample:lastsample
                %if exist(bamfiles.individual{i},'dir')
                %    disp([bamfiles.individual{i} ' dir already exists for sample ' sample ', skipping.']);
                %else
                    mkdir(bamfiles.individual{i});
                    cd(bamfiles.individual{i});
                    %support=countreads(dRangerFile,BreakPointerFile,bamfiles.bam_file{i},firstidx,lastidx,P);
                    %save_struct(support, 'support.onenormal.txt');
                    fid = fopen('readstotnum.onenormal.txt','w');
                    fprintf(fid,bamfiles.totreads{i});
                    fclose(fid);
                    cd ..
                %end
            end
        else
            error('BAM file(s) parameter must be one of .bam; .txt; or .list'); 
        end
    end
end

