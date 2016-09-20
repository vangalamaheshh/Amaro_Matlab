function job=run_glad_batch(CN,rewrite_data,run_on_fly,lsf_queue,wait_for_jobs,gladidx,params,temp_dir)
%  job =
%  run_glad_batch(CN,rewrite_data,run_on_fly,lsf_queue,wait_for_jobs,gladidx,params,temp_dir)
%
%           glad_idx = indices of CN to glad
% 

%       Revisions:
%           12 Dec 07 -- changed location of R call to
%           ~/CancerGenomeAnalysis/trunk/bin
%
%---
% $Id$
% $Date: 2007-09-18 13:20:16 -0400 (Tue, 18 Sep 2007) $
% $LastChangedBy: rameen $
% $Rev$

job=[];
if ~exist('wait_for_jobs','var')
  wait_for_jobs=0;
end

if ~exist('lsf_queue','var')
  lsf_queue='classic';
end

if ~exist('rewrite_data','var')
  rewrite_data=1;
end

if ~exist('run_on_fly','var') || isempty(run_on_fly) 
  run_on_fly=1;
end

if ~exist('params','var')
  params=[];
end

if ~exist('temp_dir','var')
  temp_dir = './';
end

if ~exist('gladidx','var')
  gladidx=[1:size(CN.dat,2)];  %gladidx defaults to all columns
end

if iscell(CN.sdesc)
  CN.sdesc=strvcat(CN.sdesc);
end

if ~isempty(gladidx)
    for i=gladidx
        Ctmp = reorder_D_cols(CN,i,'allmem');
        fname=[temp_dir deblank(CN.sdesc(i,:)) '.dat'];
        % remove unwanted chars
        % fname=regexprep(fname,'[ \\\/#*%]','_');

        if exist(fname,'file')
            disp([ fname ' already exists']);
        end
        if (exist(fname,'file') && rewrite_data) || ...
                (~exist(fname,'file'))
            write_as_dchip(fname,Ctmp);
        end

        if exist([temp_dir 'Sample' sprintf('%03d',i) '.dat'],'file')
            unix(['rm' temp_dir 'Sample' sprintf('%03d',i) '.dat']);
        end
        unixstr=['ln -s ' replace_space_for_unix(fname) ' Sample' sprintf('%03d',i) '.dat'];
        disp(unixstr);
        unix(unixstr);

        run_glad([temp_dir 'Sample' sprintf('%03d',i)],Ctmp,0,params);
        if run_on_fly
            st=[temp_dir 'Sample' sprintf('%03d',i) ];
            if ~exist([st '.seg.dat'],'file')
                unixstr=['bsub -J GLAD1 -q ' lsf_queue ' -o ' st '.log.txt -e ' st '.err.txt  -i ' st ...
                    '.R -r ~/CancerGenomeAnalysis/trunk/bin/leadR211'];
                disp(unixstr);
                [r1,r2]=unix(unixstr);
                tok=regexp(r2,'\<([0-9]+)\>','tokens');
                if ~isempty(tok)
                    job(i)=str2num(tok{1}{1});
                else
                    job(i)=-1;
                end
            else
                disp(['skipping ' st]);
            end
        end
        deleteDfiles(Ctmp)
    end
end


if ~run_on_fly
  f=fopen('run_all','w');
  for i=1:getsize(CN.dat,2)
    st=[temp_dir 'Sample' sprintf('%03d',i) ];
    if ~exist([st '.seg.dat'],'file')
      fprintf(f,['bsub -J GLAD1 -q normal -o ' st '.log.txt -e ' st '.err.txt  -i ' st ...
                 '.R -r ~/CancerGenomeAnalysis/trunk/bin/leadR211\n']);
    else
      fprintf(f,['echo skipping ' st '\n' ]);
    end
  end
  fclose(f);
end

% unixstr='bsub -J A sleep 100';
% job=[];
% for i=1:10
%       [r1,r2]=unix(unixstr); 
%       tok=regexp(r2,'\<([0-9]+)\>','tokens');
%       if ~isempty(tok)
%         job(i)=str2num(tok{1}{1});
%       else
%         job(i)=-1;
%       end
% end

if wait_for_jobs
  while(~isempty(job))
    [r1,r2]=unix('bjobs');
    tok=regexp(r2,'\n([0-9]+)','tokens');
    running_jobs=str2num(strvcat(cat(1,tok{:})));
    fin=setdiff(job,running_jobs);
    if ~isempty(fin)
      disp('finished jobs:');
      disp(fin);
    end
    job=intersect(job,running_jobs);
    pause(0.1);
  end
end

%foreach fl( Sample*.seg.dat )
%   awk '{ print substr(FILENAME,7,3),$2,$3/1E6,$4/1E6,$5,$6 }' $fl | tail +2 > $fl:r:r.s.dat
%end
%cat *.s.dat > output

end
