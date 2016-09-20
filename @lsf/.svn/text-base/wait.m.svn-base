function [l,res]=wait(l,hvec,maxsec)
if (nargin==1) | isempty(hvec)
  hvec=find(l.jobs>0);
  verbose(['waiting for ' num2str(length(hvec)) ' jobs']);
end
if nargin<3
  maxsec=Inf;
end

start_hvec=hvec;
retries=zeros(length(start_hvec),1);
tic;
  nextdisp=10;
  while ~isempty(hvec) & (toc<maxsec)
    % Poll for results
    for i=1:length(hvec)
      basename=[ l.lsf_path get_basename(l,hvec(i))];
      if exist([basename '_finished.txt'],'file')
        disp(['collecting files for ' basename ]);
        pause(1);
        output_name=[basename '_out.mat'];
        if exist(output_name,'file')
          % load results
          l.res{hvec(i)}=load(output_name);
          
          % delete files
          
          check_file_exists_and_delete(output_name);
          check_file_exists_and_delete([basename '_in.mat']);
          check_file_exists_and_delete([basename '_code.m']);
          check_file_exists_and_delete([basename '_runme.bat']);
          check_file_exists_and_delete([basename '_stdout.txt']);
          check_file_exists_and_delete([basename '_stderr.txt']);
          check_file_exists_and_delete([basename '_code_main.c']);
          check_file_exists_and_delete([basename '_code_mcc_component_data.c']);
          check_file_exists_and_delete([basename '_code.ctf']);
          check_file_exists_and_delete([basename '_finished.txt']);
        end
        l.jobs(hvec(i))=0;
        hvec(i)=0;
      else
        % run again if stderr not empty  
        d=dir([basename '*']);
        erridx=findstrings(strvcat(d(:).name),[basename '_stderr.txt']);
        if ~isempty(erridx) && d(erridx).bytes>0
          % ignore error due to not finding shopt 
          [ux_status,ux_res]=unix(['grep -v shopt ' d(erridx).name ' | wc -l']);
          if str2num(ux_res)>0
            verbose(['[lsf::wait] stderr in ' basename ' - running again. Retry #' num2str(retries(hvec(i)))]);
            [ux_status1,ux_res1]=unix(['bkill ' num2str(l.jobs(hvec(i)))]);
%            check_file_exists_and_delete([basename '_stderr.txt']);
            if exist([basename '_stderr.txt'],'file');
              unix(['mv ' basename '_stderr.txt ' basename '_stderr.take' num2str(retries(hvec(i))) '.txt']);
            end
            l.jobs(hvec(i))=execute_lsf_job(l,basename);
            retries(hvec(i))=retries(hvec(i))+1;
          end
        end
        % run again if job disappeared (status: EXIT or DONE)
        [ux_status,ux_res]=unix(['bjobs ' num2str(l.jobs(hvec(i)))]);
        ux_lines=dlmsep(ux_res,char(10));
        if isempty(grep('(PEND|RUN)',ux_lines{2},1))
          verbose(ux_lines{2});
        end
        code_idx=grep('(EXIT|DONE)',ux_lines{2},1);
        if ~isempty(code_idx) && retries(hvec(i))<15
          verbose(['[lsf::wait] killing job ' num2str(l.jobs(hvec(i))) ' and restarting. Retry #' num2str(retries(hvec(i)))]);
          verbose(ux_lines{2});
          [ux_status1,ux_res1]=unix(['bkill ' num2str(l.jobs(hvec(i)))]);
          check_file_exists_and_delete([basename '_stderr.txt']);
          l.jobs(hvec(i))=execute_lsf_job(l,basename);
          retries(hvec(i))=retries(hvec(i))+1;
        end
      end 
    end
    if toc>nextdisp
      disp(['HVEC=' num2str(hvec) '; RETRIES=' num2str(retries') ]);
      lsf_update(l);
      nextdisp=nextdisp+10;
    end
    pause(0.1);
    hvec(hvec==0)=[];
  end
  check_file_exists_and_delete([l.lsf_path l.lsf_uid '.log']);
  res=l.res(start_hvec);

function check_file_exists_and_delete(fname)
cf_='[lsf::wait]';
if exist(fname,'file')
  verbose([cf_ 'deleting ' fname]);
  delete(fname);
end
