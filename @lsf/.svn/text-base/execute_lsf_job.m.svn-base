function jobid=execute_lsf_job(l,basename)

if (l.use_compiler)
  [s,w]=unix(['source ~/.my.cshrc ; cd ' l.lsf_path '; bsub -r -o ' basename '_stdout.txt -e ' basename '_stderr.txt ' ...
              './' basename '_code']);  
else  
  [s,w]=unix(['cd ' l.lsf_path '; bsub -r -o ' basename '_stdout.txt -e ' basename '_stderr.txt ' ...
              'matlab -nodisplay -r ''' basename '_code; exit;''']);
end

pos1=find(w=='<');
pos2=find(w=='>');
if ~isempty(pos1) & ~isempty(pos2) & pos2(1)>pos1(1)
   jobid=str2num(w((pos1(1)+1):(pos2(1)-1)));
else
   jobid=-1;
end

