function dRanger_manual_review(samples,P)
% dRanger_manual_review(samples,P)
%
% for each specified sample:
% (1) loads "results" input file and extracts records for which filter==0
% (2) writes BED file with two lines per candidate rearrangement (one for each breakpoint):
%     chr start end name
% (3) prompts user to open an IGV session with the tumor+normal BAMS and the review BED file
% (4) allows interactive traversal of the list
%     and recording of user's judgements for each breakpoint (0 = keep, 1,2,3... = filter out)
% (5) when user is finished, outputs "reviewed_results" file.
%
% Mike Lawrence 2009-07-20

if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'results_name','dRanger_results');
P=impose_default_value(P,'review_bedfile_name',[]);
P=impose_default_value(P,'review_somatic_only',true);
if isfield(P,'review_somatics_only'), P = rename_field(P,'review_somatics_only','review_somatic_only'); end
P=impose_default_value(P,'use_IGV_control',true);
P=impose_default_value(P,'reviewing_workstation_IP',[]);
P=impose_default_value(P,'IGV_listening_port',60151);
P=impose_default_value(P,'IGV_window_width',2000);

import java.util.*;
import java.io.*;
import java.net.Socket;

try

fprintf('\ndRanger_manual_review\n');
if P.review_somatic_only, fprintf('  (somatic only)\n'); end

for s=1:length(samples)
  sample = samples{s};

  direc = ['/xchip/tcga_scratch/lawrence/' sample];
  if isdir(direc)
    tbam = [direc '/tumor.bam'];
    nbam = [direc '/normal.bam'];
    resultsname = [direc '/' P.results_name '.txt'];
    if isempty(P.review_bedfile_name), review_bedfile_name =...
          regexprep([P.results_name '_review.bed'],'dRanger','dR'); % so IGV won't think it's a dRanger results file
    else review_bedfile_name = P.review_bedfile_name;
    end
    rbfname = [direc '/' review_bedfile_name];
    rrfname = [direc '/' P.results_name '.txt'];
  else
    resultsname = sample;
    tbam = '<tumor_bamfile>';
    nbam = '<normal_bamfile>';
    rbfname = [resultsname '.review.bed'];
    rrfname = [resultsname '.review.txt'];
  end

  % (1) load results file and extract records that survived automatic filtering
  fprintf('Loading %s...\n',resultsname);
  if ~exist(resultsname,'file'), error('Can''t find file %s',resultsname); end
  X = load_struct(resultsname);
  if slength(X)==0, error('dRanger results file empty'); end
  if strncmpi(X.chr1{1},'chr',3), X.chr1 = convert_chr(X.chr1); X.chr2 = convert_chr(X.chr2); end
  if strncmpi(X.str1{1},'(',1), ss={'(+)','(-)'}; X.str1=listmap(X.str1,ss)-1; X.str2=listmap(X.str2,ss)-1; end
  X = make_numeric(X,{'chr1','str1','chr2','str2','tumreads','normreads'});
  if isfield(X,'filter')
    X = make_numeric(X,'filter');
  else
    X.filter = zeros(slength(X),1);
  end
  if isfield(X,'min1')
    X = make_numeric(X,{'min1','min2','max1','max2'});
  else
    X = make_numeric(X,{'pos1','pos2'});
    X.min1 = X.pos1 - 200;
    X.max1 = X.pos1 + 200;
    X.min2 = X.pos2 - 200;
    X.max2 = X.pos2 + 200;
  end
  Xall = X;
  nxall = slength(Xall);
  if P.review_somatic_only, fidx = find(Xall.filter==0 & Xall.normreads==0);
  else fidx = find(Xall.filter==0); end
  X = reorder_struct(Xall,fidx);
  nx = slength(X);
  fprintf('%d candidate rearrangements to review.\n', nx);

  % (2) write BED file for IGV
  B = [];
  B.chr = [X.chr1; X.chr2];
  B.start = [X.min1; X.min2];
  B.stop = [X.max1; X.max2];
  [B ord] = sort_struct(B,{'chr','start','stop'});
  which_bkpt = 1+(ord>nx);
  which_rearr = fidx(mod(ord-1,nx)+1);

  nb = slength(B);
  for i=1:nb, B.name{i,1} = ['Bkpt ' num2str(i)]; end
  B.review = -ones(nb,1);
  if exist(rbfname,'file')
    fprintf('Manual review was already begun for this sample.\nFile: %s\n', rbfname);
    while true
      a = upper(input('  (R)esume?  or start (N)ew review? ','s'));
      if strcmp(a,'R')
        fprintf('Resuming previous review session.\n');
        B = load_struct(rbfname,'%f%f%f%s%f',0);
        B = rename_fields(B,colx(1:5),{'chr','start','stop','name','review'});
        nb = slength(B);
        % check to ensure dRanger results file has not changed in the interim
        if (nb ~= nx*2), error('Previous review file is not a match\n'); end
        q = [Xall.chr1(which_rearr) Xall.chr2(which_rearr)];
        qq = q(:,1); qq(which_bkpt==2) = q(which_bkpt==2,2);
        if any(B.chr~=qq), error('Previous review file is not a match\n'); end
        q = [Xall.min1(which_rearr) Xall.min2(which_rearr)];
        qq = q(:,1); qq(which_bkpt==2) = q(which_bkpt==2,2);
        if any(B.start~=qq), error('Previous review file is not a match\n'); end
        q = [Xall.max1(which_rearr) Xall.max2(which_rearr)];
        qq = q(:,1); qq(which_bkpt==2) = q(which_bkpt==2,2);
        if any(B.stop~=qq), error('Previous review file is not a match\n'); end
        % passed comparison
        pos = find(B.review==-1,1);
        if isempty(pos), pos = 1; end
        break
      elseif strcmp(a,'N')
        fprintf('Discarding previous file; starting new review session.\n');
        save_struct(B,rbfname,'no_headers');
        pos = 1;
        break
      end
    end
  else
    fprintf('Starting review session.\n');
    try
      save_struct(B,rbfname,'no_headers');
    catch me
      fprintf('Error: Was unable to write bed file %s\n',rbfname);
      error('Possibly we don''t have write permission at that location?');
    end
    pos = 1;
  end

  % (3) prompt user to open IGV
  if isempty(P.reviewing_workstation_IP)
    fprintf('\nPlease enter IP address of reviewing workstation (where IGV is running)\n');
    fprintf('You can get this by looking at http://whatismyipaddress.com/\n');
    P.reviewing_workstation_IP = input('        Address? ','s');
    fprintf('Thanks.  You can specify this IP address as P.reviewing_workstation_IP\n');
  end
  while true
    try
      fprintf('\nPlease open an IGV session and load the following files:\n');
      fprintf('  %s\n  %s\n  %s\n',tbam,nbam,rbfname);
      if P.use_IGV_control
        input('Press Enter to begin.','s');
        socket = Socket(P.reviewing_workstation_IP,P.IGV_listening_port);
        IGVout = PrintWriter(socket.getOutputStream(),true);
      end
      break
    catch me
      fprintf('\nFailed to open socket to IGV.  Try closing and restarting IGV.\n');
      fprintf('         Type Ctrl-C to abort.\n');
    end
  end

  % (4) Interactive traversal of list
  modified = false;
  done = false;
  while ~done

    xidx = which_rearr(pos);
    chr = B.chr(pos); mn = B.start(pos); mx = B.stop(pos);
    if which_bkpt(pos)==1
      if Xall.str1(xidx)==0, str='+'; else str='-'; end
      ochr = Xall.chr2(xidx); omn = Xall.min2(xidx); omx = Xall.max2(xidx);
      if Xall.str2(xidx)==0, ostr='+'; else ostr='-'; end
    else
      if Xall.str2(xidx)==0, str='+'; else str='-'; end
      ochr = Xall.chr1(xidx); omn = Xall.min1(xidx); omx = Xall.max1(xidx);
      if Xall.str1(xidx)==0, ostr='+'; else ostr='-'; end
    end
    fprintf('\n%s  /  %d       chr%d:%d-%d(%s)      Current judgement: %d\n',...
      B.name{pos},nb,chr,mn,mx,str,B.review(pos));
    fprintf('  tumreads = %d  normreads = %d  pairmate = chr%d:%d-%d(%s)\n',...
      Xall.tumreads(xidx),Xall.normreads(xidx),ochr,omn,omx,ostr); 

    if P.use_IGV_control
      if chr<23, chrstring = num2str(chr);
      elseif chr==23, chrstring = 'X';
      elseif chr==24, chrstring = 'Y';
      else chrstring = '?'; end
      igvmid = (mn+mx)/2; igvleft = round(igvmid-(P.IGV_window_width/2));
      igvright = igvleft + P.IGV_window_width;
      IGVout.println(['goto chr' chrstring ':' num2str(igvleft) '-' num2str(igvright)]);
    end

    fprintf('  (F)wd, (B)kwd, (G)oto, (I)nfo, (K)eyboard, (S)ave, (D)one, or (-1,0,1,2,3,...)\n');

    while ~done
      a = upper(input('  New judgement? ','s'));

      if strcmp(a,'F')
        pos=pos+1; if pos>nb, pos = 1; end
        break

      elseif strcmp(a,'B')
        pos=pos-1; if pos==0, pos = nb; end
        break

      elseif strcmp(a,'G')
        while true
          newpos = input('  Goto breakpoint: ','s');
          try x = str2double(newpos);
          catch me, x = nan; end
          if x<1 | x>nb | round(x)~=x
            fprintf('  Please enter an integer from 1-%d\n',nb);
          else
            pos = x; break
          end
        end
        break

      elseif strcmp(a,'I')
        look(Xall,xidx);
        break

      elseif strcmp(a,'S')
        save_the_results;
        break

      elseif strcmp(a,'D')
        if modified
         while true
            a = upper(input('\nReview done.  (S)ave, (D)iscard, or (C)ontinue? ','s'));
            if strcmp(a,'S'), done = true; save_flag = true; break
            elseif strcmp(a,'D')
              while true
                aa = upper(input('  Really discard new judgements (Y/N)? ','s'));
                if strcmp(aa,'Y'), done = true; save_flag = false; break
                elseif strcmp(aa,'N'), done = false; break;
                end
              end
              break
            elseif strcmp(a,'C'), done = false; break
            end
          end
        else % not modified
          done = true; save_flag = false;
        end
        break

      elseif strcmp(a,'K')
        keyboard, break

      else % entered new judgement?
        try j = str2double(a);
        catch me, j = nan; end
        if j<-1 | round(j)~=j
          fprintf('  Please enter -1,0,1,2,3...  ');
        else
          if B.review(pos) ~= j
            B.review(pos) = j;
            modified = true;
          end
          pos = pos-1+find(B.review(pos:end)==-1,1);
          if isempty(pos), pos = 1; end
          break
        end
      end
    end
  end
        
  % done!

  if save_flag
    % (5) save results
    save_the_results
  end

end % next sample

fprintf('\ndRanger_manual_review finished.\n');

catch me; excuse(me); end

  function save_the_results
    fprintf('Saving results...\n');
    save_struct(B,rbfname,'no_headers');
    fprintf('Saved %s\n',rbfname);
    Xsave = rename_field(Xall,'filter','filterAuto');
%fprintf('Please supervise this step: possible bug\n');
%keyboard
    Xsave.filterManual1 = -ones(nxall,1);
    Xsave.filterManual2 = -ones(nxall,1);
    Xsave.filterManual1(which_rearr(which_bkpt==1)) = B.review(which_bkpt==1);
    Xsave.filterManual2(which_rearr(which_bkpt==2)) = B.review(which_bkpt==2);
    Xsave.filterManual = 1 - (Xsave.filterManual1<1 & Xsave.filterManual2<1);
    Xsave.filter = 1 - (Xsave.filterAuto==0 &  Xsave.filterManual==0);
    Xsave = sort_struct(Xsave,{'filter','normreads','tumreads'},[1 1 -1]);
    save_struct(Xsave,rrfname);
    fprintf('Saved %s\n',rrfname);
    modified = false;
  end

end

