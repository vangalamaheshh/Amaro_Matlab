function [S P F] = riker(P,F,params)
% [S P F] = riker(P,F,params)
%
% Mike Lawrence 2009-09-15

fprintf('%s\n',datestr(now));

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'freeze_file','freeze.txt');
params = impose_default_value(params,'ignore_running_file',false);
if nargin==1 && isnumeric(P)
  if isfield(params,'configno') && params.configno ~= P, error('configno conflict'); end
  params.configno = P;
  clear P;
else
  params = impose_default_value(params,'configno',1);
end

if params.configno==1
  workspace = 'trunk';
  build = 'hg18';
  fdir = '/xchip/cga1/firehose_output/trunk/Individual';
  outdir = '/xchip/cga1/lawrence/riker';
  banner = 'RIKER: trunk';
elseif params.configno==2
  workspace = 'trunk_hg19';
  build = 'hg19';
  fdir = '/xchip/cga1/firehose_output/trunk_hg19/Individual';
  outdir = '/xchip/cga1/lawrence/riker_hg19';
  banner = 'RIKER: trunk_hg19';
elseif params.configno==3
  workspace = 'mouse_mm9';
  build = 'mouse';
  fdir = '/xchip/cga1/firehose_output/mouse_mm9/Individual';
  outdir = '/xchip/cga1/lawrence/riker_mouse';
  banner = 'RIKER: mouse_mm9';
else
  error('Unknown configuration number %d',params.configno);
end

ensure_dir_exists(outdir);

if ~params.ignore_running_file
  runningfile = [outdir '/running.txt'];
  if exist(runningfile,'file')
    d = dir(runningfile);
    fprintf('Riker is already running: since %s\n',d.date);
    v = datevec(now-d.datenum);
    if v(1)>0 || v(2)>0 || v(3)>0 || v(4)>12
      fprintf('Has been running more than 12 hours: will try running again.\n');
    else
      fprintf('Will not start new run.\n');
      if ~exist('P','var'), P=[]; end
      if ~exist('F','var'), F=[]; end
      S=[];
      return;
    end
  end
  fprintf('Riker starting...\n');
  save_textfile('running',runningfile);
end

tottime = tic;

if ~exist(params.freeze_file)
  fprintf('Freeze file not found.\n');
  params.freeze_file = '';
end

rikerdir = outdir;
if ~exist('P','var') || ~exist('F','var')
  fprintf('crawl_firehose_2:  ');
  tic
%  F = crawl_firehose(fdir);
  F = crawl_firehose_2(workspace,outdir);
  toc
  fprintf('crawl_picard:  ');
  tic;
  P = crawl_picard(rikerdir);
  toc;
end

fprintf('Picard and Firehose info obtained.\n');

P_all = P;
if isfield(P_all,'ignore')
  P = reorder_struct(P_all,~P_all.ignore);
end

P.fidx = listmap(P.bampath,F.bampath);
F.pidx = listmap(F.bampath,P.bampath);

S = [];
[tmp si sj] = unique(stringsplice([P.sample P.tech],1,':'));
S.sample = P.sample(si);
S.tech = P.tech(si);
ns = slength(S);

if ~isempty(params.freeze_file)
  fname = [outdir '/' params.freeze_file];
  d = dir(fname);
  if isempty(d), error('Not found: %s\n',fname); end
  freezedate = d.date;
  Fz = load_struct(fname);
  Fz = make_numeric(Fz,{'tech'});
  freezeidx = listmap(S.sample,Fz.sample);
else
  Fz = [];
end

fprintf('Correlating Picard and Firehose...');

for i=1:ns
  % what is in Picard?
  pidx = find(sj==i);  
  S.project{i,1} = concat(unique(P.project(pidx)),',');
  S.initiative{i,1} = concat(unique(P.initiative(pidx)),',');
  S.guess_individual{i,1} = concat(unique(P.guess_individual(pidx)),',');
  S.guess_tn{i,1} = concat(unique(P.guess_tn(pidx)),',');
  [tmp tmpidx] = max(P.bamtimenum(pidx));
  pidx_l = pidx(tmpidx);
  if isempty(pidx_l)
    S.latest_version{i,1} = 'none';
  else
    S.latest_version{i,1} = P.version{pidx_l};
  end
  pidx_f = pidx(find(strcmp(P.finished(pidx),'true')));
  [tmp tmpidx] = max(P.bamtimenum(pidx_f));
  pidx_lf = pidx_f(tmpidx);
  if isempty(pidx_lf)
    S.latest_finished_version{i,1} = 'none';
    S.latest_finished_build{i,1} = 'none';
    S.latest_finished_bamdate{i,1} = 'none';
    S.latest_finished_bamsize(i,1) = -1;
    S.latest_finished_numlanes(i,1) = -1;
  else
    S.latest_finished_version{i,1} = P.version{pidx_lf};
    S.latest_finished_build{i,1} = P.build{pidx_lf};
    S.latest_finished_bamdate{i,1} = P.bamdate{pidx_lf};
    S.latest_finished_bamsize(i,1) = P.bamsize(pidx_lf);
    S.latest_finished_numlanes(i,1) = P.numlanes(pidx_lf);
  end
  % what is in the freeze?
  if ~isempty(Fz) & ~isnan(freezeidx(i))
    fridx = freezeidx(i);
  else
    fridx = [];
  end
  % what is in Firehose?
  fidx = find(ismember(F.pidx,pidx));
  if isempty(fidx);
    S.firehose_individual{i,1} = '?';
    S.firehose_version{i,1} = 'none';
    S.firehose_build{i,1} = 'none';
    S.firehose_tn{i,1} = '?';
    S.firehose_bamdate{i,1} = 'none';
    S.firehose_bamsize(i,1) = -1;
    S.firehose_numlanes(i,1) = -1;
    if isempty(pidx_lf)
      S.firehose_status{i,1} = 'PENDING';
    else
      if ~isempty(Fz)
        S.firehose_status{i,1} = 'FROZEN';
      else
        S.firehose_status{i,1} = 'MISSING';
      end
    end
  else
    tmp=[];
    tmp.in = unique(F.indiv(fidx));
    tmp.ve = unique(P.version(F.pidx(fidx)));
    tmp.bl = unique(P.build(F.pidx(fidx)));
    tmp.tn = unique(F.tn(fidx));
    tmp.bd = unique(P.bamdate(F.pidx(fidx)));
    S.firehose_individual{i,1} = concat(tmp.in,',');
    S.firehose_version{i,1} = concat(tmp.ve,',');
    S.firehose_build{i,1} = concat(tmp.bl,',');
    S.firehose_tn{i,1} = concat(tmp.tn,',');
    S.firehose_bamdate{i,1} = concat(tmp.bd,',');
    S.firehose_bamsize(i,1) = max(P.bamsize(F.pidx(fidx)));
    S.firehose_numlanes(i,1) = max(P.numlanes(F.pidx(fidx)));
    if length(tmp.in)>1 || length(tmp.ve)>1 || length(tmp.bl)>1 || length(tmp.tn)>1 || length(tmp.bd)>1
      S.firehose_status{i,1} = 'REDUNDANT';
    else
      % FH points to only one version: see if it's outdated or not
      fv = sscanf(S.firehose_version{i},'v%d');
      lv = sscanf(S.latest_finished_version{i},'v%d');
      if ~isempty(fridx)
        if strcmp(Fz.firehose_version{fridx},'none')
          frzv = 0;
        else
          frzv = sscanf(Fz.firehose_version{fridx},'v%d');
        end
      else
        frzv = [];
      end
      if isempty(frzv)
        if isempty(lv) || fv>=lv
          S.firehose_status{i,1} = 'OK';
      else
        S.firehose_status{i,1} = 'OUTDATED';
        end
      else
        if fv == frzv
          if isempty(lv) || fv >= lv
            S.firehose_status{i,1} = 'OK';
          else
            S.firehose_status{i,1} = 'FROZEN';
          end
        else % fv ~= frzv
          S.firehose_status{i,1} = 'DEFROSTED';
        end
      end
    end
  end
end

fprintf('done.\n');

% for individuals not in Firehose, use the "guessed" individual name + TN

S.individual = S.firehose_individual;
S.tn = S.firehose_tn;
S.infirehose = true(slength(S),1);
S.guessed_individual = false(slength(S),1);
idx = find((strcmp('?',S.individual)|strcmp('?',S.tn)) & ~(strcmp('?',S.guess_individual)|strcmp('?',S.guess_tn)));
S.individual(idx) = S.guess_individual(idx);
S.guessed_individual(idx) = true;
S.tn(idx) = S.guess_tn(idx);
S.infirehose(idx) = false;

fprintf('Writing webpage to %s...',outdir);

% write webpage

status = {...    % this list also defines the preferred display order
    'PENDING',...
    'MISSING',...
    'OUTDATED',...
    'REDUNDANT',...
    'DEFROSTED',...
    'FROZEN',...
    'OK',...
};

formatted_status = {...
    '<font color="grey">PENDING...</font>',...
    '<b><font color="red">MISSING</font></b>',...
    '<font color="red">OUTDATED</font>',...
    '<font color="blue">REDUNDANT LINKS</font>',...
    '<b><font color=00bbff>DEFROSTED</b></font> <img src="defrost.jpg">',...
    'FROZEN',...
    'OK',...
};

status_note = {...
    'No BAM is ready in Picard yet (based on presence of "finished.txt").',...
    '',...
    '',...
    'Firehose has multiple links to this BAM.',...
    'The freeze specifies a different version (or no version at all) of this BAM.',...
    'Firehose is OK, but newer version of this BAM is available.',...
    '',...
};

S.formatted_firehose_status = map_across(S.firehose_status,status,formatted_status);

indexname = [outdir '/index.html'];
f = fopen(indexname,'wt');
fprintf(f,'<html><head><title>%s</title></head>\n',banner);
fprintf(f,'<table><td width=120><img width=100 src="riker.jpg"></td>');
fprintf(f,'<td>\n');
fprintf(f,'  <h1>%s</h1>\n',banner);
fprintf(f,'  <h2>Picard samples in Firehose</h2>\n');
if ~isempty(Fz), fprintf(f,'  Freeze in effect: %s\n',freezedate); end
fprintf(f,'</td></table>\n');
fprintf(f,'<hr>\n');

techname = {'WGS','Hybrid selection','Other'};
hyb = grepi('Hybrid|Exome',S.initiative,1);
wgs = grep('WGS',S.initiative,1);
S.tech = 3*ones(slength(S),1);
S.tech(wgs) = 1;
S.tech(hyb) = 2;
S.tech(intersect(hyb,wgs)) = 3;
% load list of initiatives that we are forcing to be listed under "Other"
othI = load_struct([rikerdir '/initiatives_forced_to_other.txt']);
S.tech(ismember(S.initiative,othI.initiative)) = 3;

% create pages
pageidx = 0;
for techidx=1:3
  St = reorder_struct(S,S.tech==techidx);
  if slength(St)==0, continue; end
  fprintf(f,'<p><table width=1050><tr><td width=340><font size=+3><b>%s</b></font></td>',techname{techidx});
  z = unique(St.firehose_status);
  zidx = listmap(z,status);  % sort by prefered display order
  [tmp zord] = sort(zidx);
  z = z(zord);
  for b=1:length(z)+1
    if b<=length(z), name = z{b}; else name = 'TOTAL'; end
    fprintf(f,'<td width=100 align=center><b>%s</b></td>',name);
    totalsamps(b) = 0;
    totalpairs(b) = 0;
  end
  fprintf(f, '</tr><tr></tr>\n');
  [init ii ij] = unique(St.initiative);
  for i=1:length(init)
    pageidx = pageidx + 1;
    fprintf(f,'<tr><td><a href="init%d.html">%s</a></td>',pageidx,init{i});
    Sti = reorder_struct(St,ij==i);
    % find matched T/N pairs in this tech:intiative
    Sti.ispaired = false(slength(Sti),1);
    [u ui uj] = unique(Sti.individual);
    for j=1:length(u)
      idx = find(uj==j);
      tidx = idx(strcmp('tumor',Sti.tn(idx)));
      nidx = idx(strcmp('normal',Sti.tn(idx)));
      while ~isempty(tidx) && ~isempty(nidx)
        dwt = parse(Sti.sample(tidx),'TCGA-..-....-...-..(.).*',{'DW'}); dwt = dwt.DW;
        dwn = parse(Sti.sample(tidx),'TCGA-..-....-...-..(.).*',{'DW'}); dwn = dwn.DW;
        u = unique([dwt;dwn]);
        matcht = 1; matchn = 1;
        for k=1:length(u)
          idxt = find(strcmp(u{k},dwt),1);
          idxn = find(strcmp(u{k},dwn),1);
          if ~isempty(idxt) && ~isempty(idxn)
            matcht = idxt(1); matchn = idxn(1);
            break;
        end,end
        Sti.ispaired([tidx(matcht);nidx(matchn)]) = true;
        tidx(matcht)=[]; nidx(matchn)=[];
      end        
    end
    for b=1:length(z)+1
      if b<=length(z)
        bidx = find(strcmp(z{b},Sti.firehose_status));
      else
        bidx = 1:slength(Sti);
      end
      nsamps = length(bidx);
      fprintf(f,'<td align=center>%d',nsamps);
      totalsamps(b) = totalsamps(b) + nsamps;
      npairsz = sum(Sti.ispaired(bidx))/2;
      totalpairs(b) = totalpairs(b) + npairsz;
      if nsamps>0 && npairsz>0
        fprintf(f,' (');
        if floor(npairsz)>0, fprintf(f,'%d',floor(npairsz)); end
        if npairsz>floor(npairsz), fprintf(f,'&frac12;'); end
        fprintf(f,')');
      end
      fprintf(f,'</td>');
    end
    fprintf(f,'</tr>\n');
    initname = sprintf([outdir '/init%d.html'],pageidx);
    g = fopen(initname,'wt');
    fprintf(g,'<html><head><title>%s</title></head>\n',banner);
    fprintf(g,'<h1>%s</h1>\n',init{i});
    fprintf(g,'<h2>Firehose status by sample</h2>\n');
    if ~isempty(Fz) && any(Fz.tech==techidx)
      fprintf(g,'Freeze in effect: %s\n',freezedate);
    end
    fprintf(g,'<hr><table width=1200><tr>');
    fprintf(g,'<td width=65><b>project<b></td>');
    fprintf(g,'<td width=150><b>sample</b></td>');
    fprintf(g,'<td width=65><b>Picard</b></td>');
    fprintf(g,'<td width=50><b>build</b></td>');
    fprintf(g,'<td width=30><b>GB</b></td>');  
    fprintf(g,'<td width=40><b>lanes</b></td>');
    fprintf(g,'<td width=100><b>date</b></td>');
    fprintf(g,'<td width=150><b>individual</b></td>');
    fprintf(g,'<td width=30><b>T/N</b></td>');
    fprintf(g,'<td width=30><b>pair</b></td>');
    fprintf(g,'<td width=65><b>Firehose</b></td>');
    fprintf(g,'<td width=30><b>GB</b></td>');
    fprintf(g,'<td width=40><b>lanes</b></td>');
    fprintf(g,'<td width=100><b>date</b></td>');
    fprintf(g,'<td width=150><b>status</b></td>');
    fprintf(g,'</tr><tr></tr>\n');
%    [tmp sortidx] = sort_struct(Sti,{'infirehose','ispaired','individual','tn','project','sample'},[1 1 1 -1 1 1]);
    [tmp sortidx] = sort_struct(Sti,{'ispaired','infirehose','individual','tn','project','sample'},[1 1 1 -1 1 1]);

    for j=1:length(sortidx),k = sortidx(j);
      fprintf(g,'<tr>');
      % in Picard
      fprintf(g,'<td>%s</td>',Sti.project{k});
      fprintf(g,'<td>%s</td>',Sti.sample{k});
      fprintf(g,'<td>%s',Sti.latest_finished_version{k});
      if ~strcmp(Sti.latest_version{k},Sti.latest_finished_version{k})
        fprintf(g,' (%s)',Sti.latest_version{k});
      end
      fprintf(g,'</td>');
      fprintf(g,'<td>');
      wrongbuild = ~strcmp(build,Sti.latest_finished_build{k});
      if wrongbuild, fprintf(g,'<font color=red><b>'); end
      fprintf(g,'%s',Sti.latest_finished_build{k});
      if wrongbuild, fprintf(g,'</font></b>'); end
      fprintf(g,'</td>');
      if Sti.latest_finished_bamsize(k)>-1
        fprintf(g,'<td>%d</td>',round(Sti.latest_finished_bamsize(k)/1073741824));
        fprintf(g,'<td>%d</td>',round(Sti.latest_finished_numlanes(k)));
        fprintf(g,'<td>%s</td>',Sti.latest_finished_bamdate{k});
      else
        fprintf(g,'<td></td><td></td><td></td>');
      end
      % in Firehose
      fprintf(g,'<td>');
      if Sti.guessed_individual(k), fprintf(g,'<font color=grey>'); end
      fprintf(g,'%s',Sti.individual{k});
      if Sti.guessed_individual(k), fprintf(g,'</font>'); end
      fprintf(g,'</td>');
      fprintf(g,'<td>');
      if Sti.guessed_individual(k), fprintf(g,'<font color=grey>'); end
      fprintf(g,'%s',upper(Sti.tn{k}(1)));
      if Sti.guessed_individual(k), fprintf(g,'</font>'); end
      fprintf(g,'</td>');
      if strcmp(Sti.tn{k},'?')
        fprintf(g,'<td></td>');
      elseif Sti.ispaired(k)
        fprintf(g,'<td>+</td>');
      else
        fprintf(g,'<td><font color=red><b>&times;</b></font></td>');
      end
      fprintf(g,'<td>%s</td>',Sti.firehose_version{k});
      if Sti.firehose_bamsize(k)>-1
        fprintf(g,'<td>%d</td>',round(Sti.firehose_bamsize(k)/1073741824));
        fprintf(g,'<td>%d</td>',round(Sti.firehose_numlanes(k)));
        fprintf(g,'<td>%s</td>',Sti.firehose_bamdate{k});
      else
        fprintf(g,'<td></td><td></td><td></td>');
      end
      % Status
      fprintf(g,'<td>%s</td>',Sti.formatted_firehose_status{k});
      fprintf(g,'</tr>\n');
    end
    fprintf(g,'</table><hr>');
    flag = false;
    for i=1:length(status)
      if any(strcmp(Sti.firehose_status,status{i})) && ~isempty(status_note{i})
        if flag, fprintf(g,'<hr>'); end
        fprintf(g,'%s &mdash; %s\n',status{i},status_note{i});
        flag = true;
      end
    end
    if flag, fprintf(g,'<hr>\n'); end
    fprintf(g,'<p>Last updated: %s</html>\n',datestr(now));
    fclose(g);
  end

  % print TOTALS for this tech type
  if (techidx<3)
    fprintf(f,'<tr><td></td>');
    for b=1:length(z)+1, fprintf(f,'<td><hr></td>'); end
    fprintf(f,'</tr>\n');
    fprintf(f,'<tr><td></td>');
    for b=1:length(z)+1
      fprintf(f,'<td align=center><b>%d',totalsamps(b));
      if totalsamps(b)>0 && totalpairs(b)>0
        fprintf(f,' (');
        if floor(totalpairs(b))>0, fprintf(f,'%d',floor(totalpairs(b))); end
        if totalpairs(b)>floor(totalpairs(b)), fprintf(f,'&frac12;'); end
        fprintf(f,')');
      end
      fprintf(f,'</b></td>');
    end
    fprintf(f,'</tr>\n');
  end

  fprintf(f,'</table>\n');
  if (techidx<3), fprintf(f,'<hr>\n'); end

end % next techidx (WGS, Hybrid Selection, Other)

fprintf(f,'<p><table><tr><td width=340></td><td>samples (pairs)</td></tr></table>\n');

fprintf(f,'<hr>\n');
fprintf(f,'Ignore lists:');
fprintf(f,'<font color=white>&mdash;</font><a href="non-cancer_initiatives_to_ignore.txt">non-cancer initiatives</a>');
fprintf(f,'<font color=white>&mdash;</font><a href="not-this-build_initiatives_to_ignore.txt">not-this-build initiatives</a>');
fprintf(f,'<font color=white>&mdash;</font><a href="samples_to_ignore.txt">samples</a>\n');

fprintf(f,'<hr><p>Last updated: %s</html>\n',datestr(now));
fclose(f);

fprintf('Done:\nFinished at %s:   ',datestr(now)); toc(tottime);


save([outdir '/riker.mat']);

if ~params.ignore_running_file
  delete(runningfile);
end
