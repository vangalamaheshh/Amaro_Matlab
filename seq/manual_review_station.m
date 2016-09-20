function manual_review_station(P)
% manual_review_station(P)
%
% Facilitates manual review of a set of loci.
%
% Required parameters:

%   P.bedfile = filename of BED file listing loci to review
%               chr start end [name] [score]
%               "score" column records judgements of manual review
%               ** no header line **
%
% Optional parameters:
%
%   P.bamfile = can be either a bam filename or a cell array of bam filenames
%               --> these files are not actually touched;
%                   they are simply opened in IGV (or the user is told to open them in IGV)
%
%   P.infofile = file with same number of records as P.bedfile
%                --> the info can be displayed by invoking the (I)nfo command
%
%   P.infofile_has_header = true/false: default is false
%
%   P.use_IGV_control = true/false:  if true, uses control of IGV through port I/O
%
%   P.reviewing_workstation_IP = no default; if not specified, the user is prompted for it.
%
%   P.IGV_listening_port = 60151 by default
%
%   P.IGV_window_width = 2000 by default; for larger loci, increases the window size.

% Mike Lawrence 2009-07-23


if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Input parameter should be a structure.'); end

P=impose_default_value(P,'bedfile',[]);
P=impose_default_value(P,'pcbedfile',[]);
P=impose_default_value(P,'bamfile',[]);
P=impose_default_value(P,'pcbamfile',[]);
P=impose_default_value(P,'snapshotdir',[]);
P=impose_default_value(P,'infofile',[]);
P=impose_default_value(P,'infofile_has_header',true);
P=impose_default_value(P,'use_IGV_control',true);
P=impose_default_value(P,'reviewing_workstation_IP',[]);
P=impose_default_value(P,'IGV_listening_port',60151);
P=impose_default_value(P,'IGV_window_width',2000);

import java.util.*;
import java.io.*;
import java.net.Socket;

fprintf('\nmanual_review_station\n\n');

if isempty(P.bedfile), error(['This function should receive a single input parameter,\n'...
  'a structure, having at least the field "bedfile", specifying the filename of a list\n'...
  'of loci to review.  The bedfile should have no header line, and each line should\n'...
  'have the following tab-delimited format:  <chr> <start> <end> [<name> [<score>]]\n\n'...
  'Fields conform to the UCSC standard: http://genome.ucsc.edu/FAQ/FAQformat#format1\n'...
  '  <chr> is the name of the chromosome (e.g. chr3, chrY)\n'...
  '  <start> is the first base of the locus, using zero-based numbering\n'...
  '  <end> is the first base after the last base of the locus\n'...
  '  <name> optional string to display\n'...
  '  <score> numeric value representing review judgement (-1 = no judgement)\n']);
end

% LOAD BEDFILE

if ~exist(P.bedfile,'file'), error('Input bedfile %s not found.',P.bedfile); end

fprintf('Loading bedfile %s\n', P.bedfile);
B = load_struct_specify_numeric_cols(P.bedfile,[2 3],0);
if isfield(B,'col1'), B = rename_field(B,'col1','chr');
else error('First column of bedfile should be chr name'); end
if isfield(B,'col2'), B = rename_field(B,'col2','start');
else error('Second column of bedfile should be start position'); end
if isfield(B,'col3'), B = rename_field(B,'col3','end');
else error('Third column of bedfile should be end position'); end
if isfield(B,'col4'), B = rename_field(B,'col4','name'); end
if isfield(B,'col5'), B = rename_field(B,'col5','score'); end

nb = slength(B);

if ~isfield(B,'name')
  B.name = repmat({''},nb,1);
end

if isfield(B,'score')
  try B = make_numeric(B,'score');
  catch me
    fprintf('Score column should be numeric!  Resetting contents to -1\n');
    B.score = -ones(nb,1);
  end
  if any(B.score ~= -1)
    fprintf('BED file contains judgements from previous session.');
    while true
    a = upper(input('  (R)esume?  or start (N)ew review? ','s'));
      if strcmp(a,'R')
        fprintf('Resuming previous review session.\n');
        pos = find(B.score==-1,1);
        if isempty(pos), pos = 1; end
        break
      elseif strcmp(a,'N')
        fprintf('Discarding previous judgements; starting new review session.\n');
        B.score = -ones(nb,1);
        break
      end
    end
  end
else
  fprintf('Starting review session.\n');
  B.score = -ones(nb,1);
  pos = 1;
end

% LOAD INFOFILE

if ~isempty(P.infofile)
  fprintf('Loading infofile %s\n',P.infofile);
  I = load_struct_specify_numeric_cols(P.infofile,[],1*(P.infofile_has_header));
  if slength(I) ~= nb
    error('Infofile should have same number of lines as bedfile.\nSet P.infofile_has_header = (true/false)');
  end
  infofile_exists = true;
else
  infofile_exists = false;
end  

% PROMPT USER TO OPEN IGV

if P.use_IGV_control
  if isempty(P.reviewing_workstation_IP)
    P.reviewing_workstation_IP = ...
      input('\nPlease enter IP address of reviewing workstation (where IGV is running): ','s');
    fprintf('Thanks.  You can specify this IP address as P.reviewing_workstation_IP\n');
    
  end
end

fprintf('\nPlease open an IGV session\n');

if P.use_IGV_control
  input('\nPress Enter to begin.','s');
  try
    socket = Socket(P.reviewing_workstation_IP,P.IGV_listening_port);
    IGVout = PrintWriter(socket.getOutputStream(),true);
    IGVin = BufferedReader(InputStreamReader(socket.getInputStream()));
  catch me
    fprintf('Failed to establish communication with IGV.  Error was:\n');
    disp(me.message);
    fprintf('Cannot use automatic IGV control\n');
    P.use_IGV_control = false;
  end
end

fprintf('\nPlease load the following files:\n');
if ~isempty(P.bamfile)
  if ~iscell(P.bamfile), P.bamfile = {P.bamfile}; end
  for i=1:length(P.bamfile)
    fprintf('  %s\n',P.bamfile{i});
  end
end
fprintf('  %s\n', P.bedfile);
if P.use_IGV_control
  IGVout.println('new');
  get_IGV_response(IGVin);  
  if ~isempty(P.pcbamfile)
    if ~iscell(P.pcbamfile), P.pcbamfile = {P.pcbamfile}; end
    for i=1:length(P.pcbamfile)
      IGVout.println(['load ' P.pcbamfile{i}]);
      get_IGV_response(IGVin);
    end
  end
  IGVout.println(['load ' P.pcbedfile]);
  get_IGV_response(IGVin);
end

if ~isempty(P.snapshotdir)
  IGVout.println(['snapshotDirectory ' P.snapshotdir]);
  get_IGV_response(IGVin);
end
 

% TRAVERSE LIST
modified = false;
done = false;
num_snapshots=0;
while ~done

    fprintf('\nLocus %d / %d    %s    ', pos, nb, B.name{pos});
    chr = B.chr{pos}; st = B.start(pos); en = B.end(pos);
    fprintf('%s:%d-%d      Current judgement: %d\n',chr,st,en,B.score(pos));

    if P.use_IGV_control
      if length(chr)<4 | ~strncmpi(chr,'chr',3), chrstring = ['chr' chr];
      else chrstring = chr; end
      chrstring = regexprep(chrstring,'chr23','chrX');
      chrstring = regexprep(chrstring,'chr24','chrY');
      if en-st < 0.8 * P.IGV_window_width
        mid = (st+en)/2;
        left = mid - (P.IGV_window_width/2);
        right = left + P.IGV_window_width;
      else
        left = st - (0.1 * P.IGV_window_width);
        right = st + (0.1 * P.IGV_window_width);
      end
      disp(['goto ' chrstring ':' num2str(round(left)) '-' num2str(round(right))]);
      IGVout.println(['goto ' chrstring ':' num2str(round(left)) '-' num2str(round(right))]);
      get_IGV_response(IGVin);
      if ~isempty(P.snapshotdir)
        IGVout.println('snapshot');
        get_IGV_response(IGVin);
        num_snapshots=num_snapshots+1;
      end
    end

    fprintf('  (F)wd, (B)kwd, (G)oto, (I)nfo, (S)ave, (D)one, or (-1,0,1,2,3,...)\n');
    while ~done
      if isempty(P.snapshotdir)
        a = upper(input('  New judgement? ','s'));
      else
        if num_snapshots<nb
          a='F';
        else
          a='D';
          modified=0;
        end
        disp(a);
      end
      
      if strcmp(a,'F')
        pos=pos+1; if pos>nb, pos = 1; end
        break

      elseif strcmp(a,'B')
        pos=pos-1; if pos==0, pos = nb; end
        break

      elseif strcmp(a,'G')
        while true
          newpos = input('  Goto locus: ','s');
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
        look(I,pos);
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
          if B.score(pos) ~= j
            B.score(pos) = j;
            modified = true;
          end
          oldpos = pos;
          pos = oldpos+find(B.score(oldpos+1:end)==-1,1);
          if isempty(pos), pos = oldpos+1; if pos>nb, pos = 1; end, end
          break
        end
      end
    end
end
        
% done!

IGVout.close;
IGVin.close;
socket.close;

if save_flag
  % (5) save results
  save_the_results
end

fprintf('\nManual review finished.\n');

  function save_the_results
    fprintf('Saving results...\n');
    save_struct(B,P.bedfile,'no_headers');
    fprintf('Saved %s\n',P.bedfile);
    modified = false;
  end

end

%%-----------------------------------
function st=get_IGV_response(IGVin)

st=[];
tests=30;

i=0;
while ~IGVin.ready()
  pause(1);
  i=i+1;
  fprintf(1,'.');
  if i>tests
      error('waited too long for IGV');
  end
end
st=IGVin.readLine();
disp(['IGV: ' char(st) ' (' num2str(i) ')']);

end
  
