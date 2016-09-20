function find_insertions(sample)
% Mike Lawrence 2010-03-02

% chromosome scanning
jobsize = 10e6;
joboverlap = 1000;

% core statistic calculation
sw = 300;   % locality window
sw2 = 50;   % smoothing window

% peak finding
cutoff = 6;           % minimum peak height
peakmask = 600;       % minimum peak separation
chunksize = 10000;    % scan width
overlap = 200;        % scan overlap

fprintf('Sample %s\n',sample);
if ~contains(sample,'/')    % convert from Firehose individual to lawrence directory
  samps = load_struct('/xchip/tcga_scratch/lawrence/dRanger/20100301_local/all_wgs_samps.txt');
  sample = mapacross({sample},samps.indiv,samps.dir);
  sample = sample{1};
end

td = ['/xchip/cga1/lawrence/' sample '/dRl.T'];
nd = ['/xchip/cga1/lawrence/' sample '/dRl.N'];
of = ['/xchip/cga1/lawrence/' sample '/insertions3.txt'];

result = [];
for chr=1:24

  cname = num2str(chr);
  fprintf('Chromosome %s\n',cname);
  tf_file = [td '/chr' cname '.hmccF'];
  tr_file = [td '/chr' cname '.hmccR'];
  nf_file = [nd '/chr' cname '.hmccF'];
  nr_file = [nd '/chr' cname '.hmccR'];
  dtf = dir(tf_file); dtr = dir(tr_file);
  dnf = dir(nf_file); dnr = dir(nr_file);
  if isempty(dtf) || isempty(dnf) || isempty(dtr) || isempty(dnr)
    fprintf('  Chr %s not complete\n',cname);
    continue;
  end
  tflen = dtf.bytes/8;
  trlen = dtr.bytes/8;
  nflen = dnf.bytes/8;
  nrlen = dnr.bytes/8;
  minlen = min([tflen trlen nflen nrlen]);

  jobleft = 1; jobright = jobsize + joboverlap;
  doneflag = false;
  while(~doneflag)
    if jobright>minlen, jobright = minlen; doneflag = true; end
    fprintf('  chr%d:%d-%d\n',chr,jobleft,jobright);
    tic,fprintf('    Loading... ');
      tf = get_block(tf_file,'long',jobleft-1,jobright-1);
      tr = get_block(tr_file,'long',jobleft-1,jobright-1);
      nf = get_block(nf_file,'long',jobleft-1,jobright-1);
      nr = get_block(nr_file,'long',jobleft-1,jobright-1);
    toc

    tic,fprintf('    Processing... ');
      [tt ts] = calculate_insertion_stats(tf,tr,sw,sw2);
      [nn ns] = calculate_insertion_stats(nf,nr,sw,sw2);
    toc
    
    tic,fprintf('    Finding peaks... ');
      pos1 = 1; pos2 = chunksize + overlap;
      len = length(tt);
      while(true)
        [tval tidx] = max(tt(pos1:pos2));
        if tval<cutoff
          if pos2 >= len, break; end
          pos1 = pos1 + chunksize;
          pos2 = pos2 + chunksize;
          pos2 = min(pos2,len);
          continue;
        end
        tidx = tidx + (pos1-1);
        left = max(1,tidx-peakmask); right = min(len,tidx+peakmask);
        [nval nidx] = max(nn(left:right));
        result = [result;chr jobleft+tidx-1 tval nval ts(tidx) ns(tidx)]; tt(left:right) = 0;
      end
    toc

  jobleft = jobleft + jobsize;
  jobright = jobright + jobsize;
  end % next chromosome chunk
  
end % next chromosome

tic,fprintf('\nSaving... ');
save_matrix(result,of);
toc

fprintf('\nDone\n');

