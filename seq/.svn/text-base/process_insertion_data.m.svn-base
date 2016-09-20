function [t n result] = process_insertion_data(sample,chr,jobleft,jobright,flag)

% core statistic calculation
sw = 300;   % locality window
sw2 = 50;   % smoothing window

% peak finding
cutoff = 6;           % minimum peak height
peakmask = 600;       % minimum peak separation
chunksize = 1000;    % scan width
overlap = 200;        % scan overlap

fprintf('Sample %s\n',sample);
if ~contains(sample,'/')    % convert from Firehose individual to lawrence directory
  samps = load_struct('/xchip/tcga_scratch/lawrence/dRanger/20100301_local/all_wgs_samps.txt');
  sample = mapacross({sample},samps.indiv,samps.dir);
  sample = sample{1};
end

t = [];
n = [];

t.dir = ['/xchip/cga1/lawrence/' sample '/dRl.T'];
n.dir = ['/xchip/cga1/lawrence/' sample '/dRl.N'];

  cname = num2str(chr);

  fprintf('Chromosome %s\n',cname);

  t.ffile = [t.dir '/chr' cname '.hmccF'];
  t.rfile = [t.dir '/chr' cname '.hmccR'];
  n.ffile = [n.dir '/chr' cname '.hmccF'];
  n.rfile = [n.dir '/chr' cname '.hmccR'];
  dtf = dir(t.ffile); dtr = dir(t.rfile);
  dnf = dir(n.ffile); dnr = dir(n.rfile);
  if isempty(dtf) || isempty(dnf) || isempty(dtr) || isempty(dnr)
    error('  Chr %s not complete\n',cname);
  end
  t.flen = dtf.bytes/8;
  t.rlen = dtr.bytes/8;
  n.flen = dnf.bytes/8;
  n.rlen = dnr.bytes/8;
  minlen = min([t.flen t.rlen n.flen n.rlen]);

    if jobright>minlen, jobright = minlen; doneflag = true; end
    fprintf('  chr%d:%d-%d\n',chr,jobleft,jobright);
    tic,fprintf('    Loading... ');
      t.fhm = get_block(t.ffile,'long',jobleft-1,jobright-1);
      t.rhm = get_block(t.rfile,'long',jobleft-1,jobright-1);
      n.fhm = get_block(n.ffile,'long',jobleft-1,jobright-1);
      n.rhm = get_block(n.rfile,'long',jobleft-1,jobright-1);
    toc

    if flag
      tic,fprintf('    Smoothing... ');
        t.fhm_sm = smooth(diff([t.fhm;0]),sw);
        t.rhm_sm = smooth(diff([t.rhm;0]),sw);
        n.fhm_sm = smooth(diff([n.fhm;0]),sw);
        n.rhm_sm = smooth(diff([n.rhm;0]),sw);
      toc
    end

    tic,fprintf('    Processing... ');
      [t.dhmz t.shmz] = calculate_insertion_stats(t.fhm,t.rhm,sw,sw2);
      [n.dhmz n.shmz] = calculate_insertion_stats(n.fhm,n.rhm,sw,sw2);
    toc

    result = [];
    tt = t.dhmz;
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
        [nval nidx] = max(n.dhmz(left:right));
        result = [result;chr jobleft+tidx-1 tval nval t.shmz(tidx) n.shmz(tidx)]; tt(left:right) = 0;
      end
    toc

