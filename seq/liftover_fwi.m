function liftover_fwi(fwi_in,fwi_out,build1,build2)
% liftover_fwi(fwi_in,fwi_out,build1,build2)
%
% for invalid positions, outputs consecutive positions on "chromosome 10000"
%

if nargin==2
  build1=18; build2=19;
  fprintf('Assuming you want hg18->hg19\n');
elseif nargin~=4
  fprintf('Need fwi_in,fwi_out,build1,build2\n');
end

b1 = interpret_build(build1);
b2 = interpret_build(build2);
if ~((b1==18 && b2==19) || (b1==19 && b2==18))
  error('Unsupported combination of builds\n');
end

if b1==18 && b2==19
  fwb = '/cga/tcga-gsc/home/lawrence/db/hg18/liftOverToHg19/all.fwb';
elseif b1==19 && b2==18
  fwb = '/cga/tcga-gsc/home/lawrence/db/hg19/liftOverToHg18/all.fwb';
else
  error('wha?');
end

F = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwb);

if ischar(fwi_in), fwi_in={fwi_in}; end
if ischar(fwi_out), fwi_out={fwi_out}; end
if length(fwi_in)~=length(fwi_out), error('length(fwi_in)~=length(fwi_out)'); end
nfiles = length(fwi_in);

for fi=1:nfiles, fprintf('Converting %s -> %s\n',fwi_in{fi},fwi_out{fi});

  x = load_struct_noheader(fwi_in{fi},{'chr','start','end'});
  x.chr = convert_chr(x.chr);
  x = make_numeric(x,{'start','end'});
  xlen = x.end-x.start+1;
  xtotlen = sum(xlen);
  nx = slength(x);

  y = cell(nx,1);
  BAD_CHR = 10000; bad_pos = 1;
  for i=1:nx, if ~mod(i,1000), fprintf('%d/%d ',i,nx); end
    st = F.get(x.chr(i),x.start(i),x.end(i));
    len = x.end(i)-x.start(i)+1;
    chr = x.chr(i)*ones(len,1);
    bad = (st<1);
    nbad = sum(bad);
    chr(bad) = BAD_CHR;
    st(bad) = bad_pos+(0:nbad-1);
    bad_pos = bad_pos + nbad;
    dchr = diff(chr);
    dst = diff(st);
    keep = [true;(dchr~=0 | dst~=+1)];
    y{i}.chr = chr(keep);
    y{i}.st = st(keep);
    y{i}.en = st([keep(2:end);true]);
  end, fprintf('\n');

  y = concat_structs(y);
  ylen = y.en-y.st+1;
  ytotlen = sum(ylen);
  if ytotlen~=xtotlen, error('ytotlen~=xtotlen'); end
  fprintf('   total territory = %ld\n',xtotlen);

  save_struct_noheader(y,fwi_out{fi});
end

F.close();

