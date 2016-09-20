function cb = temp_mutrate_vs_cov(fwb,covmax,outfile,M)

COV = FWB(fwb);
fprintf('Getting coverage across territory\n');
t = load_matrix('/xchip/cga1/lawrence/mut/analysis/20121031_pancan/fwb2/refseq.hg19.targs.1.fwi');
for ti=1:1:size(t,1),if ~mod(ti-1,1000), fprintf('%d ',ti); end
  cov = double(COV.get(t(ti,1),t(ti,2),t(ti,3)));
  cov(cov==-1)=0;
  tmp = histc(cov,0:covmax);
  if ti==1, ch=tmp; else ch=ch+tmp; end
end,fprintf('\n');

if ~exist('M','var')
  fprintf('Loading mutations\n');
  load('/xchip/cga1/lawrence/mut/analysis/20121031_pancan/data.v2.patients.2.2.mat','M');
end
fprintf('Getting coverage across mutations\n');
mcov = COV.get(M.mut.chr,M.mut.start);

fprintf('Computing rates\n');
cb=[];
cb.min = (0:covmax)';
cb.max = cb.min;
cb.nbp = nan(slength(cb),1);
midx = find(M.mut.splicedist<=2);
cb.n = histc(mcov(midx),cb.min);
for i=1:slength(cb)
  cb.nbp(i) = sum(ch(cb.min(i)+1:cb.max(i)+1));
end
cb.N = cb.nbp*5330;
cb.mutrate = 1e6*cb.n./cb.N;
cb.mid = cb.min;
cb.cov = cb.mid/covmax;
cb.smutrate = smooth(cb.mutrate,length(cb.mutrate));

if exist('outfile','var'), save(outfile,'ch','mcov','cb'); end
