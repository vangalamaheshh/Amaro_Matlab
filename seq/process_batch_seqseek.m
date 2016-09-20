function pr = process_batch_seqseek(samps,indir)
for i=1:slength(samps), disp(samps.name{i});
  pr{i} = process_seekseq_data([indir '/' samps.name{i}]);
end
