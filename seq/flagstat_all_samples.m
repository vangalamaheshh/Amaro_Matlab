function flagstat_all_samples

fprintf('flagstat_all_samples\n');
fprintf('  --> compiling list of all samples\n');
samples = list_all_samples;
flagstat(samples);
