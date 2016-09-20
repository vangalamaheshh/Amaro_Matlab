function f = get_TN_scale_factor(sample)
[t n] = get_bam_stats(sample);
f = t.mapped / n.mapped;

