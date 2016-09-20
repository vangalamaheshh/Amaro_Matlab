function annotate_with_conservartion(q, q_hg18, out, out_hg18)

beep off

%addpath /xchip/cga1/petar/CancerGenomeAnalysis/trunk/matlab/
%addpath /xchip/cga1/petar/CancerGenomeAnalysis/trunk/matlab/seq
%addpath /xchip/cga2/lawrence/cga/trunk/matlab/seq
%addpath /xchip/cga2/petar/CancerGenomeAnalysis/trunk/matlab/seq
%addpath /xchip/cga1/petar/CancerGenomeAnalysis/trunk/matlab/mike
%addpath /xchip/cga1/petar/bmr
javaclasspath({...
'/xchip/cga2/lawrence/cga/trunk/analysis_pipeline/tools/dist/FixedWidthBinary.jar';...
});
import org.broadinstitute.cga.tools.seq.*;
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;


%import('org.broadinstitute.cga.tools.seq.FixedWidthBinary');


q = load_struct(q);
q.chr = convert_chr(q.Chromosome);
q.start = str2double(q.Start_position);
q.cons = nan(slength(q),1);

z = org.broadinstitute.cga.tools.seq.FixedWidthBinary('/xchip/cga1/lawrence/db/hg19/conservation46/all.fwb');

q.cons = z.get(q.chr,q.start);
q.cons = (double(q.cons)-50)/(25/3);      % convert back to negative to positive
q.cons(q.cons==18) = NaN; % missing data

save_struct(q, out);

q_hg18 = load_struct(q_hg18);
q_hg18.chr = convert_chr(q_hg18.Chromosome);
q_hg18.start = str2double(q_hg18.Start_position);
q_hg18.cons = nan(slength(q_hg18),1);

z = org.broadinstitute.cga.tools.seq.FixedWidthBinary('/xchip/cga1/lawrence/db/hg18/conservation44/all.fwb');

q_hg18.cons = nan(slength(q_hg18),1);

q_hg18.cons = z.get(q_hg18.chr,q_hg18.start);
q_hg18.cons = (double(q_hg18.cons)-50)/(25/3);      % convert back to negative to positive
q_hg18.cons(q_hg18.cons==18) = NaN; % missing data

save_struct(q_hg18, out_hg18);

quit
