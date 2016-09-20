function [X R mut] = discordup_filter(mut,bam,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'qualcutoff',20);
P = impose_default_value(P,'mapqcutoff',20);
P = impose_default_value(P,'build','hg19');
P = impose_default_value(P,'refdir',['/xchip/cga1/annotation/db/ucsc/' P.build]);
P = impose_default_value(P,'blacklist','none');
P = impose_default_value(P,'manual_inspection',false);
P = impose_default_value(P,'quiet',true);

P = impose_default_value(P,'duplicate_criterion','start1+start2');
dcnames = {'start1+start2','start1+start2+insertsize','start1+insertsize','start2+insertsize'};
dupcrit = decell(mapacross({P.duplicate_criterion},dcnames,{[1 2],[1 2 3],[1 3],[2 3]},nan));
if isnan(dupcrit), disp(dcnames), error('P.duplicate_criterion must be one of the above'); end

if ~ischar(bam), error('"bam" should be a path'); end
if ~isstruct(mut), error('"mut" should be a struct'); end
demand_fields(mut,{'chr','start'});
if ~isfield(mut,'i_tumor_f'), mut.i_tumor_f = nan(slength(mut),1); end

import org.broadinstitute.cga.tools.seq.*;
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;

BAM = org.broadinstitute.cga.tools.seq.BamGrasp();
if P.quiet, BAM.setQuietModeOn(); end
BAM.openFile(String(bam),String(P.blacklist),String(P.refdir));
P.include_unmapped_reads = false;
P.include_duplicate_reads = true;
P.include_aux_cols = true;

nm = slength(mut);
X = zeros(10,10,nm);

for i=1:nm, if ~mod(i,10), fprintf('%d/%d ',i,nm); end
  chr=mut.chr(i);pos=mut.start(i); [R B S] = BamGrasp_load_region(BAM,chr,pos,P);
  if P.manual_inspection, Borig=B; end
  B = B(B(:,4)==pos & B(:,2)>=P.qualcutoff,1:3); B=[R(B(:,3),[4 11 13 14 8]) B]; B=B(B(:,5)>=P.mapqcutoff,:);
  B(B<0)=nan; [u ui uj] = unique(B(:,dupcrit),'rows'); h = histc(uj,1:length(u));
  B = [h(uj) uj B]; C = B(B(:,1)>0,[2 8]); [u ui uj] = unique(C(:,1));
  for j=1:length(u), crows=find(uj==j);
    nref = sum(C(crows,2)>=1 & C(crows,2)<=4); nalt = sum(C(crows,2)>=65 & C(crows,2)<=68); pause=false;
    if P.manual_inspection
      if nalt>=1 & nref>=1
        look(mut,i),fprintf(['\n--(i=%d,  chr%d:%d   AF %0.2f)' repmat('-',1,80) '\n'],i,chr,pos,mut.i_tumor_f(i));
        brows = find(ismember(B(:,2),C(crows,1)));
        margin=5; st=pos-margin; en=pos+margin; E = nan(en-st+1,2*length(brows)+2); E(:,1) = (st:en)';
        for k=1:length(brows), D = Borig(Borig(:,3)==B(brows(k),end),:);
          D=D(D(:,4)>=st & D(:,4)<=en,:); E(D(:,4)-st+1,(k*2)+[0 1]) = D(:,1:2);
        end, E = E(~all(isnan(E(:,2:end)),2),:);
        for k=1:size(E,1), z=E(k,2:2:end-2); z(isnan(z))=[]; if length(unique(z))>1, E(k,end) = -10000; else E(k,end)=0; end; end
        disp(E);disp(B(brows,:));pause=true;
      end
      if pause,keyboard;end
    end
    if nref+1>size(X,1) || nalt+1>size(X,2), X(nref+1,nalt+1,i)=1; else X(nref+1,nalt+1,i)=X(nref+1,nalt+1,i)+1; end
  end
end, fprintf('\n');
BAM.close();

% REPORTING
if isfield(mut,'is_CC_CA') && isfield(mut,'is_CCG_CAG') && isfield(mut,'is_good')
  mut.totreads = squeeze(sum(sum(X,1),2));mut.discdups = squeeze(sum(sum(X(2:end,2:end,:),1),2));
  r2=(mut.is_CCG_CAG | ~mut.is_CC_CA);
  a1=hist2d_fast(double(mut.is_CC_CA),double(mut.discdups>=1),0,1,0,1);
  a2=hist2d_fast(double(mut.is_CC_CA),double(mut.discdups>=2),0,1,0,1);
  a3=hist2d_fast(double(mut.is_CC_CA),double(mut.discdups>=3),0,1,0,1);
  a4=hist2d_fast(double(mut.is_CCG_CAG(r2)),double(mut.discdups(r2)>=1),0,1,0,1);
  a5=hist2d_fast(double(mut.is_CCG_CAG(r2)),double(mut.discdups(r2)>=2),0,1,0,1);
  a6=hist2d_fast(double(mut.is_CCG_CAG(r2)),double(mut.discdups(r2)>=3),0,1,0,1);
  a7=hist2d_fast(double(mut.is_good),double(mut.discdups>=1),0,1,0,1);
  a8=hist2d_fast(double(mut.is_good),double(mut.discdups>=2),0,1,0,1);
  a9=hist2d_fast(double(mut.is_good),double(mut.discdups>=3),0,1,0,1);
  N=[a1(2,2)/(a1(2,2)+a1(2/1)) a2(2,2)/(a2(2,2)+a2(2/1)) a3(2,2)/(a3(2,2)+a3(2/1))
     a4(2,2)/(a4(2,2)+a4(2/1)) a5(2,2)/(a5(2,2)+a5(2/1)) a6(2,2)/(a6(2,2)+a6(2/1))
     a7(2,2)/(a7(2,2)+a7(2/1)) a8(2,2)/(a8(2,2)+a8(2/1)) a9(2,2)/(a9(2,2)+a9(2/1))];
  E=[a1(1,2)/a1(1,1)\(a1(2,2)/a1(2/1)) a2(1,2)/a2(1,1)\(a2(2,2)/a2(2/1)) a3(1,2)/a3(1,1)\(a3(2,2)/a3(2/1))
     a4(1,2)/a4(1,1)\(a4(2,2)/a4(2/1)) a5(1,2)/a5(1,1)\(a5(2,2)/a5(2/1)) a6(1,2)/a6(1,1)\(a6(2,2)/a6(2/1))
     a7(1,2)/a7(1,1)\(a7(2,2)/a7(2/1)) a8(1,2)/a8(1,1)\(a8(2,2)/a8(2/1)) a9(1,2)/a9(1,1)\(a9(2,2)/a9(2/1))];
  F=[fisher_exact_test(a1) fisher_exact_test(a2) fisher_exact_test(a3)
     fisher_exact_test(a4) fisher_exact_test(a5) fisher_exact_test(a6)
     fisher_exact_test(a7) fisher_exact_test(a8) fisher_exact_test(a9)];
  R=[];R.test={'CC->CA','CCG->CAG','good'}';R.frac1=N(:,1);R.frac2=N(:,2);R.frac3=N(:,3);
  R.eff1=E(:,1);R.eff2=E(:,2);R.eff3=E(:,3);R.pval1=F(:,1);R.pval2=F(:,2);R.pval3=F(:,3);
  pr(R)
end
