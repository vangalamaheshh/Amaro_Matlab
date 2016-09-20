function draw_bkpt_diagram(B,lanes,chr,mn,mx,binsize)

numbins = ceil(mx-mn)/binsize;

ct=[];
B.bin = round((B.start - mn)/binsize);
maxlane = max(lanes.lane);
for lane=0:maxlane, fprintf('%d/%d ',lane,maxlane);
  b = B.bin(B.lane==lane);
  ct(:,lane+1) = histc(b,0:numbins);
end, fprintf('\n');

good = {'300V0','30105','302U6','302U8','310NK','30252','302N4','304HM','3052B','314N6','3027E'};
empty = {'302NE','31555'};
mystery = {'304MG'};
alltumor = {'304M9','304ND'};

fc = regexprep(lanes.id,'AAXX.*','');
good_tums = find(lanes.istum & ismember(fc,good));
tums_in_alltum_FC = find(lanes.istum & ismember(fc,alltumor));
norms_in_alltum_FC = find(~lanes.istum & ismember(fc,alltumor));
good_norms = find(~lanes.istum & ismember(fc,good));
mystery_FC = find(ismember(fc,mystery));

ord = [...
good_tums;...
tums_in_alltum_FC;...
norms_in_alltum_FC;...
mystery_FC;...
good_norms;...
];

categ = [...
1*ones(length(good_tums),1);...
2*ones(length(tums_in_alltum_FC),1);...
3*ones(length(norms_in_alltum_FC),1);...
4*ones(length(mystery_FC),1);...
5*ones(length(good_norms),1);...
];

cov = sum(ct);
mincov = 0.3 * median(cov);
overmincov = find(cov>=mincov);
ordgood = ord(ismember(ord,overmincov));
categgood = categ(ismember(ord,overmincov));

% PLOT
ct_norm = ct ./ repmat(cov,numbins+1,1);
imagesc(ct_norm(:,ordgood));
title('coverage (normalized by lane)','fontsize',20);

% Y-AXIS
labels = {}; ticks = [];
for i = mn:2e6:mx
  ticks = [ticks; ((i-mn)/binsize)+1];
  labels = [labels; sprintf('%d',i/1e6)];
end
set(gca,'ytick',ticks,'yticklabel',labels);
ylabel(['chr' num2str(chr) ' position (Mb)'],'fontsize',20);

% X-AXIS
categ_label = {{'tumor lanes','in good FCs'};{'tumor lanes','in mixup FC'};{'normal lanes','in mixup FC'};...
  {'mystery','flowcell'};{'normal lanes','in good FCs'}};
ticks = [0.5];
for i = 1:length(categ_label);
  idx = find(categgood == i);
  tick = max(idx)+0.5;
  ticks = [ticks; tick];
  line([tick tick],[0 numbins],'color',[0 0 0]);
  text((min(idx)+max(idx))/2,numbins*1.05,categ_label{i},'fontsize',12,...
    'horizontalalign','center','clipping','off','verticalalign','middle');
end
set(gca,'xtick',ticks,'xticklabel',[]);
set(gca,'tickdir','out');
