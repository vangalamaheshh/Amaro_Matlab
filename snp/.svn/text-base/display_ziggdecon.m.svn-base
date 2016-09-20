function display_ziggdecon(D,sample,chr,Qs)

if ~exist('Qs','var') || isempty(Qs)
    Qs = D.Qs;
end

if ~isfield(Qs,'aod')
    aod = Qs.amp(:,4)>0 & Qs.amp(:,7)<=0;
    Qs.aod = Qs.amp(aod,:);
    Qs.amp = Qs.amp(~aod,:);
end

if ~isfield(Qs,'doa')
    doa = Qs.del(:,4) < 0 & Qs.del(:,7) >= 0;
    Qs.doa = Qs.del(doa,:);
    Qs.del = Qs.del(~doa,:);
end

y = D.dat(D.chrn==chr,sample);
x = 1:length(y);

plot(x,y,'k','LineWidth',2);
set(gca,'Ylim',get(gca,'YLim')+[-1,1]);
title(sprintf('sample %d / chromosome %d',sample,chr));
xlabel('marker index');
ylabel('copy number');
set(gca,'TickDir','out');
set(gca,'XLim',[1,length(y)]);
qfields = {'amp','del','aod','doa'};
qfcols = {[1,0,0],[0,0,1],[1,.5,.5],[.5,.5,1]};
qsign = [1,-1,1,-1];

chr0 = find(D.chrn==chr,1,'first');
for f = 1:length(qfields)
    fld = qfields{f};
    q = Qs.(fld)(Qs.(fld)(:,1)==chr & Qs.(fld)(:,5)==sample,:);
    for k = 1:size(q,1)
        cx = q(k,2) - chr0;
        dx = q(k,3) - q(k,2);
        if qsign(f) > 0
            cy = q(k,6);
            dy = q(k,7) - q(k,6);
        else
            cy = q(k,7);
            dy = q(k,6) - q(k,7);
        end
        if (dx>0)
            rectangle('Position',[cx,cy,dx,dy],'FaceColor',qfcols{f},'EdgeColor',[.5,0,.5]);
        end
    end
end
hold on;plot(x,y,'k','LineWidth',2);

% mark centromere if we can
if isfield(D,'armn')
    cpos = find(D.chrn == chr & D.armn==2,1,'first') - chr0;
    text(cpos,0,'*','Fontweight','bold','FontSize',24,'HorizontalAlignment','center','VerticalAlignment','top');
end

%hold('on');