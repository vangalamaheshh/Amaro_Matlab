function list_ziggdecon(D,sample,chr,Qs)

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
qfields = {'amp','del','aod','doa'};

qsign = [1,-1,1,-1];

chr0 = find(D.chrn==chr,1,'first');
fprintf('type chr:startp-endp (armlen) startCN-endCN (ampl) [arm_level]\n');
for f = 1:length(qfields)
    fld = qfields{f};
    q = Qs.(fld)(Qs.(fld)(:,1)==chr & Qs.(fld)(:,5)==sample,:);
    for k = 1:size(q,1)
        x = q(k,2) - chr0;
        dx = q(k,3) - q(k,2);
        if qsign(f) > 0
            y = q(k,6);
            dy = q(k,7) - q(k,6);
        else
            y = q(k,7);
            dy = q(k,6) - q(k,7);
        end
        fprintf('%s chr%d:%d-%d (%0.3f) %0.3f (%0.3f=>%0.3f) [%0.3f]\n',...
            fld,chr,x,x+dx,q(k,8),q(k,4),y,y+dy,q(k,10) );
%!      rectangle('Position',[x,y,dx,dy],'FaceColor',qfcols{f});
    end
end

%hold('on');