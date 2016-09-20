function T3 = reshape_TCC_to_3N_breakdown(T)
% reshape TCC matrix to "3N breakdown" form

T2 = zeros(6,8);
T2(:,1) = squeeze(sum(sum(T(1:3,:,:),2)));      % A
T2(:,2) = squeeze(sum(sum(T(10:12,:,:),2)));    % T
T2(:,3) = squeeze(sum(sum(T(4:6,[3 7 11 15],:),2)));      % C in CpG
T2(:,4) = squeeze(sum(sum(T(4:6,[13 14 16],:),2)));   % C in TpC but not CpG
T2(:,5) = squeeze(sum(sum(T(4:6,[1 2 4 5 6 8 9 10 12],:),2)));   % other C
T2(:,6) = squeeze(sum(sum(T(7:9,[5 6 7 8],:),2)));   % G in CpG
T2(:,7) = squeeze(sum(sum(T(7:9,[1 9 13],:),2)));   % G in GpA but not CpG
T2(:,8) = squeeze(sum(sum(T(7:9,[2 3 4 10 11 12 14 15 16],:),2)));   % other G

T3 = zeros(3,8);
T3(1,:) = T2(1,:);   % silent
T3(2,:) = sum(T2([2 4 5 6],:));   % missense
T3(3,:) = T2(3,:);   % nonsense
