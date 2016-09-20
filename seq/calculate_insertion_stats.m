function [r s] = calculate_insertion_stats(f,r,sw,sw2,method)
% calculate_insertion_stats(f,r,sw,sw2)
%
% parameters:
%
% f = cumulative count of half-mapped forward reads
% r = cumulative count of half-mapped reverse reads
% sw = window (left and right) for computing core statistic (typically ~300)
% sw2 = window for smoothing core statistic (typically ~50)
%
% returns:
%
% r = insertion statistic (sum of leftward and rightward evidence)
% s = sidedness statistic (difference of leftward and rightward evidence)
%
% Mike Lawrence 2010-03-02

if ~exist('method','var'), method=2; end

if method==1        % slightly faster, no extra memory required

  f = f - r;
  f([1:sw end-sw+1:end]) = 0;
  r = 2 * f;
  r(1+sw:end-sw) = r(1+sw:end-sw) - f(1:end-2*sw);
  r(1+sw:end-sw) = r(1+sw:end-sw) - f(1+2*sw:end);
  r([1:2*sw end-2*sw+1:end]) = 0;

  % method 1 doesn't yet compute sidedness

elseif method==2    % requires three extra full-length variables

  % forward: calculate at each point: the SUM in a sw ending at that point
  fw = zeros(length(f),1);
  fw(1+sw:end) = f(1+sw:end) - f(1:end-sw);
  % reverse: calculate at each point: the SUM in a sw beginning at that point
  rw = zeros(length(r),1);
  rw(1:end-sw) = r(1+sw:end) - r(1:end-sw);
  % combined: calculate at each point, dhm = [f(-) + r(+)] - [r(-) + f(+)]
  c = zeros(length(f),1);   % Note: having nan's slows down smooth ~100x
  c(1+sw:end-sw) = fw(1+sw:end-sw) + rw(1+sw:end-sw) - fw(1+2*sw:end) - rw(1:end-2*sw);
  r = c;
  % sidedness: calculate at each point, shm = [f(-) - f(+)] - [r(+) - r(-)]
  c = zeros(length(f),1);   % Note: having nan's slows down smooth ~100x
  c(1+sw:end-sw) = fw(1+sw:end-sw) - fw(1+2*sw:end) - rw(1+sw:end-sw) + rw(1:end-2*sw);
  s = c;
end

% smooth
r = smooth(r,sw2);
s = smooth(s,sw2);

% divide by median
% med = median(abs(r));
% r = r / med;

% divide by constant
r = r / 5;
s = s / 5;
