function C = load_capture_coverage(C,P)
if ~exist('P','var'), P=[]; end
C = load_coverage(C,P);
