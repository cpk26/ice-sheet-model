function [gg] = ism_nodeLabels(gg,dd,oo)
% Assign default parameters and options
% Inputs 
%   gg struct of grid
%   dd struct of topography
%   oo struct of options [optional]
% Outputs
%   gg struct of grid

% 16 October 2015 

%Buffer mask matrix by 1 node on each side
maskB = zeros(gg.nJ + 2, gg.nI + 2) -Inf;
maskB(2:gg.nJ+1,2:gg.nI+1) = dd.mask;

%Determine ice sheet margin nodes
nmgn = gg.xx.*0;
for jj = 2:gg.nJ-1; for ii = 2:gg.nI-1;
if isequal(maskB(jj,ii),0) || isequal(maskB(jj,ii),-Inf); continue; end
nhgbrs = maskB(jj-1:jj+1,ii-1:ii+1);
if ismember(0,nhgbrs(:)); nmgn(jj-1,ii-1) = 1; end
end; end;

%Assign nodes as ice margin, interior boundary, interior, or exterior
gg.nmgn = logical(nmgn);                                         %Ice Margin nodes. 
gg.nin =  (dd.mask ==  1); gg.nin(gg.nmgn == 1) = 0;     %Interior Nodes
gg.nbnd = isnan(dd.mask); gg.nbnd(gg.nmgn == 1) = 0;    %Boundary Nodes
gg.next = (dd.mask == 0);                               %Out of the ice sheet





end
