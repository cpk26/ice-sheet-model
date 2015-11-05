function [gg] = ism_mask(gg,dd,oo)
% Assign default parameters and options
% Inputs 
%   gg struct of grid
%   dd struct of topography
%   oo struct of options [optional]
% Outputs
%   gg struct of grid

% 16 October 2015 


A = (dd.H > 0);
B = (dd.mask == 1);
[C1, C2] = gradient(single(A)); C = C1 | C2;
[D1, D2] = gradient(single(B)); D = D1 | D2;

%Assign nodes as ice margin, ice interior, interior boundary, or exterior to the problem. 
%Nodes belong to only one category.
gg.nin = A & B & ~C & ~D;               
gg.nmgn = A & B & C;                    
gg.nbnd = A & B & D & ~C;
gg.next = ~(gg.nin | gg.nmgn | gg.nbnd);



end
