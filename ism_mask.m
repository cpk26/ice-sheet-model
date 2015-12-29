function [gg] = ism_mask(gg,dd,oo)
% Assign default parameters and options
% Inputs 
%   gg struct of grid
%   dd struct of topography
%   oo struct of options [optional]
% Outputs
%   gg struct of grid

% 16 October 2015 

A = (dd.h > 0);
B = (dd.mask == 1);
[C1, C2] = gradient(single(A)); C = C1 | C2;
[D1, D2] = gradient(single(B)); D = D1 | D2;
E = dd.mask > 0; E(:,2:end-1) = 0;
F = dd.mask > 0; F(2:end-1,:) = 0;

%Assign nodes as ice margin, ice interior, interior boundary, or exterior to the problem. 
%Nodes belong to only one category.
gg.nin = A & B & ~C & ~D;               
gg.nmgn = A & B & C;                    
gg.nbnd = A & B & D & ~C;
gg.next = ~(gg.nin | gg.nmgn | gg.nbnd);
clear A B C D

%Determine nodes on the bottom row and rightmost column of the u/v grids which periodic BC apply to.
AA = (E & fliplr(E)); gg.unperBC = reshape((gg.c_hu*AA(:) == 1),gg.nJ, gg.nI+1); gg.unperBC(:,1:end-1) = 0; %u-grid
BB = (F & flipud(F)); gg.vnperBC = reshape((gg.c_hv*BB(:) == 1),gg.nJ+1, gg.nI); gg.vnperBC(2:end,:) = 0;  %v-grd

%Sampling Operators
%H-grid
S_h = spdiags([dd.mask(:) > 0],0,gg.nIJ,gg.nIJ); S_h = S_h(any(S_h,2),:);  

% U-grid;
AA = conv2(dd.mask, [1 1]/2) > 0;                                           
S_u = spdiags(AA(:) > 0,0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); S_u = S_u(any(S_u,2),:);

%U-grid nodes which periodic BC apply to. Of each pair, one is positive, the other
%is negative
BB = zeros(gg.nJ, gg.nI+1);                              %First/last column u-grid nodes with periodic BC
for j = 1:gg.nJ; if abs(AA(j,1)) == abs(AA(j,end)); BB(j,1) = 1; end; end;          %First Column only
CC = spdiags(BB(:),0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); DD = CC(any(CC,2),:);    %Reformat, isolate entries
S_u_perBC = DD - circshift(DD,(gg.nI)*gg.nJ,2);                                     %Add corresponding node (-)

% V-grid;
AA = conv2(dd.mask, [1 1]'/2) > 0;                                          
S_v = spdiags(AA(:) > 0, 0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); S_v = S_v(any(S_v,2),:); 

%V-grid nodes which periodic BC apply to. Of each pair, one is positive, the other
%is negative
BB = zeros(gg.nJ+1, gg.nI);             %Top/bottom row v-grid nodes with periodic BC
for j = 1:gg.nI; if abs(AA(1,j)) == abs(AA(end,j)); BB(1,j) = 1; end; end;          %Top row only
CC = spdiags(BB(:),0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); DD = CC(any(CC,2),:);    %Reformat, isolate entries
S_v_perBC = DD - circshift(DD,gg.nJ,2);                                  %Add corresponding node (-)

%C-grid
AA = padarray(dd.mask,[1,1], 'symmetric'); BB = conv2(AA,[1 1; 1 1]/4, 'valid');  %need symmetric padding
S_c = spdiags(BB(:) == 1,[0],(gg.nI+1)*(gg.nJ+1),(gg.nI+1)*(gg.nJ+1)); S_c = S_c(any(S_c,2),:); 

gg.S_h = S_h;
gg.S_u = S_u;
gg.S_u_perBC = S_u_perBC;
gg.S_v = S_v;
gg.S_v_perBC = S_v_perBC;
gg.S_c = S_c;

end
