function [gg] = ism_mask(gg,dd,oo)
% Assign default parameters and options
% Inputs 
%   gg struct of grid
%   dd struct of topography
%   oo struct of options [optional]
% Outputs
%   gg struct of grid

% 16 Feb 2015 

if ~isfield(dd,'nfxd'), dd.nfxd = zeros(size(dd.mask)); end;

%% Assign nodes to categories, relevant for Boundary conditions (h-grid)
gg.nbnd = bwperim(dd.mask,4);                      %Boundary Cells
gg.nin = dd.mask & ~gg.nbnd;                       %Interior Cells
gg.nfxd = dd.nfxd;                                 %Fixed Cells (Direchlet BC)
gg.nmgn = zeros(size(dd.mask));                    %Ice Margin Cells
[i,j] = find(gg.nbnd); 
for p=1:numel(i), 
A1 = max(i(p)-1,1); A2 = min(i(p)+1,gg.nJ); A3 = max(j(p)-1,1); A4 = min(j(p)+1,gg.nI);
if ismember(0,dd.h(A1:A2,A3:A4)); gg.nmgn(i(p),j(p))=1;
end,end
gg.nmgn = gg.nmgn & ~gg.nfxd;                       %Fixed bc take precedence
E = dd.mask > 0; E(:,2:end-1) = 0;                  %Periodic BC nodes (must be on domain edge)
F = dd.mask > 0; F(2:end-1,:) = 0;
AA = (E & fliplr(E)); BB = (F & flipud(F));   
gg.nperbc = AA | BB;
gg.nbndr = gg.nbnd & ~(gg.nmgn | gg.nfxd | gg.nperbc);  %Remnant boundary cells
gg.next = ~dd.mask;

%Determine nodes on the bottom row and rightmost column of the u/v grids which periodic BC apply to.
gg.unperBC = reshape((gg.c_hu*AA(:) == 1),gg.nJ, gg.nI+1); gg.unperBC(:,1:end-1) = 0; %u-grid
gg.vnperBC = reshape((gg.c_hv*BB(:) == 1),gg.nJ+1, gg.nI); gg.vnperBC(2:end,:) = 0;  %v-grd

%% Sampling Operators
%H-grid
S_h = spdiags([dd.mask(:) > 0],0,gg.nIJ,gg.nIJ); S_h = S_h(any(S_h,2),:);  

% U-grid;
AA = conv2(single(dd.mask), [1 1]/2) > 0;                                           
S_u = spdiags(AA(:) > 0,0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); S_u = S_u(any(S_u,2),:);

%u-grid nodes which periodic BC apply to. Of each pair, one is positive, the other
%is negative
BB = zeros(gg.nJ, gg.nI+1);                              %First/last column u-grid nodes with periodic BC
for j = 1:gg.nJ; if and(gg.nperbc(j,1), gg.nperbc(j,end)); BB(j,1) = 1; end; end;          %First Column only
CC = spdiags(BB(:),0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); DD = CC(any(CC,2),:);    %Reformat, isolate entries
S_u_perBC = DD - circshift(DD,(gg.nI)*gg.nJ,2);   %Add corresponding node (-)


% V-grid;
AA = conv2(single(dd.mask), [1 1]'/2) > 0;                                          
S_v = spdiags(AA(:) > 0, 0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); S_v = S_v(any(S_v,2),:); 

%v-grid nodes which periodic BC apply to. Of each pair, one is positive, the other
%is negative
BB = zeros(gg.nJ+1, gg.nI);             %Top/bottom row v-grid nodes with periodic BC
for j = 1:gg.nI; if and(gg.nperbc(1,j),  gg.nperbc(end,j)); BB(1,j) = 1; end; end;          %Top row only
CC = spdiags(BB(:),0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); DD = CC(any(CC,2),:);    %Reformat, isolate entries
S_v_perBC = DD - circshift(DD,gg.nJ,2);                                  %Add corresponding node (-)

%C-grid
AA = padarray(single(dd.mask),[1,1], 'symmetric');      %Symmetric padding for mask
AA2 = padarray(single(dd.h>0),[1,1], 'replicate');       %Extend ice mask
BB = conv2(AA,[1 1; 1 1]/4, 'valid');                   %Identify margin of mask
BB2 = conv2(AA2,[1 1; 1 1]/4, 'valid');                 %Identify ice edge margin
CC = BB > 0 & BB2 == 1;
S_c = spdiags(CC(:) == 1,[0],(gg.nI+1)*(gg.nJ+1),(gg.nI+1)*(gg.nJ+1)); S_c = S_c(any(S_c,2),:); 


%% Boundary Conditions

%% Periodic Boundary Nodes (u/v grid)
%Determine whether whether periodic BC apply within the mask
perBC = any(S_u_perBC(:)) || any(S_v_perBC(:));

if perBC
%List of indices in U velocity vector for u/v and +/- velocity sets
[nperbc_u1ind,~] = find(S_u_perBC' == 1);   %half of the u-grid points with periodic bc, 
[nperbc_u2ind,~] = find(S_u_perBC' == -1);  %corresponding half to the above points

[nperbc_v1ind,~] = find(S_v_perBC' == 1);   %v-grid points with periodic bc, set A
[nperbc_v2ind,~] = find(S_v_perBC' == -1);  %corresponding half to the above points

else
nperbc_u1ind = []; nperbc_u2ind = [];
nperbc_v1ind = []; nperbc_v2ind = [];
end

% Boundary Nodes (u/v) grid

unbnd = ~eq(gg.dh_x * dd.mask(:),0);
vnbnd = ~eq(gg.dh_y * dd.mask(:),0);

A = find(gg.c_hu*dd.mask(:) == 0.5); B = union(nperbc_u1ind,nperbc_u2ind); 
nbnd_uind = union(A,B);
A = find(gg.c_hv*dd.mask(:) == 0.5); B = union(nperbc_v1ind,nperbc_v2ind); 
nbnd_vind = union(A,B);

% Margin Boundary Nodes (u/v grid)
if ~isempty(gg.nmgn)
gg.mgnBC = 1;
unmgn = ~eq(gg.c_hu*gg.nmgn(:), 0) & unbnd;
vnmgn = ~eq(gg.c_hv*gg.nmgn(:),0) & vnbnd;

A = find(gg.c_hu*gg.nmgn(:) == 0.5); nmgn_uind = intersect(A,nbnd_uind);
A = find(gg.c_hv*gg.nmgn(:) == 0.5); nmgn_vind = intersect(A,nbnd_vind);
end

% Fixed Boundary Nodes (u/v grid)
A = find(gg.c_hu*gg.nfxd(:) == 0.5); nfxd_uind = intersect(A,nbnd_uind);
A = find(gg.c_hv*gg.nfxd(:) == 0.5); nfxd_vind = intersect(A,nbnd_vind);

% Remnant Boundary Nodes (u/v grid)
nbndr_uind = setdiff(nbnd_uind, [nmgn_uind;nfxd_uind;nperbc_u1ind;nperbc_u2ind]);
nbndr_vind = setdiff(nbnd_vind, [nmgn_vind;nfxd_vind;nperbc_v1ind;nperbc_v2ind]);

%% Mask Operators

gg.c_ch = S_h*gg.c_ch*S_c';                                             %Centering operators
gg.c_vu = S_u*gg.c_vu*S_v';
gg.c_vh = S_h*gg.c_vh*S_v';
gg.c_uv = S_v*gg.c_uv*S_u';
gg.c_uh = S_h*gg.c_uh*S_u'; 
gg.c_hu = S_u*gg.c_hu*S_h';
gg.c_hv = S_v*gg.c_hv*S_h';

gg.du_x = S_h*gg.du_x*S_u';                                                 %Finite Difference Operators
gg.dv_y = S_h*gg.dv_y*S_v';

Mu = ones(gg.nu,1) - (unbnd - unmgn); Mu = spdiags(Mu, 0, gg.nu,gg.nu);          %masks (force to be zero) dh_x/dh_y
Mv = ones(gg.nv,1) - (vnbnd - vnmgn); Mv = spdiags(Mv, 0, gg.nv,gg.nv);          %across the mask boundary at appropriate u/v nodes

gg.dh_x = S_u*Mu*gg.dh_x*S_h';
gg.dh_y = S_v*Mv*gg.dh_y*S_h';

cnbnd1 = ~eq(gg.dv_x*S_v'*S_v*ones(gg.nv,1),0);                             %c-nodes where dv_x/du_y are across mask boundary (all such c-nodes)
cnbnd2 = ~eq(gg.du_y*S_u'*S_u*ones(gg.nu,1),0);

Mc1 = ones(gg.nc,1) - (cnbnd1); Mc1 = spdiags(Mc1, 0, gg.nc,gg.nc);         %masks (force to be zero) du_y/dv_x
Mc2 = ones(gg.nc,1) - (cnbnd2); Mc2 = spdiags(Mc2, 0, gg.nc,gg.nc);         %across the mask boundary at appropriate c-nodes

gg.dhv_x = gg.c_ch*S_c*Mc1*gg.dv_x*S_v';                                    %derivative of v in x-direction from v grid onto h-grid
gg.dhu_y = gg.c_ch*S_c*Mc2*gg.du_y*S_u';                                    %derivative of u in y-direction from u grid onto h-grid

gg.dvh_x = gg.c_uv*gg.dh_x ;                                                %derivative of h in x direction from h-grid onto v-grid
gg.duh_y = gg.c_vu*gg.dh_y;                                                 %derivative of h in y direction from h-grid onto u-grid



%% Assign variables to Struct

gg.S_h = S_h;
gg.S_u = S_u;
gg.S_v = S_v;
gg.S_c = S_c;

gg.S_u_perBC = S_u_perBC;
gg.S_v_perBC = S_v_perBC;


gg.nbnd_uind = nbnd_uind;
gg.nbnd_vind = nbnd_vind;
gg.nmgn_uind = nmgn_uind;
gg.nmgn_vind = nmgn_vind;
gg.nfxd_uind = nfxd_uind;
gg.nfxd_vind = nfxd_vind;
gg.nbndr_uind = nbndr_uind;
gg.nbndr_vind = nbndr_vind;
gg.nperbc_u1ind = nperbc_u1ind;
gg.nperbc_u2ind = nperbc_u2ind;
gg.nperbc_v1ind = nperbc_v1ind;
gg.nperbc_v2ind = nperbc_v2ind;
end
