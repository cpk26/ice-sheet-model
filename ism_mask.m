function [gg] = ism_mask(m,gg,dd,oo)
% Apply mask to grid. 
% UNDERSTAND: It is important to note that dh_x/dh_y
% and dhv_x/dhu_y are modified so that they equal zero across the margin of
% fixed cells (dirichlet). This is necessary to solve the adjoint eq, and also for free
% slip. It is recommended to use centering ch_x/... for boundary detection.

% Inputs 
%   gg struct of grid
%   dd struct of topography
%   oo struct of options [optional]
% Outputs
%   gg struct of grid

% 16 Feb 2015 

if ~isfield(dd,'nfxd'), dd.nfxd = zeros(size(m)); end;
if ~isfield(dd,'nntks'), dd.nntks = (m==3); end;

%% Categorize nodes (h-grid)
nbnd = bwperim(m>1,4);                      %Boundary Cells
nin = (m == 2) & ~nbnd;                     %Interior Cells
nfxd = dd.nfxd;                             %Direchlet BC Cells
nmgn = zeros(size(m));                      %Ice Margin Cells (sans direchlet bc nodes)
[i,j] = find(nbnd); 
for p=1:numel(i), 
A1 = max(i(p)-1,1); A2 = min(i(p)+1,gg.nJ); A3 = max(j(p)-1,1); A4 = min(j(p)+1,gg.nI);
if ismember(0,dd.h(A1:A2,A3:A4)); nmgn(i(p),j(p))=1;
end,end
nmgn = nmgn & ~nfxd;                       
E = m > 0; E(:,2:end-1) = 0;               %Periodic BC nodes (must be on domain edge)
F = m > 0; F(2:end-1,:) = 0;
AA = (E & fliplr(E)); BB = (F & flipud(F));   
nperbc = (AA | BB) & ~nfxd;
next = ~(m==2);                                %Cells outside of the mask


%Determine nodes on the bottom row and rightmost column of the u/v grids which periodic BC apply to.
%gg.unperBC = reshape((gg.c_hu*AA(:) == 1),gg.nJ, gg.nI+1); gg.unperBC(:,1:end-1) = 0; %u-grid
%gg.vnperBC = reshape((gg.c_hv*BB(:) == 1),gg.nJ+1, gg.nI); gg.vnperBC(2:end,:) = 0;  %v-grd

%% Sampling Operators
%H-grid
S_h = spdiags([m(:) == 2],0,gg.nIJ,gg.nIJ); S_h = S_h(any(S_h,2),:);  

% U-grid;
AA = conv2(single((m == 2) - dd.nntks), [1 1]/2) > 0;                                           
S_u = spdiags(AA(:) > 0,0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); S_u = S_u(any(S_u,2),:);

%u-grid nodes which periodic BC apply to. Of each pair, one is positive, the other
%is negative
BB = zeros(gg.nJ, gg.nI+1);                              %First/last column u-grid nodes with periodic BC
for j = 1:gg.nJ; if and(nperbc(j,1), nperbc(j,end)); BB(j,1) = 1; end; end;          %First Column only
CC = spdiags(BB(:),0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); DD = CC(any(CC,2),:);    %Reformat, isolate entries
S_u_perBC = DD - circshift(DD,(gg.nI)*gg.nJ,2);   %Add corresponding node (-)


% V-grid;
AA = conv2(single((m == 2) - dd.nntks), [1 1]'/2) > 0;                                          
S_v = spdiags(AA(:) > 0, 0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); S_v = S_v(any(S_v,2),:); 

%v-grid nodes which periodic BC apply to. Of each pair, one is positive, the other
%is negative
BB = zeros(gg.nJ+1, gg.nI);             %Top/bottom row v-grid nodes with periodic BC
for j = 1:gg.nI; if and(nperbc(1,j),  nperbc(end,j)); BB(1,j) = 1; end; end;          %Top row only
CC = spdiags(BB(:),0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); DD = CC(any(CC,2),:);    %Reformat, isolate entries
S_v_perBC = DD - circshift(DD,gg.nJ,2);                                  %Add corresponding node (-)

%C-grid
AA = padarray(single(m==2),[1,1], 'symmetric');      %Symmetric padding for mask
AA2 = padarray(single(dd.h>0),[1,1], 'replicate');      %Extend ice mask
BB = conv2(AA,[1 1; 1 1]/4, 'valid');                   %Identify margin of mask
BB2 = conv2(AA2,[1 1; 1 1]/4, 'valid');                 %Identify ice edge margin
CC = BB > 0 & BB2 == 1;
S_c = spdiags(CC(:) == 1,[0],(gg.nI+1)*(gg.nJ+1),(gg.nI+1)*(gg.nJ+1)); S_c = S_c(any(S_c,2),:); 


%% Boundary Conditions

%% Periodic Boundary Nodes (u/v grid)
%Determine whether periodic BC apply within the mask
perBC = any(S_u_perBC(:)) || any(S_v_perBC(:));

if perBC
nperbc_ugrid = reshape(sparse(sum(S_u_perBC)), gg.nJ, gg.nI+1);
nperbc_vgrid = reshape(sparse(sum(S_v_perBC)), gg.nJ+1, gg.nI);

else
nperbc_ugrid = sparse(gg.nJ,gg.nI+1);
nperbc_vgrid = sparse(gg.nJ+1,gg.nI);

end

% Boundary Nodes (u/v) grid
nbnd_ugrid = (gg.dh_x * (m(:) > 1) ~= 0); nbnd_ugrid = reshape(nbnd_ugrid, gg.nJ,gg.nI+1);
nbnd_vgrid = (gg.dh_y * (m(:) > 1) ~= 0); nbnd_vgrid = reshape(nbnd_vgrid, gg.nJ +1,gg.nI);

nbnd_ugrid = logical(nbnd_ugrid + abs(nperbc_ugrid));
nbnd_vgrid = logical(nbnd_vgrid + abs(nperbc_vgrid));

% Margin Boundary Nodes (u/v grid)
if ~all(nmgn(:) == 0)
gg.mgnBC = 1;
nmgn_ugrid = (gg.c_hu*nmgn(:) ~= 0) & nbnd_ugrid(:); nmgn_ugrid = reshape(sparse(nmgn_ugrid), gg.nJ,gg.nI+1);
nmgn_vgrid = (gg.c_hv*nmgn(:) ~= 0) & nbnd_vgrid(:); nmgn_vgrid = reshape(sparse(nmgn_vgrid), gg.nJ +1,gg.nI);
else
nmgn_ugrid = sparse(gg.nJ,gg.nI+1);
nmgn_vgrid = sparse(gg.nJ+1,gg.nI);
end

% Fixed Boundary Nodes (u/v grid)
% nfxd_ugrid = (gg.c_hu*nfxd(:) == 0.5); nfxd_ugrid = reshape(nfxd_ugrid, gg.nJ,gg.nI+1);
% nfxd_vgrid = (gg.c_hv*nfxd(:) == 0.5); nfxd_vgrid = reshape(nfxd_vgrid, gg.nJ+1,gg.nI);

nfxd_ugrid = (gg.c_hu*nfxd(:) >0); nfxd_ugrid = reshape(nfxd_ugrid, gg.nJ,gg.nI+1);
nfxd_vgrid = (gg.c_hv*nfxd(:) >0); nfxd_vgrid = reshape(nfxd_vgrid, gg.nJ+1,gg.nI);


%% Mask Operators

gg.c_ch = S_h*gg.c_ch*S_c';                                                %Centering operators
gg.c_hc = S_c*gg.c_hc*S_h';

gg.c_uc = S_c*gg.c_uc*S_u';
gg.c_vc = S_c*gg.c_vc*S_v';

gg.c_vu = S_u*gg.c_vu*S_v';
gg.c_vh = S_h*gg.c_vh*S_v';
gg.c_uv = S_v*gg.c_uv*S_u';
gg.c_uh = S_h*gg.c_uh*S_u'; 
gg.c_hu = S_u*gg.c_hu*S_h';
gg.c_hv = S_v*gg.c_hv*S_h';

gg.du_x = S_h*gg.du_x*S_u';                                                %Finite Difference Operators
gg.dv_y = S_h*gg.dv_y*S_v';

%% For Free Slip BC and Adjoint method
Mu = ones(gg.nu,1) - (nbnd_ugrid(:) & nfxd_ugrid(:)); Mu = spdiags(Mu, 0, gg.nu,gg.nu);    %masks (force to be zero) dh_x/dh_y
Mv = ones(gg.nv,1) - (nbnd_vgrid(:) & nfxd_vgrid(:)); Mv = spdiags(Mv, 0, gg.nv,gg.nv);    %across the mask boundary at all boundary u/v
                                                                                           %nodes except the ice margin. Note, for fxd nodes
                                                                                           %this is the outermost node
gg.dh_x = S_u*Mu*gg.dh_x*S_h';
gg.dh_y = S_v*Mv*gg.dh_y*S_h';

%gg.dh_x = S_u*gg.dh_x*S_h';                                             %Sans masking at border
%gg.dh_y = S_v*gg.dh_y*S_h';

cnbnd1 = (gg.dv_x*S_v'*S_v*ones(gg.nv,1) ~= 0);                            %c-nodes where dv_x/du_y are across mask boundary
cnbnd2 = (gg.du_y*S_u'*S_u*ones(gg.nu,1) ~= 0);

%% For Free Slip BC and Adjoint method
Mc1 = ones(gg.nc,1) - (cnbnd1); Mc1 = spdiags(Mc1, 0, gg.nc,gg.nc);        %masks (force to be zero) du_y/dv_x
Mc2 = ones(gg.nc,1) - (cnbnd2); Mc2 = spdiags(Mc2, 0, gg.nc,gg.nc);        %across the mask boundary at all c-nodes determined above

gg.dhv_x = gg.c_ch*S_c*Mc1*gg.dv_x*S_v';                                   %derivative of v in x-direction from v grid onto h-grid
gg.dhu_y = gg.c_ch*S_c*Mc2*gg.du_y*S_u';                                   %derivative of u in y-direction from u grid onto h-grid

%gg.dhv_x = gg.c_ch*S_c*gg.dv_x*S_v';                                      %Sans masking at border
%gg.dhu_y = gg.c_ch*S_c*gg.du_y*S_u'; 


gg.dvh_x = gg.c_uv*gg.dh_x ;                                               %derivative of h in x direction from h-grid onto v-grid
gg.duh_y = gg.c_vu*gg.dh_y;                                                %derivative of h in y direction from h-grid onto u-grid



%% Assign variables to output structure
gg.m = m; 

gg.nha = full(sum(S_h(:)));                               %number of active h/u/v grid nodes
gg.nua = full(sum(S_u(:)));
gg.nva = full(sum(S_v(:)));

gg.nbnd = nbnd;                      
gg.nin = nin;                       
gg.nfxd = nfxd;                                 
gg.nmgn = nmgn;                    
gg.nperbc = nperbc;
gg.next = next;  
gg.nntks = dd.nntks;

gg.S_h = double(S_h);
gg.S_u = double(S_u);
gg.S_v = double(S_v);
gg.S_c = double(S_c);

gg.S_u_perBC = S_u_perBC;
gg.S_v_perBC = S_v_perBC;


gg.nperbc_ugrid = nperbc_ugrid;
gg.nperbc_vgrid = nperbc_vgrid;
gg.nbnd_ugrid = nbnd_ugrid;
gg.nbnd_vgrid = nbnd_vgrid;
gg.nmgn_ugrid = nmgn_ugrid;
gg.nmgn_vgrid = nmgn_vgrid;
gg.nfxd_ugrid = nfxd_ugrid;
gg.nfxd_vgrid = nfxd_vgrid;


end
