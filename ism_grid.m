function [gg] = ism_grid(nI,nJ,xl,xr,yb,yt,oo)
% Setup Grid
% inputs
%   nI,nJ size of grid [x,y]
%   xl,xr,yb,yt left and right, bottom and top boundary coordinates
%   oo [optional] option structure
% outputs
%   gg grid struct [see below for contents]
%
% Dec 8, 2015: based on nevis_grid

nIJ = nI*nJ; %Number of nodes
x = linspace(xl, xr, nI); %Array of xCoords
y = linspace(yb, yt, nJ); %Array of yCoords
[xx, yy] = meshgrid(x, y); %Matrices of coordinates
dx = abs(xr-xl)/(nI-1);   %Grid spacing
dy = abs(yt-yb)/(nJ-1); 

% grid point labels [ labeled inline with Matlab indexing convention, so element r,c becomes r+(c-1)*nJ ]
hgrid = (1:nI*nJ)';
ugrid = (1:(nI+1)*nJ)';
vgrid = (1:nI*(nJ+1))';
cgrid = (1:(nI+1)*(nJ+1))';

%One dimensional centering operator
c_ij = @(k) (sparse([1:k,1:k]', [1:k,2:k+1]', [ones(k,1); ones(k,1)]/2, k,k+1));
c_ij_per = @(k) (sparse([1:k+1,1:k+1]', [[1:k,1],[k,1:k]]', [ones(k+1,1); ones(k+1,1)]/2, k+1,k));

%One dimension finite difference operators
dx_ij = @(k) (sparse([1:k,1:k]', [1:k,2:k+1]', [-ones(k,1); ones(k,1)]/dx, k,k+1));
dy_ij = @(k) (sparse([1:k,1:k]', [1:k,2:k+1]', [ones(k,1); -ones(k,1)]/dy, k,k+1));

dx_ij_per = @(k) (sparse([1:k+1,1:k+1]', [[1:k,1],[k,1:k]]', [ones(k+1,1); -ones(k+1,1)]/dx, k+1,k));
dy_ij_per = @(k) (sparse([1:k+1,1:k+1]', [[1:k,1],[k,1:k]]', [-ones(k+1,1); ones(k+1,1)]/dy, k+1,k));

%Two dimensional centering operators
c_ch = kron(c_ij(nI), c_ij(nJ));
c_vu = kron(c_ij_per(nI), c_ij(nJ));
c_uv = kron(c_ij(nI),c_ij_per(nJ));
c_uh = kron(c_ij(nI),speye(nJ)); 
c_vh = kron(speye(nI),c_ij(nJ));
c_hu = kron(c_ij_per(nI),speye(nJ));
c_hv = kron(speye(nI),c_ij_per(nJ));

%Two dimensional finite difference operators 
du_x = kron(dx_ij(nI),speye(nJ)); %derivative of u in x direction from u-grid onto h-grid
dv_y = kron(speye(nI),dy_ij(nJ)); %derivative of v in y direction from v-grid onto h-grid

du_y = c_ch*kron(speye(nI+1),dy_ij_per(nJ)); %derivative of u in y direction from u-grid onto h-grid 
dv_x = c_ch*kron(dx_ij_per(nI),speye(nJ+1)); %derivative of v in x direction from v-grid onto h-grid 


dh_x = kron(dx_ij_per(nI),speye(nJ)); %derivative of h in x direction from h-grid onto u-grid
dh_y = kron(speye(nI),dy_ij_per(nJ)); %derivative of h in y direction from h-grid onto v-grid

dhu_y = c_vu*dh_y; %derivative of h in y direction from h-grid onto u-grid
dhv_x = c_uv*dh_x; %derivative of h in x direction from h-grid onto v-grid



gg.nIJ = nIJ; %Grid Details
gg.nI = nI;
gg.nJ = nJ;
gg.x = x;
gg.y = y;
gg.xx = xx;
gg.yy = yy;
gg.dx = dx;
gg.dy = dy;

gg.du_x = du_x; %Finite Difference Operators
gg.dv_y = dv_y;
gg.du_y = du_y;
gg.dv_x = dv_x;
gg.dh_x = dh_x;
gg.dh_y = dh_y;
gg.dhu_y = dhu_y;
gg.dhv_x = dhv_x;


gg.c_ch = c_ch; %Centering (Interpolation) Operators
gg.c_vu = c_vu;
gg.c_uv = c_uv;
gg.c_uh = c_uh; 
gg.c_vh = c_vh; 
gg.c_hu = c_hu;
gg.c_hv = c_hv;


end



% dh_x = -(du_x)'; %derivative of h in x direction from h-grid onto u-grid
% for r = 1:nJ     %Apply periodic boundary conditions
%     for c = 1:nI
%         if c == 1; tmp = zeros(1,nI*nJ); tmp([r, r + (nI-1)*nJ]) = [1,-1]/dx; dh_x(r,:) = tmp; end
%         if c == nI; tmp = zeros(1,nI*nJ); tmp([r, r + (nI-1)*nJ]) = [1,-1]/dx; dh_x(r + nI*nJ,:) = tmp; end
%     end
% end


% dh_y = -(dv_y)'; %derivative of h in y direction from h-grid onto v-grid
% for r = 1:nJ     %Apply periodic boundary conditions
%     for c = 1:nI
%         if r == 1; tmp = zeros(1,nI*nJ); tmp([1 + nJ*(c-1), nJ*c]) = [-1,1]/dy; dh_y(1+(c-1)*(nJ+1),:) = tmp; end
%         if r == nJ; tmp = zeros(1,nI*nJ); tmp([1 + nJ*(c-1), nJ*c]) = [-1,1]/dy; dh_y(c*(nJ+1),:) = tmp; end
%     end
% end
