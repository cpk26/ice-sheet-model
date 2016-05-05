function [gg] = ism_grid(nI,nJ,xl,xr,yb,yt,oo)
% Setup Grid for model domain
% inputs
%   nI,nJ size of grid [x,y]
%   xl,xr,yb,yt left and right, bottom and top boundary coordinates
%   oo [optional] option structure
% outputs
%   gg grid struct [see below for contents]
%

nIJ = nI*nJ;            %Number of nodes

Lx = abs(xr-xl);        %Domain dimensions (u,v grids)
Ly = abs(yt-yb);
dx = abs(Lx)/(nI-1);    %Grid spacing
dy = abs(Ly)/(nJ-1); 

xl_h = xl + dx/2;       %min x coordinate (h-grid)
xr_h = xr-dx/2;         %max x coordinate
yb_h = yb + dy/2;       %min y coordinate
yt_h = yt - dy/2;       %max y coordinate

Lx_h = abs(xr_h-xl_h);   %Domain dimensions (h grid)
Ly_h = abs(yt_h-yb_h);   %Domain dimensions (h grid)

x = linspace(xl_h,xr_h,nI);  %array of x coords (h-grid)
y = linspace(yt_h,yb_h,nJ);  %array of y coords

[xx, yy] = meshgrid(x, y); %Matrices of coordinates

x_u = linspace(xl,xr,nI+1); [xx_u,yy_u] = meshgrid(x_u, y);
y_v = linspace(yt,yb,nJ+1); [xx_v,yy_v] = meshgrid(x, y_v);



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

du_y = kron(speye(nI+1),dy_ij_per(nJ)); %derivative of u in y direction from u-grid onto c-grid 
dv_x = kron(dx_ij_per(nI),speye(nJ+1)); %derivative of v in x direction from v-grid onto c-grid 

dh_x = kron(dx_ij_per(nI),speye(nJ)); %derivative of h in x direction from h-grid onto u-grid
dh_y = kron(speye(nI),dy_ij_per(nJ)); %derivative of h in y direction from h-grid onto v-grid


gg.nIJ = nIJ; %Grid Details
gg.nI = nI;
gg.nJ = nJ;
gg.nu = gg.nJ*(gg.nI+1);
gg.nv = (gg.nJ+1)*gg.nI;
gg.nc = (gg.nJ+1)*(gg.nI+1);
gg.x = x;
gg.y = y;
gg.xx = xx;
gg.yy = yy;
gg.xx_u = xx_u;
gg.yy_u = yy_u;
gg.xx_v = xx_v;
gg.yy_v = yy_v;
gg.dx = dx;
gg.dy = dy;
gg.Lx = Lx;
gg.Ly = Ly;
gg.Lx_h = Lx_h;
gg.Ly_h = Ly_h;


gg.du_x = du_x; %Finite Difference Operators
gg.dv_y = dv_y;
gg.du_y = du_y;
gg.dv_x = dv_x;
gg.dh_x = dh_x;
gg.dh_y = dh_y;


gg.c_ch = c_ch; %Centering (Interpolation) Operators
gg.c_vu = c_vu;
gg.c_uv = c_uv;
gg.c_uh = c_uh; 
gg.c_vh = c_vh; 
gg.c_hu = c_hu;
gg.c_hv = c_hv;


end

