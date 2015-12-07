function [gg] = ism_grid(nI,nJ,xl,xr,yb,yt,oo)
% Setup Grid
% inputs
%   nI,nJ size of grid [x,y]
%   xl,xr,yb,yt left and right, bottom and top boundary coordinates
%   oo [optional] option structure
% outputs
%   gg grid struct [see below for contents]
%
% 16 October, 2015: based on nevis_grid

nIJ = nI*nJ; %Number of nodes
x = linspace(xl, xr, nI); %Array of xCoords
y = linspace(yt, yb, nJ); %Array of yCoords
[xx, yy] = meshgrid(x, y); %Matrices of coordinates
dx = abs(xr-xl)/(nI-1); dxgrid = ones(nIJ,1)*dx; %Grid spacing
dy = abs(yt-yb)/(nJ-1); dygrid = ones(nIJ,1)*dy;

% grid point labels [ labeled inline with Matlab indexing convention, so element r,c becomes r+(c-1)*nJ ]
hgrid = (1:nI*nJ)';
ugrid = (1:(nI+1)*nJ)';
vgrid = (1:nI*(nJ+1))';
cgrid = (1:(nI-1)*(nJ-1))';

%One dimension finite difference operators
dx_ij = @(k) (sparse([1:k,1:k]', [1:k,2:k+1]', [-ones(k,1); ones(k,1)]/dx, k,k+1));
dy_ij = @(k) (sparse([1:k,1:k]', [1:k,2:k+1]', [ones(k,1); -ones(k,1)]/dy, k,k+1));

%One dimensional centering operator
c_ij = @(k) (sparse([1:k,1:k]', [1:k,2:k+1]', [ones(k,1); ones(k,1)]/2, k,k+1));


%Two dimensional finite difference operators [u/v grids onto h/c]
du_x = kron(dx_ij(nI),speye(nJ)); %derivative of u in x direction from u-grid onto h-grid
dv_y = kron(speye(nI),dy_ij(nJ)); %derivative of v in y direction from v-grid onto h-grid

du_y = kron(speye(nI+1),dy_ij(nJ-1)); %derivative of u in y direction  from u-grid onto c-grid
dv_x = kron(dx_ij(nI-1), speye(nJ+1)); %derivative of v in x direction from v-grid  onto c-grid


dhu_x = -(du_x)'; %derivative of u in x direction from h-grid onto u-grid
for r = 1:nJ     %Apply periodic boundary conditions
    for c = 1:nI
        tmp = zeros(1,nI*nJ);
        if c == 1; tmp([r, r + (nI-1)*nJ]) = [1,-1]/dx; dh_x(r,:) = tmp; end
        if c == nI; tmp([r, r + (nI-1)*nJ]) = [1,-1]/dx; dh_x(r + nI*nJ,:) = tmp; end
    end
end

dhu_y = -(du_y)'; %derivative of u in y direction from h-grid onto u-grid [in progress]
for r = 1:nJ     %Apply periodic boundary conditions
    for c = 1:nI
        tmp = zeros(1,nI*nJ);
        if c == 1; tmp([r, r + (nI-1)*nJ]) = [1,-1]/dx; dh_x(r,:) = tmp; end
        if c == nI; tmp([r, r + (nI-1)*nJ]) = [1,-1]/dx; dh_x(r + nI*nJ,:) = tmp; end
    end
end

dhv_y = -(dv_y)'; %derivative of v in y direction from h-grid onto v-grid
for r = 1:nJ     %Apply periodic boundary conditions
    for c = 1:nI
        tmp = zeros(1,nI*nJ);
        if r == 1; tmp([1 + nJ*(c-1), nJ*c]) = [1,-1]/dy; dh_y(1+(c-1)*(nJ+1),:) = tmp; end
        if r == nJ; tmp([1 + nJ*(c-1), nJ*c]) = [-1,1]/dy; dh_y(c*(nJ+1),:) = tmp; end
    end
end







%Node Connects ugrid
nconnect = zeros((nI+1)* nJ,4);
for c = 1:nI
    for r = 1:nJ
        nconnect(r+(c-1)*nJ,:) = [ r+(c-2)*nJ, r+(c)*nJ r-1+(c-1)*nJ  r+1+(c-1)*nJ  ]; 
        if c == 1; nconnect(r+(c-1)*nJ,1) = r+(nI-1)*nJ; end
        if c == nI; nconnect(r+(c-1)*nJ,2) = r; end
        if r == 1; nconnect(r+(c-1)*nJ,3) = c*nJ; end
        if r == nJ; nconnect(r+(c-1)*nJ,4) = 1+(c-1)*nJ; end
    end
end



%x/y derivate operators
nddx = sparse([hgrid; hgrid],[nconnect(:,1); nconnect(:,2)],[-(2*dxgrid).^(-1); (2*dxgrid).^(-1)],nIJ,nIJ);
nddy = sparse([hgrid; hgrid],[nconnect(:,3); nconnect(:,4)],[-(2*dygrid).^(-1); (2*dygrid).^(-1)],nIJ,nIJ);


gg.nNodes = nIJ;
gg.nI = nI;
gg.nJ = nJ;
gg.x = x;
gg.y = y;
gg.xx = xx;
gg.yy = yy;
gg.dx = dx;
gg.dy = dy;
gg.nddx = nddx;
gg.nddy = nddy;


end
