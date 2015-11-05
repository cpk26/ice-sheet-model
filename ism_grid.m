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

% grid point labels [ labeled inline with Matlab indexing convention, so element r,c becomes r+(c-1)*nJ ]
ns = (1:nI*nJ)';

%Node Connections w/ reflective boundary conditions [left node, right node, top node, bottom node]
nconnect = zeros(nIJ,4);
for c = 1:nI
    for r = 1:nJ
        nconnect(r+(c-1)*nJ,:) = [ r+(c-2)*nJ, r+(c)*nJ r-1+(c-1)*nJ  r+1+(c-1)*nJ  ]; 
        if c == 1; nconnect(r+(c-1)*nJ,1) = r+(nI-1)*nJ; end
        if c == nI; nconnect(r+(c-1)*nJ,2) = r; end
        if r == 1; nconnect(r+(c-1)*nJ,3) = c*nJ; end
        if r == nJ; nconnect(r+(c-1)*nJ,4) = 1+(c-1)*nJ; end
    end
end

%Grid spacing
dx = abs(xr-xl)/(nI-1); dxgrid = ones(nIJ,1)*dx;
dy = abs(yt-yb)/(nJ-1); dygrid = ones(nIJ,1)*dy;

%x/y derivate operators
nddx = sparse([ns; ns],[nconnect(:,1); nconnect(:,2)],[-(2*dxgrid).^(-1); (2*dxgrid).^(-1)],nIJ,nIJ);
nddy = sparse([ns; ns],[nconnect(:,3); nconnect(:,4)],[-(2*dygrid).^(-1); (2*dygrid).^(-1)],nIJ,nIJ);


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
