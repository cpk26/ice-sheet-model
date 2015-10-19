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

nNodes = nI*nJ; %Number of nodes
x = linspace(xl, xr, nI); %Array of xCoords
y = linspace(yb, yt, nJ); %Array of yCoords
[xx, yy] = meshgrid(x, y); %Matrices of coordinates

gg.nNodes = nNodes;
gg.nI = nI;
gg.nJ = nJ;
gg.x = x;
gg.y = y;
gg.xx = xx;
gg.yy = yy;

end
