function [ u, v ] = ism_gridtrans( h )
%velInterp transfer values from h-grid onto u and v grids, using nanmean.
%   Detailed explanation goes here

nJ = size(h,1); nI = size(h,2);
u = NaN(nJ,nI+1);
v = NaN(nJ+1,nI);

u(1:nJ,1:nI) = h;
for i=1:nI; u(:,i+1) = nanmean([u(:,i+1),h(:,i)],2); end;

v(1:nJ,1:nI) = h;
for j=1:nJ; v(j+1,:) = nanmean([v(j+1,:);h(j,:)],1); end;


end

