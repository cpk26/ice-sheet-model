function [ gg ] = ism_label(mask,gg,dd,oo);

if ~isfield(dd,'nfxd'), dd.nfxd = zeros(size(m)); end;
if ~isfield(dd,'nntks'), dd.nntks = (m==3); end;



AA = (dd.h>0) + mask;
[fx,fy] = gradient(AA);
BB = (abs(fx) == 1) | (abs(fy)==1);
nfxd_outer = (BB.*mask);

[fx,fy] = gradient(nfxd_outer);
BB = (abs(fx) == 1) | (abs(fy)==1);
nfxd_inner = (BB.*mask) > 0;

dd.nfxd = nfxd_inner | nfxd_outer;





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






end

