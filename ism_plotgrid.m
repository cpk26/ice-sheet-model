function [] = ism_plotgrid( gg )
%ism_plotgrid(gg): Plot grid nodes
%   INPUT:
%       gg: grid with masked nodes.


figure()
imagesc(gg.nbnd);
title('Boundary Nodes');

figure()

imagesc(gg.nmgn);
title('Ice margin Nodes');

end

