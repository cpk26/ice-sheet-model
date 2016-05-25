function [ F ] = ism_falpha(alpha,vv,aa,pp,gg,oo )
%ism_falpha Calculate F_alpha integrals from Arthern, 2015
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   F     F integral

vl = 14;                                %Number of layers, must even so vl+1 is odd.
Fz = zeros(gg.nha,vl+1);                %value in each layer
sp = gg.S_h *aa.h(:)/vl;                 %Depth of each layer

for k =[0:vl]
tmpa = gg.S_h * (aa.s(:) - aa.b(:)) + k*sp;
Fz(:,k+1) = (pp.vis_i*vv.nEff).^(-1) .* (tmpa./(gg.S_h *aa.h(:))).^alpha;  
end
    

F = (sp/3) .* sum([Fz(:,[1,end]), 2*Fz(:,[2:2:end-1]), 4*Fz(:,[3:2:end-2])],2);


end

