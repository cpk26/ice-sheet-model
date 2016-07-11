function [ F ] = ism_falpha(alpha,nEff,vv,aa,pp,gg,oo )
%ism_falpha Calculate F_alpha integrals from Arthern, 2015
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   F     F integral

nl = 50;                                %Number of layers, must even so vl+1 is odd.
Frun = zeros(gg.nha,1);                 %Value in each layer
sp = gg.S_h *aa.h(:)/nl;                %Depth of each layer

for k =[0:nl]
tmpa = gg.S_h * aa.h(:) - k*sp;
F_l = (pp.vis_i*nEff).^(-1) .* (tmpa./(gg.S_h *aa.h(:))).^alpha; 

%% Running simpsons rule integration
if k==0 
    Frun = Frun + F_l;      %First Layer
elseif k==nl, 
    Frun = Frun + F_l;      %Last Layers
elseif mod(k,2), 
    Frun = Frun + 4*F_l;    %Odd layers
else
    Frun = Frun + 2*F_l;    %Even Layers
end;


end
    
F = (pp.z)*(sp/3) .* Frun;


end

