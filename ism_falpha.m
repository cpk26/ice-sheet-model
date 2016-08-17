function [ F ] = ism_falpha(alpha,U,nEff_lyrs,vv,aa,pp,gg,oo )
%ism_falpha Calculate F_alpha integrals from Arthern, 2015
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   F     F integral

if ~isfield(oo,'nl'), oo.nl = 50; end                   %%Number of layers, must even so vl+1 is odd.
nl = oo.nl;      

Frun = zeros(gg.nha,1);                 %Value in each layer
sp = gg.S_h *aa.h(:)/nl;                %Depth of each layer

u = U(1:gg.nua); u_h = gg.c_uh*u;       %Setup velocity,topographic parameters
v = U(gg.nua+1:end); v_h = gg.c_vh*v;
s = gg.S_h*aa.s(:);


for k =[0:nl]

%% Viscosity of Current Layer
tmpz = gg.S_h*aa.b(:) + (k)*sp;

nEff_l = gg.c_ch*(gg.c_hc*nEff_lyrs(:,k+1));
F_l = (nEff_l).^(-1) .* ((s-tmpz)./(gg.S_h *aa.h(:))).^alpha; 


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
    
F = (sp/3) .* Frun;


end

