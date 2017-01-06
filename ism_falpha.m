function [ F ] = ism_falpha(alpha,uv,nEff_lyrs,vv,aa,pp,gg,oo )
%ism_falpha Calculate F_alpha integrals from Arthern, 2015
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   F     F integral

nl = oo.nl;      

Frun = zeros(gg.nha,1);                 %Value in each layer
sp = gg.S_h *aa.h(:)/nl;                %Depth of each layer

u = uv(1:gg.nua);                   %Setup velocity,topographic parameters
v = uv(gg.nua+1:end); 
s = gg.S_h*aa.s(:);
b = gg.S_h*aa.b(:);
h = gg.S_h*aa.h(:);

parfor k =[0:nl]

%% Viscosity of Current Layer
tmpz = b + (k)*sp;

nEff_l = (nEff_lyrs(:,k+1));

F_l = (nEff_l).^(-1) .* ((s-tmpz)./(h)).^alpha; 


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

