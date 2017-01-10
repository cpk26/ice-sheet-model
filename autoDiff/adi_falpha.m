function [ F ] = adi_falpha(alpha,uv,C,vv,aa,pp,gg,oo )
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

u = uv(1:gg.nua); 
u_h = gg.c_uh*u;       %Setup velocity,topographic parameters
v = uv(gg.nua+1:end); 
v_h = gg.c_vh*v;
s = gg.S_h*aa.s(:);
b = gg.S_h*aa.b(:);
h = gg.S_h*aa.h(:);

exx = gg.du_x*u;                         %Calculate Longitudinal Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.dhu_y*u + gg.dhv_x*v);




for k =[0:nl]



%% Viscosity of Current Layer
tmpz = b + (k)*sp;


M  = exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2;
L1 = 0.25*(pp.c11*C.*u_h.*(s-tmpz)./(h)).^2;
L2 = 0.25*(pp.c11*C.*v_h.*(s-tmpz)./(h)).^2;

a1 = (L1 + L2)./(M);
a0 = -1./M;

Q = a1/3;
R = -a0/2;

D = Q.^3 + R.^2;

Sa = R + sqrt(D); 
Sb = ones(size(Sa)); 
Sb(Sa<0) = -1;

Ta = R - sqrt(D);
Tb = ones(size(Ta)); 
Tb(Ta<0) = -1;

S = Sb.*abs(Sa).^(1/3);
T = Tb.*abs(Ta).^(1/3);

nEff_l = S+T;

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

