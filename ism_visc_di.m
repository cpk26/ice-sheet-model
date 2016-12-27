function [nEff_di, nEff_lyrs] = ism_visc_di(uv,nEff_lyrs,C,aa,pp,gg,oo)
%% Calculate Ice Viscosity 
% Inputs:
%   z       Depth at each grid cell at which to evaluate viscosity
%   U       [u;v] ice velocities in the x,y directions
%   visEff  Current effective viscosity
%   C       Basal Slip
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: Ax = b
%   LHS     Left hand side. 
%   RHS     Right hand side.

if ~isfield(oo,'nl'), oo.nl = 50; end                   %%Number of layers, must even so vl+1 is odd.
nl = oo.nl; 

nEffrun = zeros(gg.nha,1);

u = uv(1:gg.nua);        u_h = gg.c_uh*u;       %Setup velocity,topographic parameters
v = uv(gg.nua+1:end);    v_h = gg.c_vh*v;
h = gg.S_h*aa.h(:);
s = gg.S_h*aa.s(:);
b = gg.S_h*aa.b(:);
sp = gg.S_h*aa.h(:)/nl;                 %Depth of each layer


n = pp.n_Glen;


exx = gg.du_x*u;                         %Calculate Longitudinal Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.dhu_y*u + gg.dhv_x*v);





%% Depth Integrated Viscosity
for k =[0:nl]


%% Viscosity of Current Layer
tmpz = b(:) + (k)*sp;

nEff_l = nEff_lyrs(:,k+1);

for ii = [1:2]
exz = pp.c11*C.*u_h.*(s-tmpz)./(nEff_l.*h);
eyz = pp.c11*C.*v_h.*(s-tmpz)./(nEff_l.*h);

edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + (1/4)*exz.^2 + (1/4)*eyz.^2 + pp.n_rp.^2);

nEff_l = edeff.^((1-n)/n);
end

nEff_lyrs(:,k+1) = nEff_l;

% M  = exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2;
% L1 = 0.25*(pp.c11*C.*u_h.*(s-tmpz)./(h)).^2;
% L2 = 0.25*(pp.c11*C.*v_h.*(s-tmpz)./(h)).^2;
% 
% a1 = (L1 + L2)./(M);
% a0 = -1./M;
% 
% Q = a1/3;
% R = -a0/2;
% 
% D = Q.^3 + R.^2;
% 
% Sa = R + sqrt(D); 
% Sb = ones(size(Sa)); 
% Sb(Sa<0) = -1;
% 
% Ta = R - sqrt(D);
% Tb = ones(size(Ta)); 
% Tb(Ta<0) = -1;
% 
% S = Sb.*abs(Sa).^(1/3);
% T = Tb.*abs(Ta).^(1/3);
% 
% nEff_l = S+T;
% nEff_lyrs(:,k+1) = nEff_l;




%% Running simpsons rule integration
if k==0 
    nEffrun = nEffrun + nEff_l;      %First Layer
elseif k==nl, 
    nEffrun = nEffrun + nEff_l;      %Last Layers
elseif mod(k,2), 
    nEffrun = nEffrun + 4*nEff_l;    %Odd layers
else
    nEffrun = nEffrun + 2*nEff_l;    %Even Layers
end;


end

nEff_di = (1./h).*(sp/3) .* nEffrun;

end


