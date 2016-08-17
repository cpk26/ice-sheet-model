function [nEff_di, nEff_lyrs2] = ism_visc_di(U,nEff_lyrs,C,aa,pp,gg,oo)
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
nEff_lyrs2 = zeros(gg.nha,nl+1);
sp = gg.S_h*aa.h(:)/nl;                 %Depth of each layer

u = U(1:gg.nua); u_h = gg.c_uh*u;       %Setup velocity,topographic parameters
v = U(gg.nua+1:end); v_h = gg.c_vh*v;
s = gg.S_h*aa.s(:);
b = gg.S_h*aa.b(:);
h = gg.S_h*aa.h(:);
n = pp.n_Glen;


exx = gg.du_x*u;                         %Calculate Longitudinal Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.dhu_y*u + gg.dhv_x*v);


u_c = gg.c_uc*u;
v_c = gg.c_vc*v;
s_c = gg.S_c*aa.s_c(:);
b_c = gg.S_c*aa.b_c(:);
h_c = s_c-b_c;
sp_c = h_c/nl; 


%% Depth Integrated Viscosity
for k =[0:nl]


%% Viscosity of Current Layer
%tmpz = b(:) + (k)*sp;
tmpz = b_c + (k)*sp_c;

nEff_l = nEff_lyrs(:,k+1);
for ii = [1:15]
%exz = pp.c11*C.*u_h.*(s-tmpz)./(nEff_l.*h);
%eyz = pp.c11*C.*v_h.*(s-tmpz)./(nEff_l.*h);

exz = pp.c11.*u_c.*(gg.c_hc*C(:)).*(s_c-tmpz)./((gg.c_hc*nEff_l).*h_c);         %C-grid
eyz = pp.c11.*v_c.*(gg.c_hc*C(:)).*(s_c-tmpz)./((gg.c_hc*nEff_l).*h_c);

exz = gg.c_ch*exz(:);                %H-grid 
eyz = gg.c_ch*eyz(:);

edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + (1/4)*exz.^2 + (1/4)*eyz.^2 + pp.n_rp.^2);

nEff_l = edeff.^((1-n)/n);
end

nEff_lyrs2(:,k+1) = nEff_l;


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


