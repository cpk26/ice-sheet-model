function [nEff_di] = ism_visc_di(U,nEff,C,aa,pp,gg,oo)
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

vl = 14;                                %Number of layers, must even so vl+1 is odd.
nEffz = zeros(gg.nha,vl+1);             %Viscosity in each layer
sp = gg.S_h*aa.h(:)/vl;                 %Depth of each layer

u = U(1:gg.nua); u_h = gg.c_uh*u;       %Setup velocity,topographic parameters
v = U(gg.nua+1:end); v_h = gg.c_vh*v;
s = gg.S_h*aa.s(:);
h = gg.S_h*aa.h(:);
n = pp.n_Glen;


exx = gg.du_x*u;                         %Calculate Longitudinal Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.dhu_y*u + gg.dhv_x*v);


for k =[0:vl]
tmpz = gg.S_h*aa.b(:) + k*sp;

exz = pp.c11*C.*u_h.*(s-tmpz)./(nEff.*h);
eyz = pp.c11*C.*v_h.*(s-tmpz)./(nEff.*h);
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + (1/4)*exz.^2 + (1/4)*eyz.^2 + pp.n_rp.^2);

nEffz(:,k+1) =  edeff.^((1-n)/n);                       %Effective Viscosity [dimensionless]
end

tmpA = nEffz(:,[1,end]);
tmpB = 2*nEffz(:,(2:2:end-1));
tmpC = 4*nEffz(:,(3:2:end-2));
tmpvec = [tmpA,tmpB,tmpC];
tmpvec = sum(tmpvec,2);
nEff_di = (1./(gg.S_h*aa.h(:))) .*(sp/3) .* tmpvec;


end


