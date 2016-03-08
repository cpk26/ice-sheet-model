function [vv] = ism_cost_jac(vv, pp, gg, oo)
%% Calculate gradients of alpha coefficients
% Inputs:
%   n_x,n_y   number of alpha coefficients in the x,y directions
%   vv      solution variables
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv       solution variables


cJac = zeros(size(vv.acoeff));      %Jacobian for cost function                               


%Tikanov smoothing contribution
a_xx = (pp.c10.^-2)*gg.du_x*((gg.dh_x*gg.S_h*vv.acoeff(:)).*((gg.c_hu*gg.S_h*ones(gg.nIJ,1)) == 1));
a_yy = (pp.c10.^-2)*gg.dv_y*((gg.dh_y*gg.S_h*vv.acoeff(:)).*((gg.c_hv*gg.S_h*ones(gg.nIJ,1)) == 1));
a_lap = a_xx + a_yy;
tikCntrb = -pp.L_smooth*pp.c10*reshape(gg.S_h'*a_lap,gg.nJ,gg.nI)*gg.dx*gg.dy;

DP = (gg.c_uh*(vv.u.*vv.lambda) + gg.c_vh*(vv.v.*vv.mu));           %Dot product of velocities and lagrange multiplier
Cslip = ism_cslip_field(vv, pp, gg, oo); Cslip = gg.S_h*Cslip(:);   %Basal Slipperiness

if isequal(oo.Cdisc, 'dct2')
Nx = pp.acoeff_nx;                  %Number of coefficients (x/y dirs)
Ny = pp.acoeff_ny;
m = 0:1:gg.nI-1; n = 0:1:gg.nJ-1;   %switch to generic domain
[xx,yy] = meshgrid(m,n); 
elseif isequal(oo.Cdisc, 'grid')
Nx = gg.nI;
Ny = gg.nJ;
end

for j = 0:Ny-1;
for k = 0:Nx-1;

if isequal(oo.Cdisc, 'dct2')
if j == 0; a1 = 1/sqrt(gg.nJ); else a1 = sqrt(2/gg.nJ); end 
if k == 0; a2 = 1/sqrt(gg.nI); else a2 = sqrt(2/gg.nI); end

AA = a1*a2*cos(pi*(2*yy+1)*j/(2*gg.nJ)).*cos(pi*(2*xx+1)*k/(2*gg.nI));
AA = gg.S_h*AA(:);
cJac(j+1,k+1) = -pp.c8*sum(AA.*Cslip.*DP*gg.dx*gg.dy);


elseif isequal(oo.Cdisc, 'grid')

AA = zeros(Ny,Nx); AA(j+1,k+1) = vv.acoeff(j+1,k+1);
AA = gg.S_h*AA(:);

cJac(j+1,k+1) = -pp.c8*sum(AA.*Cslip.*DP*gg.dx*gg.dy) + tikCntrb(j+1,k+1);

end

         
end
end

vv.cJac = cJac;

end