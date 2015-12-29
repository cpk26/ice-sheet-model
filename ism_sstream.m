function [vv2] = ism_sstream(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables
%   J       Jacobian matrix

numIter = 15;
iter_norm = zeros(numIter,1);
n = pp.n_Glen;                              %Ice Flow 


%% Variables (Non-Dimensionalized)
s = aa.s; s_diag = spdiags(s(:),0,gg.nIJ,gg.nIJ);      %Topography
h = aa.h; h_diag = spdiags(h(:),0,gg.nIJ,gg.nIJ);
Cslip = aa.C; Cslip_diag = spdiags(Cslip(:),0,gg.nIJ,gg.nIJ);

%Use gradient instead of gg.nddx/y 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);        
Sx = Sx(:); Sy = Sy(:); 

u = vv.u;                                       %Initial iterate velocity 
v = vv.v;


for j = 1:numIter

exx = gg.du_x*u;                                %Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.du_y*u + gg.dv_x*v);
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2);

nEff =  edeff.^((1-n)/n);        %Effective Viscosity [dimensionless]
nEff_diag = spdiags(nEff(:),0,gg.nIJ,gg.nIJ);                                  

%% Field equations for velocities


A1 = gg.S_u*gg.dh_x*4*nEff_diag*h_diag*gg.du_x*gg.S_u';     %LHS SSA 
A2 = gg.S_u*gg.dhu_y*nEff_diag*h_diag*gg.du_y*gg.S_u';
A3 = pp.c3*gg.S_u*gg.c_hu*Cslip_diag*gg.c_uh*gg.S_u';
AA = A1 + A2 - A3;


B1 = gg.S_u*gg.dh_x*2*nEff_diag*h_diag*gg.dv_y*gg.S_v';
B2 = gg.S_u*gg.dhu_y*nEff_diag*h_diag*gg.dv_x*gg.S_v';
BB =  B1+B2;

C1 = gg.S_v*gg.dh_y*2*nEff_diag*h_diag*gg.du_x*gg.S_u';
C2 = gg.S_v*gg.dhv_x*nEff_diag*h_diag*gg.du_y*gg.S_u';
CC = C1 + C2;


D1 = gg.S_v*gg.dh_y*4*nEff_diag*h_diag*gg.dv_y*gg.S_v';
D2 = gg.S_v*gg.dhv_x*nEff_diag*h_diag*gg.dv_x*gg.S_v';
D3 = pp.c3*gg.S_v*gg.c_hv*Cslip_diag*gg.c_vh*gg.S_v';
DD = D1 + D2 - D3;

LL = [AA BB; CC DD];

E1 = gg.S_u_perBC;                                          %Add Periodic Boundary Conditions
E4 = gg.S_v_perBC; 
E2 = zeros(size(E1,1), (gg.nJ+1)*gg.nI);
E3 = zeros(size(E4,1), gg.nJ*(gg.nI+1));

EE = [E1 E2; E3 E4];

LL2 = [LL;EE];                                              %Combine

f1a = pp.c4*gg.c_hu*h_diag*Sx;                              %RHS SSA
f1b = pp.c4*gg.c_hv*h_diag*Sy;
f3 = zeros(size(EE,1),1);                                   %Periodic Boundary Conditions

FF = [f1a;f1b];
FF2 = [f1a;f1b;f3];                                         %Combine

U = LL2\FF2;                                                %Solve
u = gg.S_u*U(1:(gg.nI+1)*gg.nJ);
v = gg.S_v*U((gg.nI+1)*gg.nJ+1:end);

iter_norm(j) = norm(FF2-LL2*U,2);
    
end

vv.u = u;
vv.v = v;
vv.iter_norm = iter_norm;

vv2=vv;

end


function [tbx, tby] = slidingLaw()

switch oo.sL
    case 'ismip'
        tbx = C.*u;
        tby = C.*v;
        
    case 'weertman'
        %To be implemented...
        %w = delta*(u.*Sx + v.*Sy);
        %Ub = sqrt(u.^2 + v.^2 + w.^2);
end


end

