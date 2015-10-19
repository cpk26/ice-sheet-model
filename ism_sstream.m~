function [res, jac ] = ism_sstream(vv,aa,pp,gg,oo )
% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables
%   J       Jacobian matrix

%Variables/Parameters
u = vv.u;   
v = vv.v;
[exx, eyx] = gradient(u);
[exy, eyy] = gradient(v);
n = pp.n_Glen;  
A = pp.A;
S = aa.S;       
H = aa.H;
C = aa.C;


%Surface Gradients
[Sx,Sy] = gradient(S);

%Effective Viscosity
vEff = 0.5 * A^(-1/n) * e^(1-n)/n;

%Sliding Law
BB = 1;

%Field equations for velocities
TMP1 = 4*H*vEff*exx + 2*H*vEff*eyy; [TMP1,~] = gradient(A);
TMP2 = H*vEff*(exy + eyx); [~,TMP2] = gradient(B);
TMP3 = BB^2 *u;
TMP4 = pp.rho_i*pp.g*aa.h*Sx;







FE1 = sum([TMP1,TMP2,TMP3,TMP4]);





end

