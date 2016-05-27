function [ps,pp] = ism_nondimension(pd,ps,oo)
% Nondimensionalise parameters in pd with scales in ps
% Inputs 
%   pd struct of dimensional parameters
%   ps struct of scales to use for non-dimensionalization [optional]
%   oo struct of options [optional]
% Outputs
%   ps struct of scales
%   pp struct of scaled parameters
%   sc scales for plot [optional; for compatability with older versions]
%
% 16 October 2015 : taken from nevis_nondimension
  
if nargin<2 || isempty(ps), ps = struct; end 
if nargin<3, oo = struct; end

%% Default scales
if ~isfield(ps,'x'), ps.x = 10^4; end                   %Lateral scale
if ~isfield(ps,'z'), ps.z = 10^3; end                   %Vertical scale
if ~isfield(ps,'e'), ps.e = ps.z/ps.x; end              %Epsilon
if ~isfield(ps,'u'), ps.u = 1/pd.ty; end                %Velocity scale
if ~isfield(ps,'t'), ps.t = ps.x/ps.u; end              %Time scale
if ~isfield(ps,'B'), ps.B = pd.B; end                   %Ice Stiffness parameter
if ~isfield(ps,'A'), ps.A = pd.A; end                   %Rate Factor
if ~isfield(ps,'sigma'), ps.sigma = pd.rho_i*pd.g*ps.z;
if ~isfield(ps,'vis_i'),...                             %Ice viscosity 
    ps.vis_i = 0.5*ps.B*(ps.u/ps.x)^((1-pd.n_Glen)/pd.n_Glen); end; 


pp = struct;
pp.n_Glen = pd.n_Glen;                                  
pp.u = ps.u;
pp.A = pd.A;
pp.B = pd.B;
pp.g = pd.g;
pp.x = ps.x;
pp.z = ps.z;
pp.vis_i = ps.vis_i;
pp.rho_i = pd.rho_i;
pp.n_rp = pd.n_rp/(ps.u/ps.x);
pp.C_rp = pd.C_rp;
pp.U_rp = pd.U_rp/ps.u;
pp.L_vel = pd.L_vel;
pp.L_smooth = pd.L_smooth;  
pp.acoeff_nx = pd.acoeff_nx; 
pp.acoeff_ny = pd.acoeff_ny; 

%pp.c1 = (ps.sigma/ps.u) * (ps.e).^3;                         %SIA
pp.c1 = (pd.rho_i*pd.g*ps.z/ps.u) * (ps.z/ps.x);              %SIA
pp.c2 = ((2*pp.A)/(pp.n_Glen+1))*((pd.rho_i*pd.g).^(pp.n_Glen)) ...
    *(ps.z^(pp.n_Glen+1))*((ps.z/ps.x)^(pp.n_Glen));

pp.c3 = (ps.x^2)/(ps.z*ps.vis_i);                             %SSTREAM
pp.c4 = (ps.x*ps.z*pp.rho_i*pp.g)/(ps.vis_i*ps.u);

pp.c5 = NaN;                               %Inversion
pp.c6 = (ps.vis_i * ps.z)/(ps.x^2);
pp.c7 = ps.u;
pp.c8 = (ps.x^2)*ps.u;
pp.c9 = (ps.x^2)*ps.u^2;
pp.c10 = (ps.x^2);
pp.c11 = (ps.x./ps.vis_i);
pp.c12 = (ps.z*ps.vis_i*ps.u)/(ps.x^2);


end

