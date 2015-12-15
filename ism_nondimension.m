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
    ps.vis_i = ps.B/((ps.u/ps.x)^(1-1/pd.n_Glen)); end; 


pp = struct;
pp.n_Glen = pd.n_Glen;                                  
pp.u = ps.u;
pp.A = pd.A;
pp.B = pd.B;
pp.g = pd.g;
pp.rho_i = pd.rho_i;
pp.n_rp = pd.n_rp/(ps.u/ps.x);
pp.u_max = pd.u_max/ps.u;
                      
pp.c1 = ps.e*ps.sigma/(ps.u);                         %SIA
pp.c2 = pp.A*ps.z*(ps.sigma^pd.n_Glen)*(ps.e^(pd.n_Glen))/(2*ps.u);

pp.c3 = (ps.x^2)/(ps.z*ps.vis_i);                                %SSTREAM
pp.c4 = (ps.x*ps.z*pp.rho_i*pp.g)/(ps.vis_i*ps.u);



end

