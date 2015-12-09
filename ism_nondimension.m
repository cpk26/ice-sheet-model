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
if ~isfield(ps,'u'), ps.u = 1/pd.ty; end              %Velocity scale
if ~isfield(ps,'t'), ps.t = ps.x/ps.u; end              %Time scale
if ~isfield(ps,'B'), ps.B = pd.B; end                   %Ice Stiffness parameter
if ~isfield(ps,'A'), ps.A = pd.A; end 
if ~isfield(ps,'vis_i'),...                             %Ice viscosity 
    ps.vis_i = ps.B/(2*(ps.u/ps.x)^(1-1/pd.n_Glen)); end;  
if ~isfield(ps,'C'),...                                 %Basal slipperiness
    ps.C = (ps.vis_i*ps.z)/(ps.x^2); end          


pp = struct;
pp.n_Glen = pd.n_Glen;                                  
pp.u = ps.u;
pp.A = pd.A;
pp.g = pd.g;
pp.rho_i = pd.rho_i;
pp.n_rp = pd.n_rp;

pp.c1 = (-2*pp.A*(pp.g*pp.rho_i)^pp.n_Glen)/(pp.n_Glen+2);          %SIA
pp.c2 = ps.z^(pp.n_Glen+1) * (ps.z/ps.x)^pp.n_Glen * 1/ps.u;
pp.c3 = pp.c1*pp.c2;

pp.c4 = pd.B/ps.B;                                                  %SSTREAM
pp.c5 = pd.rho_i*pd.g*ps.z/(ps.B*(ps.u/ps.x)^(1/3));



end

