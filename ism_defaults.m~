function [pd,oo] = ism_defaults(pd,oo)
% Assign default parameters and options
% Inputs 
%   pd struct of pre-assigned dimensional parameters [optional]
%   oo struct of scales to use for non-dimensionalization [optional]
% Outputs
%   pd struct of parameters
%   oo struct of scaled parameters
%
% 16 October 2015 : taken from nevis_defaults

if nargin<1 || isempty(pd), pd = struct; end
if nargin<2, oo = struct; end

if ~isfield(pd,'ts'), pd.ts = 24*60*60; end                        % default time step [s]
if ~isfield(pd,'td'), pd.td = 24*60*60; end                        % seconds in day [s/d]
if ~isfield(pd,'ty'), pd.ty = 24*60*60*365; end                    % seconds in year [s/y]
if ~isfield(pd,'rho_i'), pd.rho_i = 917; end                       % ice density [kg/m^3]
if ~isfield(pd,'g'), pd.g = 9.81; end                              % gravitational acceleration [m/s^2]
if ~isfield(pd,'n_Glen'), pd.n_Glen = 3; end                       % exponent for ice rheology
if ~isfield(pd,'A'), pd.A = 6.8*10^(-24); end                      % ice rheological parameter [Pa^(-n)/s]
if ~isfield(pd,'K_s'), pd.K_s = 2*pd.A*pd.n_Glen^(-pd.n_Glen); end % ice sheet rheological parameter [Pa^(-n)/s]
if ~isfield(pd,'u_b'), pd.u_b = 0; end                             % initial ice velocity [m/s]

end