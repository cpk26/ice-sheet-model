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

if ~isfield(pd,'td'), pd.td = 24*60*60; end                        % seconds in day [s/d]
if ~isfield(pd,'ty'), pd.ty = 24*60*60*365; end                    % seconds in year [s/y]
if ~isfield(pd,'rho_i'), pd.rho_i = 917; end                       % ice density [kg/m^3]
if ~isfield(pd,'g'), pd.g = 9.81; end                              % gravitational acceleration [m/s^2]
if ~isfield(pd,'n_Glen'), pd.n_Glen = 3; end                       % exponent for ice rheology
if ~isfield(pd,'A'), pd.A = 3.5*10^(-25); end                      % ice rheological parameter [Pa^(-n)/s], see Cuffey/Paterson 2010
if ~isfield(pd,'B'), pd.B = 0.5*pd.A^(-1/pd.n_Glen); end;          % Ice stiffness parameter (associated rate factor) 


if ~isfield(oo,'fm'), oo.pT = 1; end                               % problem type: [1] for forward model, [2] for inverse problem
if ~isfield(oo,'sL'), oo.sL = 'ismip'; end                         %Sliding Law: ismip                                           
                                                                   

end