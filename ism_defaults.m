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
if ~isfield(pd,'ty'), pd.ty = 365*24*60*60; end                    % seconds in year [s/y]
if ~isfield(pd,'rho_i'), pd.rho_i = 917; end                       % ice density [kg/m^3]
if ~isfield(pd,'g'), pd.g = 9.81; end                              % gravitational acceleration [m/s^2]
if ~isfield(pd,'n_Glen'), pd.n_Glen = 3; end                       % exponent for ice rheology
if ~isfield(pd,'E'), pd.E = 1; end                                 % Enhancement Factor, see Cuffey/Paterson 2010
if ~isfield(pd,'A'), pd.A = pd.E*3.5*10^(-25); end                 % ice rheological parameter [Pa^(-n)/s], see Cuffey/Paterson 2010
if ~isfield(pd,'B'), pd.B = pd.A^(-1/pd.n_Glen); end;              % Ice stiffness parameter (associated rate factor)
if ~isfield(pd,'n_rp'), pd.n_rp = 10^-5/(pd.ty); end;              % Effective Viscosity regularization parameter (m/s) (Arthern et al, 2015)
if ~isfield(pd,'C_rp'), pd.C_rp = 10^6; end;                       % Basal Slipperiness regularation parameter for SIA (to avoid vel=Inf)
if ~isfield(pd,'U_rp'), pd.U_rp = .001/pd.ty; end;                    % Velocity regularization parameter for inversion s.t. there are no zero vels
if ~isfield(pd,'mdR'), pd.mdR = 0.8; end;                           % Limit maximum deformational velocity to a fraction of surface velocity in ism_inverse_sia.
                                                                   % This acts to smooth the initial guess of C
if ~isfield(pd,'L_vel'), pd.L_vel = 1; end                         % Velocity mismatch regularization parameter (Tikhonov reg)
if ~isfield(pd,'L_smooth'), pd.L_smooth = 1; end                   % Smoothness cost regularization parameter (Tikhonov reg)
if ~isfield(pd,'acoeff_nx'), pd.acoeff_nx = 45; end                % Number of terms in x,y directions to keep in DCT2
if ~isfield(pd,'acoeff_ny'), pd.acoeff_ny = 45; end         

if ~isfield(oo,'mdR'), oo.mdR = 0.8; end;                           % Limit maximum deformational velocity to a fraction of surface velocity in ism_inverse_sia.

if ~isfield(oo,'pT'), oo.pT = 'forward'; end                       %problem type: 'forward' or 'inverse'
if ~isfield(oo,'inv_iter'), oo.pic_iter = 10; end;                 %number of picard iterations iterations
if ~isfield(oo,'hybrid'), oo.hybrid = 1; end                       %Approximation: 1 for hybrid, else default to SSA
if ~isfield(oo,'nl'), oo.nl = 50; end                              %Number of vertical layers to use for Simpson's Rule
if ~isfield(oo,'inv_meth'), oo.inv_meth = 'AD'; end                %Inversion Method, either LM or AD
if ~isfield(oo,'savePicIter'), oo.savePicIter = 0; end             %Flag to save intermediate arrays in picard iterations
if ~isfield(oo,'inv_opt'), oo.inv_opt = 'dg'; end                  %Inversion acoeff optimization algorithm: 'dg'/'lbfgs': downgradient/lBFGS 
if ~isfield(oo,'inv_cst'), oo.inv_msft = 'abs'; end                %Least square solution using absolute or relative error
if ~isfield(oo,'norm'), oo.norm = 2; end                           %solution norm for testing convergence. see norm()
if ~isfield(oo,'Cdisc'), oo.Cdisc = 'dct2'; end                    %Discritization of C slip coefficient. Either 'dct2' or 'grid'

                                               

end