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
if ~isfield(pd,'rho_w'), pd.rho_w = 1000; end                      % water density [kg/m^3]
if ~isfield(pd,'g'), pd.g = 9.81; end                              % gravitational acceleration [m/s^2]
if ~isfield(pd,'L'), pd.L = 3.35*10^5; end                         % latent heat [J/kg]
if ~isfield(pd,'G'), pd.G = 0.063; end                             % geothermal heat flux [J/s/m^2]
if ~isfield(pd,'n_Glen'), pd.n_Glen = 3; end                       % exponent for ice rheology
if ~isfield(pd,'E'), pd.E = 1; end                                 % Enhancement Factor, see Cuffey/Paterson 2010
if ~isfield(pd,'A'), pd.A = pd.E*3.5*10^(-25); end                 % ice rheological parameter [Pa^(-n)/s], see Cuffey/Paterson 2010
if ~isfield(pd,'B'), pd.B = pd.A^(-1/pd.n_Glen); end;              % Ice stiffness parameter (associated rate factor)
if ~isfield(pd,'p'), pd.p = 1/3; end;                              % Power of effective pressure in weertman sliding law
if ~isfield(pd,'q'), pd.q = 1/3; end;                              % Power of velocity in weertman sliding law
if ~isfield(pd,'lambda_b'), pd.lambda_b = 1; end;                  % Bed roughness length for schoof sliding law

if ~isfield(pd,'n_rp'), pd.n_rp = 10^-5/(pd.ty); end;              % Effective Viscosity regularization parameter (m/s) (Arthern et al, 2015)
if ~isfield(pd,'C_rp'), pd.C_rp = 10^2; end;                       % Basal Slipperiness regularation parameter for SIA (to avoid vel=Inf)
if ~isfield(pd,'U_rp'), pd.U_rp = .01/pd.ty; end;                  % Velocity regularization parameter for inversion s.t. there are no zero vels
if ~isfield(pd,'N_rp'), pd.N_rp = 0.1*pd.g*pd.rho_i; end;          % Effective pressure regularization 
if ~isfield(pd,'mdR'), pd.mdR = 0; end;                            % Limit maximum deformational velocity to a fraction of surface velocity in ism_inverse_sia.
                                                                   % This acts to smooth the initial guess of C
if ~isfield(pd,'L_vel'), pd.L_vel = 1; end                         % Velocity mismatch regularization parameter (Tikhonov reg)
if ~isfield(pd,'L_smooth'), pd.L_smooth = 1; end                   % Smoothness cost regularization parameter (Tikhonov reg)
if ~isfield(pd,'acoeff_nx'), pd.acoeff_nx = 45; end                % Number of terms in x,y directions to keep in DCT2
if ~isfield(pd,'acoeff_ny'), pd.acoeff_ny = 45; end         


if ~isfield(oo,'pT'), oo.pT = 'forward'; end                       %problem type: 'forward' or 'inverse'
if ~isfield(oo,'pic_iter'), oo.pic_iter = 15; end;                 %number of picard iterations iterations
if ~isfield(oo,'pic_tol'), oo.pic_tol = 1e-6; end;                 %termination condition on tolerance of norm changes in picard iter
if ~isfield(oo,'hybrid'), oo.hybrid = 1; end                       %Approximation: 1 for hybrid, else default to SSA
if ~isfield(oo,'slidinglaw'), oo.slidinglaw = 'linear'; end        %linear, weertman, or schoof. 
if ~isfield(oo,'nl'), oo.nl = 50; end                              %Number of vertical layers to use for Simpson's Rule
if ~isfield(oo,'inv_meth'), oo.inv_meth = 'AD'; end                %Inversion Method, either LM or AD
if ~isfield(oo,'adj_iter'), oo.adj_iter = 0; end                   %Flag to save intermediate arrays in picard iterations
if ~isfield(oo,'adjAD_uv'), oo.adjAD_uv = 0; end                   %Flag to calculate adjoint of uv
if ~isfield(oo,'adjAD_alpha'), oo.adjAD_alpha = 0; end             %Flag to calculate adjoint of alpha
if ~isfield(oo,'adjAD_AFPI'), oo.adjAD_afpi = 0; end               %Apply Adjoint fixed point iteration
if ~isfield(oo,'adjAD_AFPI_iter'), oo.adjAD_afpi_iter = 75; end    %Number of AFPI iterations
if ~isfield(oo,'adjAD_AFPI_tol'),oo.adjAD_afpi_tol=oo.pic_tol; end %Tolerance termination condition
if ~isfield(oo,'inv_opt'), oo.inv_opt = 'lbfgs'; end               %Inversion acoeff optimization algorithm: 'gd'/'lbfgs': gradientDescent/lBFGS 
if ~isfield(oo,'inv_cst'), oo.inv_msft = 'abs'; end                %Least square solution using absolute or relative error
if ~isfield(oo,'inv_iter'), oo.inv_iter = 75; end                  %Maximum number of inversion iterations
if ~isfield(oo,'inv_funcEval'), oo.inv_funcEval = 150; end         %Maximum number of function evaluations allowed for minFunc()
if ~isfield(oo,'inv_progTolFrac'), oo.inv_progTolFrac = .01; end   %Progress Tolerance for minFunc() as a fraction of original cost
if ~isfield(oo,'norm'), oo.norm = 2; end                           %solution norm for testing convergence. see norm()
if ~isfield(oo,'Cdisc'), oo.Cdisc = 'grid'; end                    %Discritization of C slip coefficient. Either 'dct2' or 'grid'

                                               

end