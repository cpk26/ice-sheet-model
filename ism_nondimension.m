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
if ~isfield(ps,'u_b'), ps.u_b = 60/pd.ty; end           %Basal sliding speed
if ~isfield(ps,'x'), ps.x = 10^4; end                   %Lateral scale
if ~isfield(ps,'z'), ps.z = 10^2; end                   %Vertical scale
if ~isfield(ps,'t'), ps.t = pd.td; end                  %Time scale

pp = struct;
pp.ts = pd.ts/ps.t;                                     %Time step
pp.u_b = pd.u_b/ps.u_b;
end

