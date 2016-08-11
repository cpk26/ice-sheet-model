function [ R ] = ism_deism_res(LHS,RHS,U,aa,pp,gg,oo )
%ism_deism_res Calculate solution residual
% Inputs:
%   LHS     Matrix multiplying U
%   RHS     Driving Stress
%   U     Velocities
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   R     R residual


R = RHS - LHS*U;
R = sum(R.^2).^(-0.5);

end

