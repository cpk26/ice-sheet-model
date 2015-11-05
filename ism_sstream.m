function [vv2, FE1, FE2] = ism_sstream(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables
%   J       Jacobian matrix

%% Variables (Non-Dimensionalized)
s = aa.s;                                   %Topography
h = aa.h;
%Use gradient instead of gg.nddx/y 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);        
Sx = Sx(:); Sy = Sy(:); h = h(:);

u = vv.u;                                   %Velocities
v = vv.v;             

exx = gg.nddx*u;                             %Strain Rates
eyx = gg.nddy*u;
eyy = gg.nddy*v;
exy = gg.nddx*v;
edeff = sqrt(0.5*(exx.^2 + eyy.^2) + exy.^2);

n = pp.n_Glen;                              %Ice Flow 
B = pp.c4;
C = aa.C;
vEff = B * edeff.^(1-n)/n;                                     

%% Field equations for velocities
TMP1 = 2*h.*vEff.*exx + h.*vEff.*eyy; TMP1 = gg.nddx*TMP1;
TMP2 = h.*vEff.*exy; TMP2 = gg.nddy*TMP2;
[TMP3, ~] = slidingLaw(u,v,C,Sx,Sy,oo);
TMP4 = -pp.c5*h.*Sx;
FE1 = TMP1+TMP2+TMP3-TMP4;


TMP5 = 2*h.*vEff.*eyy + h.*vEff.*exx; TMP5 = gg.nddy*TMP5;
TMP6 = h.*vEff.*eyx; TMP6 = gg.nddx*TMP6;
[~, TMP7] = slidingLaw(u,v,C,Sx,Sy,oo);
TMP8 = -pp.c2*h.*Sy;
FE2 = TMP5+TMP6+TMP7-TMP8;

vv2=vv;
end


function [tbx, tby] = slidingLaw(u,v,C,Sx,Sy,oo)

switch oo.sL
    case 'ismip'
        tbx = C.*u;
        tby = C.*v;
        
    case 'weertman'
        %To be implemented...
        %w = delta*(u.*Sx + v.*Sy);
        %Ub = sqrt(u.^2 + v.^2 + w.^2);
end


end

