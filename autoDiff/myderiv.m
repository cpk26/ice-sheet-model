% This code was generated using ADiGator version 1.2
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function cst = myderiv(U,vv,aa,pp,gg,oo)
global ADiGator_myderiv
if isempty(ADiGator_myderiv); ADiGator_LoadData(); end
Gator1Data = ADiGator_myderiv.myderiv.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %% Inversion cost function
%User Line: % Inputs:
%User Line: %   vv      struct containing initial solution variables
%User Line: %   aa      prescribed fields, including inputs and boundary conditions
%User Line: %   oo      options
%User Line: % Outputs:
%User Line: %   cst     Inversion cost
cada1f1 = 1:gg.nua;
u.dU = U.dU(Gator1Data.Index1);
u.f = U.f(cada1f1);
%User Line: u = U(1:gg.nua);
cada1f1 = gg.nua + 1;
cada1f2 = length(U.f);
cada1f3 = cada1f1:cada1f2;
v.dU = U.dU(Gator1Data.Index2);
v.f = U.f(cada1f3);
%User Line: v = U(gg.nua+1:end);
%User Line: %Velocity Misfit
cadaconditional1 = strcmp(oo.inv_cst,'abs');
%User Line: cadaconditional1 = strcmp(oo.inv_cst,'abs');
cadaconditional2 = strcmp(oo.inv_cst,'rel');
%User Line: cadaconditional2 = strcmp(oo.inv_cst,'rel');
    cada1f1 = pp.c9*0.5;
    cada1f2dU = u.dU;
    cada1f2 = u.f - aa.u;
    cada1f3dU = 2.*cada1f2(:).^(2-1).*cada1f2dU;
    cada1f3 = cada1f2.^2;
    cada1td1 = sparse(Gator1Data.Index3,Gator1Data.Index4,cada1f3dU,5700,5700);
    cada1td1 = gg.c_uh*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f4dU = full(cada1td1(Gator1Data.Index5));
    cada1f4 = gg.c_uh*cada1f3;
    cada1f5dU = v.dU;
    cada1f5 = v.f - aa.v;
    cada1f6dU = 2.*cada1f5(:).^(2-1).*cada1f5dU;
    cada1f6 = cada1f5.^2;
    cada1td1 = sparse(Gator1Data.Index6,Gator1Data.Index7,cada1f6dU,5700,5700);
    cada1td1 = gg.c_vh*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f7dU = full(cada1td1(Gator1Data.Index8));
    cada1f7 = gg.c_vh*cada1f6;
    cada1td1 = zeros(22500,1);
    cada1td1(Gator1Data.Index9) = cada1f4dU;
    cada1td1(Gator1Data.Index10) = cada1td1(Gator1Data.Index10) + cada1f7dU;
    cada1f8dU = cada1td1;
    cada1f8 = cada1f4 + cada1f7;
    cada1td1 = sum(sparse(Gator1Data.Index11,Gator1Data.Index12,cada1f8dU,5625,11400),1);
    cada1f9dU = full(cada1td1(:));
    cada1f9 = sum(cada1f8);
    cada1f10dU = cada1f1.*cada1f9dU;
    cada1f10 = cada1f1*cada1f9;
    cada1f11dU = gg.dx.*cada1f10dU;
    cada1f11 = cada1f10*gg.dx;
    cst.dU = gg.dy.*cada1f11dU;
    cst.f = cada1f11*gg.dy;
    %User Line: cst = pp.c9*0.5*sum( gg.c_uh*(u-aa.u).^2 + gg.c_vh*(v-aa.v).^2 )*gg.dx*gg.dy;
    %User Line: cst = pp.c10*0.5*sum( gg.c_uh*((u-aa.u)./aa.u).^2 + gg.c_vh*((v-aa.v)./aa.v).^2 )*gg.dx*gg.dy;
cst.dU = pp.L_vel.*cst.dU;
cst.f = pp.L_vel*cst.f;
%User Line: cst = pp.L_vel*cst;
%User Line: %Add Cost function enforcing smoothness for 'grid' discretization
cada1f1 = isequal(pp.L_smooth,0);
cada1f2 = not(cada1f1);
cadaconditional1 = and(1,cada1f2);
%User Line: cadaconditional1 = strcmp(oo.Cdisc,'grid') & ~isequal(pp.L_smooth,0);
    %User Line: C1 = (1/pp.x)*gg.dh_x*gg.S_h*vv.acoeff(:).*((gg.c_hu*gg.S_h*ones(gg.nIJ,1)) == 1);
    %User Line: C2 = (1/pp.x)*gg.dh_y*gg.S_h*vv.acoeff(:).*((gg.c_hv*gg.S_h*ones(gg.nIJ,1)) == 1);
    %User Line: tik = 0.5*pp.c10*(sum(C1(:).^2) + sum(C2(:).^2))*gg.dx*gg.dy;
    %User Line: cst = cst + pp.L_smooth*tik;
cst.dU_size = 11400;
cst.dU_location = Gator1Data.Index13;
end


function ADiGator_LoadData()
global ADiGator_myderiv
ADiGator_myderiv = load('myderiv.mat');
return
end