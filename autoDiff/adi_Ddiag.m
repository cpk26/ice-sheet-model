


X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; speye(gg.nua,gg.nua) sparse(gg.nua,gg.nva); sparse(gg.nva,gg.nua) speye(gg.nva,gg.nva)];
X2 = [gg.dh_x gg.dh_x gg.duh_y speye(gg.nua,gg.nua) sparse(gg.nva,gg.nua)'; gg.dh_y -gg.dh_y gg.dvh_x sparse(gg.nua,gg.nva)' speye(gg.nva,gg.nva)];

mm = struct();
mm.X = X;
mm.X2 = X2;

% n_d = speye(gg.nha,gg.nha);
% h_d = speye(gg.nha,gg.nha);
% Cu_d = speye(gg.nua,gg.nua);
% Cv_d = speye(gg.nva,gg.nva);
% 
% D = blkdiag(3*n_d*h_d, n_d*h_d, n_d*h_d, -pp.c3*Cu_d, -pp.c3*Cv_d);
% mm.D = D;

C = adigatorCreateDerivInput([gg.nha,1], 'C');
nEff = adigatorCreateAuxInput([gg.nha,1]);


adigator('ism_dim_Ddiag', {C,nEff,aa,pp,gg,oo,mm},'myderiv')





%Testing Routine
% U = [u;v];
% z = gg.S_h*aa.s(:);
% z = gg.S_h*aa.b(:);
% C = gg.S_h*aa.C(:);
% ism_visc(z,U,nEff,C,aa,pp,gg,oo);