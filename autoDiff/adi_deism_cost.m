function [ J ] = adi_deism_cost(uv,C,alpha,vv,aa,pp,gg,oo )


F1 = [];
F2 = [];

if oo.hybrid
F1 = adi_falpha(1,uv,C,vv,aa,pp,gg,oo );
F2 = adi_falpha(2,uv,C,vv,aa,pp,gg,oo );
end

Cb = adi_slidinglaw(alpha,uv,C,F2,vv,aa,pp,gg,oo);


[J] = adi_inv_cost(uv,Cb,alpha,F1,F2,vv,aa,pp,gg, oo);    



end

