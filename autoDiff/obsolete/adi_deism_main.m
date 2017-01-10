function [uv2] = adi_deism_main(uv1,C,nEff,vv,aa,pp,gg,oo )


[uv_f] = adi_deism_backbone(uv1,C.f,nEff.f,vv,aa,pp,gg,oo );


[LHS] = ism_deism_fieldeq(C.f,nEff.f, aa,pp,gg,oo);              %Field Equations
[LHS_C] = ism_deism_fieldeq(zeros(gg.nha,1),nEff.dacoeff, aa,pp,gg,oo);              %Field Equations
[LHS_nEff] = ism_deism_fieldeq(C.dacoeff,zeros(gg.nha,1), aa,pp,gg,oo);              %Field Equations

dLHS_dacoeff = LHS_C + LHS_nEff;

duv = LHS\(-dLHS_dacoeff*uv_f);

uv2.dacoeff = duv;
uv2.f = uv_f;

end

