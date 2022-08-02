function u = static_resolution(PrimalSub,F)
%STATIC_RESOLUTION 
%Solve a static problem. Fext is a (1,nSubs) cell

if iscell(F)
    F = L_to_global(PrimalSub,F);
end

Fc = PrimalSub.constrain_vector(F);
uc = PrimalSub.DATA.Kc\Fc;
u = PrimalSub.unconstrain_vector(uc);
end

