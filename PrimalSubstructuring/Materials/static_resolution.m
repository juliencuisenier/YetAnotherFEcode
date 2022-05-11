function u = static_resolution(PrimalSub,Fext)
%STATIC_RESOLUTION Summary of this function goes here
%   Detailed explanation goes here

fg = zeros(PrimalSub.nDOFglobal,1);

for iSub=1:PrimalSub.nSubs
    fg = fg + PrimalSub.L{iSub}'*Fext{iSub};
end

fgc = PrimalSub.constrain_vector(fg);
uc = PrimalSub.DATA.Kc\fgc;
u = PrimalSub.unconstrain_vector(uc);
end

