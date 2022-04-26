function u = static_resolution(self,Fext)
%STATIC_RESOLUTION Summary of this function goes here
%   Detailed explanation goes here

fg = [];

for iSub=1:self.nSubs
    if isempty(fg)
        fg = self.L{iSub}'*Fext{iSub};
    else
        fg = fg + self.L{iSub}'*Fext{iSub};
    end
end

fgc = self.constrain_vector(fg);
uc = self.DATA.Kc\fgc;
u = self.unconstrain_vector(uc);
end

