function us = deformation2substructs(Substructuring,u)
% This function rearranges the deformation of the global system to
% each substructure so that it can be plotted
nSubs = Substructuring.nSubs;
us = {};
for i = 1:nSubs
    Ls = Substructuring.L{i};
    us{i} = Ls * u;
end
end

