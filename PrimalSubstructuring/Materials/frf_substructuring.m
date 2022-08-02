function [FRF] = frf_substructuring(PrimalSub,F,DOFs,fmin,fmax,fn)
%FRF_SUBSTRUCTURING Summary of this function goes here
%   Calculate the FRF for a given force applied on each substructure and
%   plot the module and phase of given DOFs
%   F : can be a a cell containing the forces on each substructures
%       or a vector already converted by L_to_global
%   DOFs : array of DOFs whose we want to know the FRF
%   fmax : scalar, the upper bound of the frequency range
%   fn : scalar, length of the frequency array

FRF = zeros(PrimalSub.nDOFglobal,fn);

if iscell(F)
    F = L_to_global(PrimalSub,F);
end

Fc = PrimalSub.constrain_vector(F);

freq = fmin:(fmax-fmin)/(fn-1):fmax;
om = 2*pi*freq;

for i=1:fn
    G = -om(i)^2*PrimalSub.DATA.Mc +1i*om(i)*PrimalSub.DATA.Cc +PrimalSub.DATA.Kc;
    FRF(:,i) = PrimalSub.unconstrain_vector(G\Fc);
end

for xi=DOFs
   figure
   subplot(2,1,1);
   plot(freq,abs(FRF(xi,:)));
   title(strcat('Module of the DOF ',num2str(xi),"'s FRF"));
   xlabel('Frequency');
   set(gca,'yscale','log');
   subplot(2,1,2);
   plot(freq,unwrap(angle(FRF(xi,:))))
   title(strcat('Phase of the DOF ',num2str(xi),"'s FRF [rad]"));
   xlabel('Frequency');
end




end

