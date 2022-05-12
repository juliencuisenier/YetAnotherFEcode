function [FRF] = frf_substructuring(PrimalSub,Fs,DOFs,fmax,fn)
%FRF_SUBSTRUCTURING Summary of this function goes here
%   Calculate the FRF for a given force applied on each substructure and
%   plot the module and phase of given DOFs
%   Fs : cell of length PrimalSub.nSubs
%   DOFs : array of DOFs whose we want to know the FRF
%   fmax : scalar
%   fn : scalar, length of the frequency array

FRF = zeros(length(PrimalSub.globalFreeDOFs),fn);

F = L_to_global(PrimalSub,Fs);

Fc = PrimalSub.constrain_vector(F);

freq = 0:fmax/(fn-1):fmax;
om = 2*pi*freq;

for i=1:fn
    G = -om(i)^2*PrimalSub.DATA.Mc +1i*om(i)*PrimalSub.DATA.Cc +PrimalSub.DATA.Kc;
    FRF(:,i) = G\Fc;
end

for xi=DOFs
   figure
   hold on
   subplot(2,1,1);
   plot(freq,abs(FRF(xi,:)));
   title(strcat('Module of the DOF',num2str(xi),"'s FRF"));
   xlabel('Frequency');
   subplot(2,1,2);
   plot(freq,unwrap(angle(FRF(xi,:))))
   title(strcat('Phase of the DOF',num2str(xi),"'s FRF [rad]"));
   xlabel('Frequency');
end




end

