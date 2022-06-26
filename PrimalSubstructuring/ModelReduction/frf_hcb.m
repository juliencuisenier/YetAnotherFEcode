function [FRF] = frf_hcb(PrimalSub,M_hcb,C_hcb,K_hcb,T_hcb,L_hcb,Fs,DOFs,fmax,fn)
%FRF_SUBSTRUCTURING Summary of this function goes here
%   Calculate the FRF for a given force applied on each substructure and
%   plot the module and phase of given DOFs
%   Fs : cell of length PrimalSub.nSubs
%   DOFs : array of DOFs whose we want to know the FRF
%   fmax : scalar
%   fn : scalar, length of the frequency array

FRF = zeros(PrimalSub.nDOFglobal,fn);

%FRF_hcb = zeros(size(M_hcb,1),fn);

F_hcb = applying_force_hcb(PrimalSub,Fs,T_hcb,L_hcb);

freq = 0:fmax/(fn-1):fmax;
om = 2*pi*freq;

for i=1:fn
    G = -om(i)^2*M_hcb +1i*om(i)*C_hcb +K_hcb;
    %FRF_hcb(:,i) = G\F_hcb;
    FRF_hcb = G\F_hcb;
    FRFs = converter_reducted_vector(PrimalSub,T_hcb,L_hcb,FRF_hcb);
    
    
    FRF(:,i) = L_to_global(PrimalSub,FRFs);
end




for xi=DOFs
   figure
   hold on
   subplot(2,1,1);
   plot(freq,abs(FRF(xi,:)));
   title(strcat('Module of the DOF',num2str(xi),"'s FRF (from reduction)"));
   xlabel('Frequency');
   set(gca,'yscale','log');
   subplot(2,1,2);
   plot(freq,unwrap(angle(FRF(xi,:))))
   title(strcat('Phase of the DOF',num2str(xi),"'s FRF (from reduction) [rad]"));
   xlabel('Frequency');
end



end
