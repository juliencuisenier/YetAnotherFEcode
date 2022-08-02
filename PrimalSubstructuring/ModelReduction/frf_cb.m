function [FRF] = frf_cb(PrimalSub,M_cb,C_hb,K_cb,T_cb,L_cb,Fs,DOFs,fmin,fmax,fn)
%FRF_SUBSTRUCTURING Summary of this function goes here
%   Calculate the FRF for a given force applied on each substructure and
%   plot the module and phase of given DOFs
%   Fs : cell of length PrimalSub.nSubs
%   DOFs : array of DOFs whose we want to know the FRF
%   fmax : scalar
%   fn : scalar, length of the frequency array

FRF = zeros(PrimalSub.nDOFglobal,fn);

%FRF_hcb = zeros(size(M_hcb,1),fn);

F_cb = applying_force_cb(PrimalSub,Fs,T_cb,L_cb);

freq = fmin:(fmax-fmin)/(fn-1):fmax;
om = 2*pi*freq;

for i=1:fn
    G = -om(i)^2*M_cb +1i*om(i)*C_hb +K_cb;
    FRF_cb = G\F_cb;
    FRF(:,i) = converter_reducted_vector(PrimalSub,T_cb,L_cb,FRF_cb);
end




for x_i=DOFs
   figure
   hold on
   subplot(2,1,1);
   plot(freq,abs(FRF(x_i,:)));
   title(strcat('Module of the DOF',num2str(x_i),"'s FRF (from reduction)"));
   xlabel('Frequency');
   set(gca,'yscale','log');
   subplot(2,1,2);
   plot(freq,unwrap(angle(FRF(x_i,:))))
   title(strcat('Phase of the DOF',num2str(x_i),"'s FRF (from reduction) [rad]"));
   xlabel('Frequency');
end



end
