function [psi_jk] = wavel_trf(j, k, t)
%create wavelet trancform

psi_jk = 1/sqrt(2^j)*basic_wavel((t-(2^j)*k)/2^j);
end