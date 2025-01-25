function [psi] = basic_wavel(t)
% create basic wavelet

% mhat wavelet
psi = (1-t^2)*exp(-t^2/2);
end