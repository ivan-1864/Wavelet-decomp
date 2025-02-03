function [psi] = basic_wavel(t)
% create basic wavelet

% mhat wavelet
psi = (1-t.^2).*exp(-t.^2/2);

% haar wavelet
% if t >= 0 && t < 1/2
%     psi = 1;
% elseif t >= 1/2 && t < 1
%     psi = -1;
% else
%     psi = 0;
% end

end