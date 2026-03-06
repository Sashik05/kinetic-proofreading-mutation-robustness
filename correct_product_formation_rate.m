
function [fx] = correct_product_formation_rate(m, kx)
kxp = 250;            
lx = 0.2;             
lxp = 0.01;           
w = 0;                
% --- Range of m values ---
% m = logspace(-10, 10, 5000);

% --- Single kx value ---
% kx = 15;

% --- Compute formation terms ---
fx = ((m .* kxp) ./ ((kx + m) .* (lx + w))) + (lxp / (lx + w));
end
