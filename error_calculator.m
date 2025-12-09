% function [fz] = error_calculator(m, ky)
% kxp = 100;
% kyp = kxp;
% alpha = 0.0001;
% % ky = 2;
% kx = alpha*ky;
% ly = 5;
% lx = alpha*ly;
% lxp = 0.0097;
% lyp = lxp;
% w = 0;
% fx = ((m*kxp)/((kx + m)*(lx + w)))+ (lxp/(lx + w));
% fy = ((m*kyp)/((ky + m)*(ly + w)))+ (lyp/(ly + w));
% fz = fy/fx;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% expt values
function [fz] = error_calculator(m, kx)
kxp = 250;
kyp = kxp;
alpha = 0.0001;
% ky = 2;
ky = kx/alpha;
lx = 0.2;
ly = lx/alpha;
lxp = 0.01;
lyp = lxp;
w = 0;
fx = ((m*kxp)/((kx + m)*(lx + w)))+ (lxp/(lx + w));
fy = ((m*kyp)/((ky + m)*(ly + w)))+ (lyp/(ly + w));
fz = fy/fx;
end