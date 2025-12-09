
% function [f_min, m_o, fz, m] = err_fraction(ky)
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
% 
% m = [linspace(0.00000000001,0.002,30000000),linspace(0.002+0.200001990025879*10^-7,1000,1000000)];
% 
% fx = ((m.*kxp)./((kx + m).*(lx + w)))+ (lxp/(lx + w));
% fy = ((m.*kyp)./((ky + m).*(ly + w)))+ (lyp/(ly + w));
% 
% fz = fy./fx;
% f_min = min(fz);
% y = find(fz == f_min);
% m_o = m(y(1));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Expt values

function [f_min, m_o, fz, m] = err_fraction(kx)
kxp = 250;
kyp = kxp;
alpha = 0.0001;
ky = kx/alpha;
lx = 0.2;
ly = lx/alpha;
lxp = 0.01;
lyp = lxp;
w = 0;

m = [linspace(0.00000000001,0.002,30000000),linspace(0.002+0.200001990025879*10^-7,1000,1000000)];

fx = ((m.*kxp)./((kx + m).*(lx + w)))+ (lxp/(lx + w));
fy = ((m.*kyp)./((ky + m).*(ly + w)))+ (lyp/(ly + w));

fz = fy./fx;
f_min = min(fz);
y = find(fz == f_min);
m_o = m(y(1));
end