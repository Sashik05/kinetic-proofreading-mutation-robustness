% function[m_val] = driving_function(f_target)
% syms m 
% assume(m>0)
% 
% kxp = 100;
% kyp = kxp;
% alpha = 1e-4;
% ky = 10;
% kx = alpha*ky;
% ly = 5;
% beta = 1e-4;
% lx = beta*ly;
% lxp = 0.0097;
% lyp = lxp;
% w = 1e-5;
% 
% fx = ((m*kxp)/((kx + m)*(lx + w)))+ (lxp/(lx + w));
% fy = ((m*kyp)/((ky + m)*(ly + w)))+ (lyp/(ly + w));
% 
% fz = fy/fx;
% 
% equation = fz == f_target;
% m_val = solve(equation, m, [1e-5, 1e-3]);
% % m_val = solve(equation, m, 'IgnoreAnalyticConstraints', true);
% m_val = vpa(m_val, 3);
% end

function [m_val, fz_val] = driving_function(f_target, ky)
    syms m
    assume(m > 0)

    % Parameters
    kxp = 100;
    kyp = kxp;
    alpha = 1e-4;
%     ky = 10;
    kx = alpha * ky;
    ly = 5;
    beta = 1e-4;
    lx = beta * ly;
    lxp = 0.0097;
    lyp = lxp;
    w = 1e-5;

    % Define fx and fy
    fx = ((m * kxp) / ((kx + m) * (lx + w))) + (lxp / (lx + w));
    fy = ((m * kyp) / ((ky + m) * (ly + w))) + (lyp / (ly + w));
    fz = simplify(fy / fx);   % Trait ratio

    % Try symbolic solve
    equation = fz == f_target;
    m_val_sym = vpasolve(equation, m, [1e-10, 1e3]);

    if isempty(m_val_sym)
        warning('No exact symbolic solution found. Using numerical minimization.');

        % Use fminbnd for accurate numeric solution
        fz_fun = matlabFunction(fz, 'Vars', m);
        cost_fun = @(mv) abs(fz_fun(mv) - f_target);
        m_val = fminbnd(cost_fun, 1e-10, 1e3);  % High-precision minimization

        fz_val = fz_fun(m_val);
        m_val = vpa(m_val, 3);       % Increase precision
        fz_val = vpa(fz_val, 3);
    else
        m_val = vpa(m_val_sym, 3);
        fz_val = vpa(subs(fz, m, m_val), 3);
    end
end


