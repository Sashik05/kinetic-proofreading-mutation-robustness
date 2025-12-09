%%%% variation of f_c/f_min with K(T);
k_old = 1;
T = [1:50];
k_new = (1/293).*((T+293).*exp(-8.4225e03.*((1./(T+293))-(1/293))));
for i = 1:length(k_new)
[~, m0, ~, ~] = err_fraction(k_old);
f_c(i) = error_calculator(m0, k_new(i));
f_min = error_calculator(m0, k_old);
end

del_f = f_c./f_min;
figure
plot(T, del_f, '--')
xlabel('\Delta T = T_{new} - T_{old}')
ylabel('\Delta f = f_{increased}/f_{min}')