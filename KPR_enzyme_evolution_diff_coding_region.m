clear; clc;
if isempty(gcp('nocreate')), parpool; end

%% Parameters
N0 = 1000;
k_old = 1;
% k_new = 42.0098;
k_new = 65;
f_th = 1e-6;
sigma = 1e-3;
T_max = 5*1e03;
run = 16;
N_cap = 2 * N0;
gen_len_values = [1e3, 1e6, 1e8]; % three different genome lengths

%% Trait bounds and optima
m_values = double(vpa(sort(driving_function(f_th, k_old)), 3));
[~, m0, ~, ~] = err_fraction(k_old);
[~, m0_new, ~, ~] = err_fraction(k_new);
m_opt = m0_new;

%% Preallocate result storage for three genome lengths
num_cases = length(gen_len_values);
results = struct();

for g = 1:num_cases
    gen_len = gen_len_values(g);
    fprintf('Running simulation for gen_len = %.1e\n', gen_len);

    mean_m = NaN(run, T_max);
    var_m  = NaN(run, T_max);
    N_hist = NaN(run, T_max);
    beta_grad = NaN(run, T_max);
    FinalPop = cell(run, 1);

    parfor r = 1:run
        M_pop = m0 + sigma * randn(N0, 1);  % Initial population
        local_mean = NaN(1, T_max);
        local_var  = NaN(1, T_max);
        local_hist = NaN(1, T_max);
        local_beta = NaN(1, T_max);

        for t = 1:T_max
            if isempty(M_pop), continue; end
            N = length(M_pop);
            mean_val = mean(M_pop);

            % Environmental switch: old -> new proofreading constant
            if t <= 500
                mean_error = error_calculator(mean_val, k_old);
            else
                mean_error = error_calculator(mean_val, k_new);
            end
%             dm = gen_len * mean_error;
            % Store population stats
            local_mean(t) = mean_val;
            local_var(t)  = var(M_pop);
            local_hist(t) = N;

            % --- Compute β (selection gradient) ---
            eps_fd = 1e-8;
            f_plus  = error_calculator(mean_val + eps_fd, k_new);
            f_minus = error_calculator(mean_val - eps_fd, k_new);
            df_dm = (f_plus - f_minus) / (2 * eps_fd);

            fbar = mean_error;
            p_kill_bar = fbar / f_th;
            S_bar = 1 - p_kill_bar;
            p_off_bar = exp(5*S_bar) / (1 + exp(5*S_bar));
            W_bar = S_bar * p_off_bar;
            dlnW_df = -(1/(f_th*(1 - fbar/f_th))) - (5/f_th) * (1 - p_off_bar);
            local_beta(t) = dlnW_df * df_dm;

            % --- Birth–death–mutation process ---
            death_mask = false(N,1);
            offspring_list = [];

            for i = 1:N
                m_i = M_pop(i);

                % Use old or new proofreading constant
                if t <= 500
                    f_i = error_calculator(m_i, k_old);
                else
                    f_i = error_calculator(m_i, k_new);
                end

                dm = gen_len * f_i;

                % Death check
                p_kill = f_i / f_th;
                if p_kill >= rand || f_i < 0
                    death_mask(i) = true;
                    continue;
                end

                % Three offspring variants
                m1 = m_i - dm; 
                m2 = m_i;
                m3 = m_i + dm;
                if any([m1, m2, m3] <= 0), continue; end

                % Reproduction probability
                a = 5;
                s = 1 - p_kill;
                p_off = exp(a*s) / (1 + exp(a*s));

                % Adaptation-coupled reproduction (slow evolution → fewer offspring)
%                 eta = abs(dm) / (abs(dm) + 1);
%                 p_eff = eta * p_off;

                % Determine offspring number (same rule)
                if (p_off - 0.25) > rand
                    n_off = 2;
                else
                    n_off = 1;
                end

                % Randomly choose offspring
                r_pick = rand;
                if n_off == 1
                    if r_pick < 1/3
                        chosen_off = m1;
                    elseif r_pick < 2/3 && r_pick > 1/3
                        chosen_off = m2;
                    elseif r_pick > 2/3 && r_pick < 1
                        chosen_off = m3;
                    end
                else
                    if r_pick < 1/3
                        chosen_off = [m1; m3];
                    elseif r_pick < 2/3 && r_pick > 1/3
                        chosen_off = [m1; m2];
                    elseif r_pick > 2/3 && r_pick < 1
                        chosen_off = [m2; m3];
                    end
                end

                offspring_list = [offspring_list; chosen_off];
                death_mask(i) = true;
            end

            % Apply deaths & births
            M_pop = M_pop(~death_mask);
            if ~isempty(offspring_list)
                M_pop = [M_pop; offspring_list];
            end

            % Cap population after update
            if length(M_pop) > N_cap
                M_pop = M_pop(randi(length(M_pop), N_cap, 1));
            end
        end

        mean_m(r,:) = local_mean;
        var_m(r,:)  = local_var;
        N_hist(r,:) = local_hist;
        beta_grad(r,:) = local_beta;
        FinalPop{r} = M_pop;
    end

    %% Aggregate stats
    trait_m = mean(mean_m, 1, 'omitnan');
    var_trait = mean(var_m, 1, 'omitnan');
    N_hist_mean = mean(N_hist, 1, 'omitnan');
    delta_trait = [NaN, diff(trait_m)];
    results(g).trait_m = trait_m;
    results(g).N_hist = N_hist_mean;
    results(g).rate = delta_trait;
    results(g).gen_len = gen_len;
end

%% ======= Visualization =======
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultLineLineWidth', 1.8);
x = 1:T_max;
colors = lines(num_cases);

figure('Color','w');
subplot(3,1,1); hold on;
for g = 1:num_cases
    plot(x, results(g).N_hist, 'Color', colors(g,:), 'DisplayName', sprintf('l=%.0e', results(g).gen_len));
end
xlabel('Generation'); ylabel('Population size');
title('Population growth under different genome lengths');
legend show; grid on;

subplot(3,1,2); hold on;
for g = 1:num_cases
    plot(x, results(g).trait_m, 'Color', colors(g,:), 'DisplayName', sprintf('l=%.0e', results(g).gen_len));
end
yline(m_opt, '--k', 'Optimum m');
xlabel('Generation'); ylabel('Mean trait (m)');
title('Trait evolution across genome lengths');
legend show; grid on;

subplot(3,1,3); hold on;
for g = 1:num_cases
    plot(x, results(g).rate, 'Color', colors(g,:), 'DisplayName', sprintf('l=%.0e', results(g).gen_len));
end
xlabel('Generation'); ylabel('\Delta m');
title('Adaptation rate');
legend show; grid on;
