clear; clc;

%% Start parallel pool
if isempty(gcp('nocreate')), parpool; end

%% Parameters
N_list = [1000, 10000, 20000];  % three different population sizes
k_old = 1;
k_new = 65;
f_th = 1e-6;
sigma = 1e-3;
T_max = 2*1e03;
run = 16; 
gen_len = 1e07;

%% Trait bounds and optima
m_values = double(vpa(sort(driving_function(f_th, k_old)), 3));
[~, m0, ~, ~] = err_fraction(k_old);
[~, m0_new, ~, ~] = err_fraction(k_new);
m_opt = m0_new;

%% Run simulations for each population size
results = struct();

for n_idx = 1:length(N_list)
    N0 = N_list(n_idx);
    N_cap = 2 * N0;

    %% Preallocate global storage
    mean_m = NaN(run, T_max);
    var_m  = NaN(run, T_max);
    N_hist = NaN(run, T_max);
    beta_grad = NaN(run, T_max);
    FinalPop = cell(run, 1);

    %% PARALLEL SIMULATIONS
    parfor r = 1:run
        M_pop = m0 + sigma * randn(N0, 1);
        local_mean = NaN(1, T_max);
        local_var  = NaN(1, T_max);
        local_hist = NaN(1, T_max);
        local_beta = NaN(1, T_max);

        for t = 1:T_max
            if isempty(M_pop), continue; end
            N = length(M_pop);
            mean_val = mean(M_pop);

            % Environmental shift: old ? new proofreading constant
            if t <= 500
                mean_error = error_calculator(mean_val, k_old);
            else
                mean_error = error_calculator(mean_val, k_new);
            end
            dm = gen_len * mean_error;
            % Store stats before reproduction
            local_mean(t) = mean_val;
            local_var(t)  = var(M_pop);
            local_hist(t) = N;

            % --- Compute ? gradient ---
            eps_fd = 1e-8;
            f_plus  = error_calculator(mean_val + eps_fd, k_new);
            f_minus = error_calculator(mean_val - eps_fd, k_new);
            df_dm = (f_plus - f_minus) / (2 * eps_fd);

            fbar = mean_error;
            p_kill_bar = fbar / f_th;
            S_bar = 1 - p_kill_bar;
            p_off_bar = exp(5*S_bar) / (1 + exp(5*S_bar));

            % effective fitness
            dlnW_df = -(1/(f_th*(1 - fbar/f_th))) - (5/f_th) * (1 - p_off_bar);
            local_beta(t) = dlnW_df * df_dm;

            % --- Birth–death and mutation ---
            death_mask = false(N,1);
            offspring_list = [];

            for i = 1:N
                m_i = M_pop(i);

                if t <= 500
                    f_i = error_calculator(m_i, k_old);
                else
                    f_i = error_calculator(m_i, k_new);
                end

%                 dm = gen_len * f_i;
                p_kill = f_i / f_th;

                if p_kill >= rand || f_i < 0
                    death_mask(i) = true;
                    continue;
                end

                % Three offspring possibilities
                m1 = m_i - dm;
                m2 = m_i;
                m3 = m_i + dm;

                if any([m1, m2, m3] <= 0), continue; end

                % Reproduction probability
                a = 5;
                s = 1 - p_kill;
                p_off = exp(a*s) / (1 + exp(a*s));

                % Determine number of offspring
                if (p_off - 0.25) > rand
                    n_off = 2;
                else
                    n_off = 1;
                end

                % Random offspring selection
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
                death_mask(i) = true; % parent dies after reproduction
            end

            % Update population
            M_pop = M_pop(~death_mask);
            if ~isempty(offspring_list)
                M_pop = [M_pop; offspring_list];
            end

            % Cap population at N_cap
            if length(M_pop) > N_cap
                M_pop = M_pop(randi(length(M_pop), N_cap, 1));
            end
        end

        % Store replicate data
        mean_m(r,:) = local_mean;
        var_m(r,:)  = local_var;
        N_hist(r,:) = local_hist;
        beta_grad(r,:) = local_beta;
        FinalPop{r} = M_pop;
    end

    %% Aggregate results
    trait_m = mean(mean_m, 1, 'omitnan');
    delta_trait = [NaN, diff(trait_m)];
    var_trait = mean(var_m, 1, 'omitnan');
    N_hist_mean = mean(N_hist, 1, 'omitnan');
    beta_mean = mean(beta_grad, 1, 'omitnan');
    predicted_rate = var_trait .* beta_mean;

    results(n_idx).N0 = N0;
    results(n_idx).trait_m = trait_m;
    results(n_idx).delta_trait = delta_trait;
    results(n_idx).N_hist_mean = N_hist_mean;
    results(n_idx).predicted_rate = predicted_rate;
end

%% === Plot comparative results ===
colors = lines(length(N_list));
x = 1:T_max;

figure('Color','w');
subplot(3,1,1); hold on;
for n_idx = 1:length(N_list)
    plot(x, results(n_idx).N_hist_mean, 'Color', colors(n_idx,:), ...
        'DisplayName', sprintf('N_0 = %d', N_list(n_idx)));
end
xlabel('Generation'); ylabel('Population Size');
title('Population Growth'); legend('Location','best'); grid on;

subplot(3,1,2); hold on;
for n_idx = 1:length(N_list)
    plot(x, results(n_idx).trait_m, 'Color', colors(n_idx,:), ...
        'DisplayName', sprintf('N_0 = %d', N_list(n_idx)));
end
xlabel('Generation'); ylabel('Mean Trait');
title('Evolution of Trait'); legend('Location','best'); grid on;

subplot(3,1,3); hold on;
for n_idx = 1:length(N_list)
    plot(x, results(n_idx).delta_trait, 'Color', colors(n_idx,:), ...
        'DisplayName', sprintf('N_0 = %d', N_list(n_idx)));
end
xlabel('Generation'); ylabel('\Delta m');
title('Adaptation Rate (Observed)'); legend('Location','best'); grid on;
