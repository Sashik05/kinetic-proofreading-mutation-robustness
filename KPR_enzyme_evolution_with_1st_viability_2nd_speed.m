clear; clc;

%% Start parallel pool
if isempty(gcp('nocreate')), parpool; end

%% Parameters
N0 = 1000;
k_old = 1;
k_new = 42.0098;
f_th = 1e-6;
sigma = 1e-3;
T_max = 2500;
run = 16;
N_cap = 2 * N0;
gen_len = 1e06;

%% Trait bounds and optima
m_values = double(vpa(sort(driving_function(f_th, k_old)), 3));
[~, m0, ~, ~] = err_fraction(k_old);
[~, m0_new, ~, ~] = err_fraction(k_new);
m_opt = m0_new;

%% Preallocate global storage
mean_m = NaN(run, T_max);
var_m  = NaN(run, T_max);
N_hist = NaN(run, T_max);
beta_grad = NaN(run, T_max);
FinalPop = cell(run, 1);

%% PARALLEL REPLICATE SIMULATIONS
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

        if t <= 500
            mean_error = error_calculator(mean_val, k_old);
        else
            mean_error = error_calculator(mean_val, k_new);
        end

        %% store stats
        local_mean(t) = mean_val;
        local_var(t)  = var(M_pop);
        local_hist(t) = N;

        %% theoretical beta gradient
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

        %% birth death mutation
        death_mask = false(N,1);
        offspring_list = [];

        for i = 1:N
            m_i = M_pop(i);

            if t <= 500
                f_i = error_calculator(m_i, k_old);
            else
                f_i = error_calculator(m_i, k_new);
            end

            dm = gen_len * f_i;

            p_kill = f_i / f_th;
            if p_kill >= rand || f_i < 0
                death_mask(i) = true;
                continue;
            end

            m1 = m_i - dm;
            m2 = m_i;
            m3 = m_i + dm;

            if any([m1, m2, m3] <= 0), continue; end

            a = 5;
            s = 1 - p_kill;
            p_off = exp(a*s) / (1 + exp(a*s));

            if (p_off - 0.25) > rand
                n_off = 2;
            else
                n_off = 1;
            end

            r_pick = rand;
            chosen_off = [];

            if n_off == 1
                if r_pick < 1/3
                    chosen_off = m1;
                elseif r_pick < 2/3
                    chosen_off = m2;
                else
                    chosen_off = m3;
                end
            else
                if r_pick < 1/3
                    chosen_off = [m1; m3];
                elseif r_pick < 2/3
                    chosen_off = [m1; m2];
                else
                    chosen_off = [m2; m3];
                end
            end

            offspring_list = [offspring_list; chosen_off];
            death_mask(i) = true;
        end

        %% apply births & deaths
        M_pop = M_pop(~death_mask);
        if ~isempty(offspring_list)
            M_pop = [M_pop; offspring_list];
        end

        %% ============================
        %% CAP POPULATION (MODIFIED)
        %% viability FIRST, speed SECOND
        %% ============================

        if length(M_pop) > N_cap

            if t <= 500
                kx_current = k_old;
            else
                kx_current = k_new;
            end

            % compute survival
            f_vals = arrayfun(@(m) error_calculator(m, kx_current), M_pop);
            s_vals = max(0, 1 - f_vals./f_th);

            % compute speed
            speeds = correct_product_formation_rate(M_pop, kx_current);

            % lexicographic ranking
            ranking_matrix = [s_vals(:), speeds(:)];

            % highest viability first, break ties with speed
            [~, order] = sortrows(ranking_matrix, [-1 -2]);

            keep_idx = order(1:N_cap);
            M_pop = M_pop(keep_idx);
        end

    end

    mean_m(r,:) = local_mean;
    var_m(r,:)  = local_var;
    N_hist(r,:) = local_hist;
    beta_grad(r,:) = local_beta;
    FinalPop{r} = M_pop;
end

%% statistics
trait_m = mean(mean_m, 1, 'omitnan');
stderr_m = std(mean_m, 0, 1, 'omitnan') / sqrt(run);
z = 3.291;

ci_upper = trait_m + z * stderr_m;
ci_lower = trait_m - z * stderr_m;

var_trait = mean(var_m, 1, 'omitnan');
N_hist_mean = mean(N_hist, 1, 'omitnan');

delta_trait = [NaN, diff(trait_m)];
norm_rate = delta_trait;

beta_mean = mean(beta_grad, 1, 'omitnan');
predicted_rate = var_trait .* beta_mean;

%% plots
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineLineWidth',1.5);
x = 1:T_max;

figure('Color','w');
subplot(4,1,1);
plot(N_hist_mean,'k');
xlabel('Generation'); ylabel('Pop Size');

subplot(4,1,2); hold on;
fill([x fliplr(x)],[ci_upper fliplr(ci_lower)],[0.8 0.8 1],'EdgeColor','none','FaceAlpha',0.4);
plot(trait_m,'b'); yline(m_opt,'--r','Optimum m');
xlabel('Generation'); ylabel('Mean m');

subplot(4,1,3);
plot(delta_trait,'g');
xlabel('Generation'); ylabel('\Delta m');

subplot(4,1,4);
plot(norm_rate,'m'); hold on;
plot(predicted_rate,'--k');
xlabel('Generation'); ylabel('Rate');
legend('Observed','Predicted');

figure('Color','w');
histogram(FinalPop{1},100,'FaceColor',[0.2 0.2 0.8]);
xlabel('Trait m'); ylabel('Count');
