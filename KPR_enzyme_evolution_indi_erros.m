clear; clc;

%% Start parallel pool
if isempty(gcp('nocreate')), parpool; end

%% Parameters
N0 = 1000;
k_old = 1;
k_new = 42.0098;
f_th = 1e-6;
sigma = 1e-3;
T_max = 6*1e05;
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
beta_grad = NaN(run, T_max); % theory-based selection gradient
FinalPop = cell(run, 1);

%% PARALLEL REPLICATE SIMULATIONS
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
        if t <= 500
            mean_error = error_calculator(mean_val, k_old);
        else
            mean_error = error_calculator(mean_val, k_new);
        end
        
        % Store stats before reproduction
        local_mean(t) = mean_val;
        local_var(t)  = var(M_pop);
        local_hist(t) = N;
        
        % --- Compute theory-based β gradient ---
        eps_fd = 1e-8;
        f_plus  = error_calculator(mean_val + eps_fd, k_new);
        f_minus = error_calculator(mean_val - eps_fd, k_new);
        df_dm = (f_plus - f_minus) / (2 * eps_fd);
        
        fbar = mean_error;
        p_kill_bar = fbar / f_th;
        S_bar = 1 - p_kill_bar;
        p_off_bar = exp(5*S_bar) / (1 + exp(5*S_bar));
        
        % effective fitness W(m)
        W_bar = S_bar * p_off_bar;
        
        % derivative d(ln W)/df
        dlnW_df = -(1/(f_th*(1 - fbar/f_th))) - (5/f_th) * (1 - p_off_bar);
        
        % β = d(ln W)/df * df/dm
        local_beta(t) = dlnW_df * df_dm;
        
        % --- Birth–death and mutation ---
        death_mask = false(N,1);
        offspring_list = [];  % dynamic storage for all offspring
        
        for i = 1:N
            m_i = M_pop(i);
            
            % Apply k_old until t <= 100, then k_new
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
            
            % Mutation step: three possible offspring
            m1 = m_i - dm;   % new variant (new driving)
            m2 = m_i;        % neutral (parental trait)
            m3 = m_i + dm;   % new variant (new driving)
            
            if any([m1, m2, m3] <= 0), continue; end
            
            % Reproduction probability
            a = 5;
            s = 1 - p_kill;
            p_off = exp(a*s) / (1 + exp(a*s));
            
            % Determine number of offspring
            if (p_off - 0.25) > rand
                n_off = 2;  % two offspring
            else
                n_off = 1;  % one offspring
            end
            
            % Randomly choose offspring combinations
            r_pick = rand;
            if n_off == 1
                if r_pick < 1/3
                    chosen_off = m1;
                elseif r_pick < 2/3 && r_pick > 1/3
                    chosen_off = m2;
                elseif r_pick > 2/3 && r_pick < 1
                    chosen_off = m3;
                end
                offspring_list = [offspring_list; chosen_off];
            else % n_off == 2
                if r_pick < 1/3
                    chosen_off = [m1; m3];
                elseif r_pick < 2/3 && r_pick > 1/3
                    chosen_off = [m1; m2];
                elseif r_pick > 2/3 && r_pick < 1
                    chosen_off = [m2; m3];
                end
                offspring_list = [offspring_list; chosen_off];
            end
            
            % Parent dies after reproduction
            death_mask(i) = true;
        end
        
        % Apply deaths and births
        M_pop = M_pop(~death_mask);
        if ~isempty(offspring_list)
            M_pop = [M_pop; offspring_list];
        end
        
        % --- ✅ CAP POPULATION AFTER ALL UPDATES ---
        if length(M_pop) > N_cap
            M_pop = M_pop(randi(length(M_pop), N_cap, 1));
        end
    end
    
    % Store per replicate
    mean_m(r,:) = local_mean;
    var_m(r,:)  = local_var;
    N_hist(r,:) = local_hist;
    beta_grad(r,:) = local_beta;
    FinalPop{r} = M_pop;
end


%% Final statistics
trait_m = mean(mean_m, 1, 'omitnan');
stderr_m = std(mean_m, 0, 1, 'omitnan') / sqrt(run);
z = 3.291;  % 99.9% CI
ci_upper = trait_m + z * stderr_m;
ci_lower = trait_m - z * stderr_m;
var_trait = mean(var_m, 1, 'omitnan');
N_hist_mean = mean(N_hist, 1, 'omitnan');

% Raw & normalized adaptation rates
delta_trait = [NaN, diff(trait_m)];
% norm_rate = delta_trait ./ var_trait;
norm_rate = delta_trait;
beta_mean = mean(beta_grad, 1, 'omitnan');
predicted_rate = var_trait .* beta_mean; % Fisher-Lande prediction
% predicted_rate = beta_mean;
%% Publication Plots
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultLineLineWidth', 1.5);
x = 1:T_max;

figure('Color', 'w');
subplot(4,1,1);
plot(N_hist_mean, 'k');
xlabel('Generation'); ylabel('Pop. Size'); title('Population Size vs Time');

subplot(4,1,2); hold on;
fill([x fliplr(x)], [ci_upper fliplr(ci_lower)], [0.8 0.8 1], 'EdgeColor','none','FaceAlpha',0.4);
plot(trait_m, 'b'); yline(m_opt, '--r', 'Optimum m');
xlabel('Generation'); ylabel('Mean m'); title('Trait Evolution');

subplot(4,1,3);
plot(delta_trait, 'g');
xlabel('Generation'); ylabel('\Delta m'); title('Raw Adaptation Rate');

subplot(4,1,4);
plot(norm_rate, 'm', 'DisplayName','Observed \beta');
hold on;
plot(predicted_rate, '--k', 'DisplayName','Predicted (Var \times \beta)');
xlabel('Generation'); ylabel('Rate');
title('Observed vs Predicted Rate of Adaptation');
legend('Location','best');

figure('Color', 'w');
histogram(FinalPop{1}, 100, 'FaceColor', [0.2 0.2 0.8]);
xlabel('Trait m'); ylabel('Count'); title('Final Trait Distribution (Sample Replicate)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

%% ---- Subplot 1: Observed beta ----
subplot(1,2,1); % side by side
plot(norm_rate, 'm', 'DisplayName', 'Observed \beta');
xlabel('Generation'); ylabel('Rate');
title('Observed \beta');
grid on;

% Zoomed region for subplot 1
x1 = 60:300;
y1 = norm_rate(x1);

% Inset for subplot 1
ax1 = axes('Position',[0.25 0.55 0.2 0.3]); % adjust position
box on; hold on;
plot(x1, y1, 'm');
xlim([60 300]);
ylim([min(y1) max(y1)]);
title('Zoom','FontSize',8);

%% ---- Subplot 2: Predicted rate ----
subplot(1,2,2);
plot(predicted_rate, '--k', 'DisplayName', 'Predicted (Var \times \beta)');
xlabel('Generation'); ylabel('Rate');
title('Predicted Rate');
grid on;

% Zoomed region for subplot 2
x2 = 1:30;
y2 = predicted_rate(x2);

% Inset for subplot 2
ax2 = axes('Position',[0.7 0.55 0.2 0.3]); % adjust position
box on; hold on;
plot(x2, y2, '--k');
xlim([1 30]);
ylim([min(y2) max(y2)]);
title('Zoom','FontSize',8);
