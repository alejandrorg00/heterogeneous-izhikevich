% inhibition_efficiency.m
% Experimento: eficiencia energética vs proporción inhibitoria p_I
% Modelo: Izhikevich RS (E) y FS (I), sin plasticidad.
% Salidas: 3 figuras + CSV con resultados.

clear; clc; close all;

%% --- Parámetros del experimento ---
pI_list      = 0.00:0.05:0.50;   % barrido de proporción inhibitoria
N_total      = 1000;             % neuronas totales
D_max        = 5;                % retardo máximo (no se usa explícitamente aquí)
T_sec        = 20;                % duración (s)
dt_ms        = 1;                % paso (ms)
seed0        = 2024;             % semilla base
energy_mode  = 'driving';         % 'simple' o 'driving'
save_csv     = true;
out_csv      = sprintf('inhibition_efficiency_%s_%s.csv', energy_mode, datestr(now,'yyyymmdd_HHMMss'));

% Drive externo (tálamo/ruido)
drive_amp    = 13.0;             % amplitud de ruido uniforme ~U(-amp/2, amp/2)
tonic_E      = 1.0;              % DC en E (para evitar apagado con pI alto)
tonic_I      = 0.0;              % DC en I

% Conectividad
p_conn       = 0.10;             % probabilidad de conexión
wE_mu        = 0.6;  wE_sd = 0.2;  wE_min = 0.05; wE_max = 1.5;
wI_mu        = 1.0;  wI_sd = 0.2;  wI_min = 0.05; wI_max = 2.0;

% Sinapsis (traza exponencial)
tau_syn_E    = 8.0;              % ms
tau_syn_I    = 5.0;              % ms

% Reversales (para energy_mode='driving')
Eexc_mV      = 0;
Einh_mV      = -70;

%% --- Plasticidad ---
learning_rule = 'stdp';  % 'none' | 'stdp' | 'da-stdp'

% Parámetros STDP (aditivo)
stdp.tau_plus   = 20.0;    % ms, traza presináptica
stdp.tau_minus  = 20.0;    % ms, traza postsináptica
stdp.A_plus     = 0.01;    % magnitud LTP
stdp.A_minus    = 0.012;   % magnitud LTD
stdp.eta        = 0.1;     % factor global de aprendizaje
stdp.clip_min   = 0.05;    % límites para WE (ya los tienes arriba)
stdp.clip_max   = 1.50;

% (Solo para DA-STDP; se ignora si learning_rule='stdp' o 'none')
stdp.DA_decay   = 0.99;    % decaimiento por paso (≈10 ms si dt=1 ms y 10 steps?)
stdp.DA_pulse   = 0.5;     % liberación cuando hay “reward”
stdp.pair_interval_ms = 20; % ventana de coincidencia pre→post que dispara reward
stdp.reward_delay_ms  = 200; % retardo del “reward” tras coincidencia

% Elige un par sináptico E→E al azar para generar “rewards” (estilo Izhikevich)
stdp.use_reward_pair = true;  % si false, no habrá eventos de reward


%% --- Resultados ---
Sp = zeros(size(pI_list));   % spikes excitatorios
En = zeros(size(pI_list));   % energía acumulada
Ef = zeros(size(pI_list));   % eficiencia = Sp/En
Ne_list = zeros(size(pI_list));
Ni_list = zeros(size(pI_list));

fprintf('--- Running sweep (energy_mode=%s) ---\n', energy_mode);
for k = 1:numel(pI_list)
    pI  = pI_list(k);
    seed = seed0 + k;
    out = run_izh_efficiency_once( ...
        pI, N_total, T_sec, dt_ms, seed, ...
        p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
        wI_mu, wI_sd, wI_min, wI_max, ...
        tau_syn_E, tau_syn_I, ...
        drive_amp, tonic_E, tonic_I, ...
        energy_mode, Eexc_mV, Einh_mV, ...
        learning_rule, stdp);


    Sp(k) = out.spikesE;
    En(k) = out.energyE;
    Ef(k) = out.efficiency;
    Ne_list(k) = out.Ne;
    Ni_list(k) = out.Ni;

    fprintf('pI=%.2f  Ne=%4d Ni=%4d  spikesE=%7d  energyE=%10.1f  eff=%.4g\n', ...
        pI, out.Ne, out.Ni, out.spikesE, out.energyE, out.efficiency);
end

%% --- Plots ---
figure('Color','w'); 
plot(pI_list, Sp, 'o-','LineWidth',1.2); grid on;
xlabel('Inhibition proportion p_I'); ylabel('Total excitatory spikes');
title('Excitatory activity vs inhibition');

figure('Color','w'); 
plot(pI_list, En, 's-','LineWidth',1.2); grid on;
xlabel('Inhibition proportion p_I'); 
if strcmpi(energy_mode,'simple')
    ylabel('Energy proxy \Sigma_E |I| \cdot dt');
else
    ylabel('Energy proxy \Sigma_E (|I_{exc}||V-E_{exc}| + |I_{inh}||V-E_{inh}|) \cdot dt');
end
title('Energy consumption vs inhibition');

figure('Color','w'); 
plot(pI_list, Ef, '^-','LineWidth',1.2); grid on;
xlabel('Inhibition proportion p_I'); ylabel('Efficiency (spikes / energy)');
title('Energy efficiency vs inhibition');

%% --- Guardar CSV ---
% if save_csv
%     T = table(pI_list(:), Ne_list(:), Ni_list(:), Sp(:), En(:), Ef(:), ...
%         'VariableNames', {'pI','Ne','Ni','spikesE','energyE','efficiency'});
%     writetable(T, out_csv);
%     fprintf('Resultados guardados en: %s\n', out_csv);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           FUNCIONES LOCALES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = run_izh_efficiency_once( ...
    pI, N_total, T_sec, dt_ms, seed, ...
    p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
    wI_mu, wI_sd, wI_min, wI_max, ...
    tau_syn_E, tau_syn_I, ...
    drive_amp, tonic_E, tonic_I, ...
    energy_mode, Eexc_mV, Einh_mV, ...
    learning_rule, stdp)
% Simula red E/I tipo Izhikevich y calcula métricas energéticas.
% Añade plasticidad STDP (y opcionalmente DA-STDP) en E→*.

    rng(seed);

    % --- tamaños ---
    Ni = round(pI * N_total);
    Ne = N_total - Ni;
    N  = N_total;

    % --- parámetros Izhikevich ---
    a = [0.02*ones(Ne,1); 0.1*ones(Ni,1)];  % RS / FS
    b = [0.20*ones(Ne,1); 0.20*ones(Ni,1)];
    c = [-65*ones(Ne,1); -65*ones(Ni,1)];
    d = [ 8.0*ones(Ne,1);  2.0*ones(Ni,1)];

    v = -65 * ones(N,1);
    u = b .* v;

    % --- Conectividad aleatoria: WE (N×Ne) y WI (N×Ni) ---
    % E -> *
    mask_E = rand(N, Ne) < p_conn;
    [rowE, colE] = find(mask_E);
    wE = wE_mu + wE_sd .* randn(numel(rowE), 1);
    wE = max(wE_min, min(wE_max, wE));
    WE = sparse(rowE, colE, wE, N, Ne);

    % I -> *
    mask_I = rand(N, Ni) < p_conn;
    [rowI, colI] = find(mask_I);
    wI = -(wI_mu + wI_sd .* randn(numel(rowI), 1));
    wI = -max(wI_min, min(wI_max, -wI));
    WI = sparse(rowI, colI, wI, N, Ni);

    % Matriz completa postsyn×presyn (solo para la corriente)
    W = [WE, WI];

    % Quitar autoconexiones
    W = W - spdiags(diag(W), 0, N, N);

    % --- Sinapsis: traza exponencial por presináptica (corrientes) ---
    alpha_E = exp(-dt_ms / tau_syn_E);
    alpha_I = exp(-dt_ms / tau_syn_I);
    s = zeros(N,1); % traza por neurona presináptica
    alpha = [alpha_E*ones(Ne,1); alpha_I*ones(Ni,1)];

    % --- Tiempo ---
    steps = round(1000*T_sec / dt_ms);

    % --- Métricas ---
    energyE = 0;
    spikesE_total = 0;

    %% STDP: Estados
    do_plastic = ~strcmpi(learning_rule,'none');

    if do_plastic
        % Trazas STDP
        tau_plus  = stdp.tau_plus;
        tau_minus = stdp.tau_minus;
        A_plus    = stdp.A_plus;
        A_minus   = stdp.A_minus;
        eta       = stdp.eta;
        clip_min  = stdp.clip_min;
        clip_max  = stdp.clip_max;

        decay_plus  = exp(-dt_ms / tau_plus);
        decay_minus = exp(-dt_ms / tau_minus);

        % Trazas presin (solo para E) y postsin (para todas las neuronas)
        x_preE  = zeros(Ne,1);   % traza presináptica para excitadoras
        x_post  = zeros(N,1);    % traza postsináptica para todas

        % Dopamina y “rewards” (solo si DA-STDP)
        DA = 0;
        if strcmpi(learning_rule,'da-stdp')
            DA_decay = stdp.DA_decay;
            DA_pulse = stdp.DA_pulse;
            pair_window = stdp.pair_interval_ms;
            reward_delay = stdp.reward_delay_ms;

            % Selecciona un par E→E aleatorio si procede
            if stdp.use_reward_pair && Ne >= 2
                pre_pair  = randi(Ne);
                % elegimos un postsin excitatorio distinto del pre
                post_pair = randi(Ne);
                if post_pair == pre_pair, post_pair = mod(post_pair,Ne)+1; end
                % Programador simple de rewards
                pending_rewards = false(steps,1);
                last_pre_spike  = -Inf;
            else
                pre_pair = [];
                post_pair = [];
                pending_rewards = false(steps,1);
                last_pre_spike  = -Inf;
            end
        end
    end

    for t = 1:steps
        % Drive externo + DC
        I_bg = drive_amp * (rand(N,1) - 0.5);
        I_bg(1:Ne)     = I_bg(1:Ne)     + tonic_E;
        I_bg(Ne+1:end) = I_bg(Ne+1:end) + tonic_I;

        % Corriente sináptica
        I_syn = W * s;

        % Total input
        I = I_bg + I_syn;

        % Energía
        if strcmpi(energy_mode, 'simple')
            IE = I(1:Ne);
            energyE = energyE + sum(abs(IE)) * dt_ms;
        else
            I_postE = I_syn(1:Ne);
            Iexc = max(I_postE, 0);
            Iinh = min(I_postE, 0);
            V_E  = v(1:Ne);
            P_step = sum( abs(Iexc).*abs(V_E - Eexc_mV) + ...
                          abs(Iinh).*abs(V_E - Einh_mV) );
            energyE = energyE + P_step * dt_ms;
        end

        % Izhikevich update (dos medios pasos para estabilidad)
        v = v + 0.5*((0.04*v + 5).*v + 140 - u + I);
        v = v + 0.5*((0.04*v + 5).*v + 140 - u + I);
        u = u + a .* (0.2*v - u);

        % Spikes
        fired = find(v >= 30);
        if ~isempty(fired)
            v(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
            spikesE_total = spikesE_total + sum(fired <= Ne);
        end

        % Trazas sinápticas (corriente)
        s = alpha .* s;
        if ~isempty(fired)
            s(fired) = s(fired) + 1.0;
        end

        %% STDP: actualización de trazas y pesos
        if do_plastic
            % Decaimiento de trazas STDP
            x_preE = x_preE * decay_plus;
            x_post = x_post * decay_minus;

            if ~isempty(fired)
                % Actualiza trazas con spikes de este paso
                firedE = fired(fired <= Ne);  % indices presin excitatorios
                x_post(fired) = x_post(fired) + 1;     % posts para todos
                if ~isempty(firedE)
                    x_preE(firedE) = x_preE(firedE) + 1; % pres E
                end

                % LTP en spikes presin E: Δw_ij += eta*A_plus * x_post(i) * M
                for jj = firedE.'   % jj: neurona presin E
                    % Sinapsis de salida desde jj: WE(:,jj)
                    [rows_i, ~, vals] = find(WE(:,jj)); %#ok<ASGLU>
                    if ~isempty(rows_i)
                        modFactor = 1.0;
                        if strcmpi(learning_rule,'da-stdp')
                            modFactor = (1.0 + DA);
                        end
                        dWij = eta * A_plus * x_post(rows_i) * modFactor;
                        % aplica y clip
                        new_w = full(WE(rows_i,jj)) + dWij;
                        new_w = min(max(new_w, clip_min), clip_max);
                        WE(rows_i,jj) = new_w;
                    end
                end

                % LTD en spikes postsin: Δw_ij -= eta*A_minus * x_preE(j) * M
                for ii = fired.'     % ii: neurona postsin (todas)
                    % Sinapsis excitatorias que llegan a ii: WE(ii,:)
                    [~, cols_j, ~] = find(WE(ii,:));
                    if ~isempty(cols_j)
                        modFactor = 1.0;
                        if strcmpi(learning_rule,'da-stdp')
                            modFactor = (1.0 + DA);
                        end
                        dWij = - eta * A_minus * x_preE(cols_j).' * modFactor;
                        new_w = full(WE(ii,cols_j)) + dWij;
                        new_w = min(max(new_w, clip_min), clip_max);
                        WE(ii,cols_j) = new_w;
                    end
                end
            end

            % Reconstruye W con WE actualizado (WI fijo)
            W = [WE, WI];
            W = W - spdiags(diag(W), 0, N, N);

            % DA-STDP: lógica de “reward” sencilla estilo Izhikevich
            if strcmpi(learning_rule,'da-stdp')
                % Programa reward si hay coincidencia n_pre -> n_post en ventana
                if stdp.use_reward_pair && ~isempty(fired)
                    t_ms = t * dt_ms;
                    if any(fired == pre_pair)
                        last_pre_spike = t_ms;
                    end
                    if any(fired == post_pair)
                        if (t_ms - last_pre_spike) > 0 && (t_ms - last_pre_spike) <= stdp.pair_interval_ms
                            idx_reward = round((t_ms + stdp.reward_delay_ms) / dt_ms);
                            if idx_reward >= 1 && idx_reward <= steps
                                pending_rewards(idx_reward) = true;
                            end
                        end
                    end
                end

                % Aplica rewards pendientes
                if pending_rewards(t)
                    DA = DA + DA_pulse;
                end
                % Decaimiento DA
                DA = DA * DA_decay;
            end
        end
    end

    out.Ne = Ne; out.Ni = Ni;
    out.spikesE = spikesE_total;
    out.energyE = energyE;
    out.efficiency = spikesE_total / max(energyE, eps);

    % (Opcional) devuelve pesos finales y, si procede, el par de reward
    out.WE_final = WE;
    if do_plastic && strcmpi(learning_rule,'da-stdp')
        out.reward_pair = struct('pre',pre_pair,'post',post_pair);
    end
end
