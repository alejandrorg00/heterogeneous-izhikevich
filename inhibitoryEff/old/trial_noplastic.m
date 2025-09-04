% inhibition_efficiency.m
% Experimento: eficiencia energética vs proporción inhibitoria p_I
% Modelo: Izhikevich RS (E) y FS (I), sin plasticidad.
% Salidas: 3 figuras + CSV con resultados.

clear; clc; close all;

%% --- Parámetros del experimento ---
pI_list      = 0.00:0.05:0.50;   % barrido de proporción inhibitoria
N_total      = 1000;             % neuronas totales
D_max        = 5;                % retardo máximo (no se usa explícitamente aquí)
T_sec        = 2;                % duración (s)
dt_ms        = 1;                % paso (ms)
seed0        = 2024;             % semilla base
energy_mode  = 'simple';         % 'simple' o 'driving'
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
        energy_mode, Eexc_mV, Einh_mV);

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
if save_csv
    T = table(pI_list(:), Ne_list(:), Ni_list(:), Sp(:), En(:), Ef(:), ...
        'VariableNames', {'pI','Ne','Ni','spikesE','energyE','efficiency'});
    writetable(T, out_csv);
    fprintf('Resultados guardados en: %s\n', out_csv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           FUNCIONES LOCALES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = run_izh_efficiency_once( ...
    pI, N_total, T_sec, dt_ms, seed, ...
    p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
    wI_mu, wI_sd, wI_min, wI_max, ...
    tau_syn_E, tau_syn_I, ...
    drive_amp, tonic_E, tonic_I, ...
    energy_mode, Eexc_mV, Einh_mV)
% Simula red E/I tipo Izhikevich y calcula métricas energéticas (sin plasticidad)

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

    % --- Conectividad aleatoria: construir WE (N×Ne) y WI (N×Ni) con sparse ---
    % E -> *
    mask_E = rand(N, Ne) < p_conn;                        % logical N×Ne
    [rowE, colE] = find(mask_E);                          % índices (fila, col)
    wE = wE_mu + wE_sd .* randn(numel(rowE), 1);
    wE = max(wE_min, min(wE_max, wE));                    % recorte
    WE = sparse(rowE, colE, wE, N, Ne);                   % N×Ne

    % I -> *
    mask_I = rand(N, Ni) < p_conn;                        % logical N×Ni
    [rowI, colI] = find(mask_I);
    wI = -(wI_mu + wI_sd .* randn(numel(rowI), 1));       % NEGATIVOS
    wI = -max(wI_min, min(wI_max, -wI));                  % recorte conservando signo -
    WI = sparse(rowI, colI, wI, N, Ni);                   % N×Ni

    % Matriz completa postsyn×presyn
    W = [WE, WI];                                         % N×N

    % Quitar autoconexiones
    W = W - spdiags(diag(W), 0, N, N);

    % --- Sinapsis: traza exponencial por neurona presináptica ---
    alpha_E = exp(-dt_ms / tau_syn_E);
    alpha_I = exp(-dt_ms / tau_syn_I);
    s = zeros(N,1);                                       % traza por presin
    alpha = [alpha_E*ones(Ne,1); alpha_I*ones(Ni,1)];

    % --- Tiempo ---
    steps = round(1000*T_sec / dt_ms);

    % --- Métricas ---
    energyE = 0;
    spikesE_total = 0;

    for t = 1:steps
        % Drive externo + DC
        I_bg = drive_amp * (rand(N,1) - 0.5);
        I_bg(1:Ne)       = I_bg(1:Ne)       + tonic_E;
        I_bg(Ne+1:end)   = I_bg(Ne+1:end)   + tonic_I;

        % Corriente sináptica
        I_syn = W * s;

        % Total input
        I = I_bg + I_syn;

        % Energía
        if strcmpi(energy_mode, 'simple')
            IE = I(1:Ne);
            energyE = energyE + sum(abs(IE)) * dt_ms;
        else
            % Descomponer por signo en población E (aprox. driving force)
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

        % Trazas sinápticas
        s = alpha .* s;
        if ~isempty(fired)
            s(fired) = s(fired) + 1.0;
        end
    end

    out.Ne = Ne; out.Ni = Ni;
    out.spikesE = spikesE_total;
    out.energyE = energyE;
    out.efficiency = spikesE_total / max(energyE, eps);
end
