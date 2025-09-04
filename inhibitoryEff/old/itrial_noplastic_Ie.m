% inhibition_efficiency_drive_sweep.m
% Compara Eficiencia vs pI para distintos niveles de drive_amp
% Basado en el experimento sin plasticidad (Izhikevich RS/FS, sinapsis exponencial)

clear; clc; close all;

%% --- Parámetros del experimento ---
pI_list       = 0.10:0.01:0.50;   % barrido de proporción inhibitoria
drive_amp_list= 10:10:100;     % niveles de drive a comparar
N_total       = 1000;             % neuronas totales
T_sec         = 2;                % duración (s)
dt_ms         = 1;                % paso (ms)
seed0         = 2024;             % semilla base
energy_mode   = 'simple';         % 'simple' o 'driving'
save_csv      = true;
out_csv       = sprintf('drive_sweep_eff_%s_%s.csv', energy_mode, datestr(now,'yyyymmdd_HHMMss'));

% drive tónico (corriente DC)
tonic_E       = 1.0;
tonic_I       = 0.0;

% Conectividad
p_conn        = 0.10;
wE_mu=0.6; wE_sd=0.2; wE_min=0.05; wE_max=1.5;
wI_mu=1.0; wI_sd=0.2; wI_min=0.05; wI_max=2.0;

% Sinapsis (traza exponencial)
tau_syn_E     = 8.0;              % ms
tau_syn_I     = 5.0;              % ms

% Reversales (para energy_mode='driving')
Eexc_mV       = 0;
Einh_mV       = -70;

%% --- Buffers de resultados ---
nP  = numel(pI_list);
nD  = numel(drive_amp_list);
Sp  = zeros(nP, nD);   % spikes excitatorios
En  = zeros(nP, nD);   % energía
Ef  = zeros(nP, nD);   % eficiencia
NeM = zeros(nP, nD);   % tamaño E
NiM = zeros(nP, nD);   % tamaño I

fprintf('=== Sweeping pI for multiple drive_amp (energy_mode=%s) ===\n', energy_mode);

for j = 1:nD
    drive_amp = drive_amp_list(j);
    fprintf('--- drive_amp = %.1f ---\n', drive_amp);
    for k = 1:nP
        pI   = pI_list(k);
        seed = seed0 + j*100 + k;    % distinta semilla por condición
        out  = run_izh_efficiency_once( ...
            pI, N_total, T_sec, dt_ms, seed, ...
            p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
            wI_mu, wI_sd, wI_min, wI_max, ...
            tau_syn_E, tau_syn_I, ...
            drive_amp, tonic_E, tonic_I, ...
            energy_mode, Eexc_mV, Einh_mV);

        Sp(k,j)  = out.spikesE;
        En(k,j)  = out.energyE;
        Ef(k,j)  = out.efficiency;
        NeM(k,j) = out.Ne;
        NiM(k,j) = out.Ni;

        fprintf('pI=%.2f | Ne=%4d Ni=%4d | spikesE=%7d | energy=%.1f | eff=%.4g\n', ...
            pI, out.Ne, out.Ni, out.spikesE, out.energyE, out.efficiency);
    end
end

%% --- Plots comparativos ---
% 1) Eficiencia vs pI, múltiples curvas (una por drive_amp)
figure('Color','w'); hold on; grid on;
marks = {'o-','s-','^-','d-','x-','+-'};  % estilos diferentes
for j = 1:nD
    mk = marks{1+mod(j-1,numel(marks))};
    plot(pI_list, Ef(:,j), mk, 'LineWidth', 1.2, 'DisplayName', sprintf('drive\\_amp = %.0f', drive_amp_list(j)));
end
xlabel('Inhibition proportion p_I');
ylabel('Efficiency (spikes / energy)');
title(sprintf('Energy efficiency vs p_I (mode = %s)', energy_mode));
legend('Location','best');

% 2) Spikes E vs pI (diagnóstico)
figure('Color','w'); hold on; grid on;
for j = 1:nD
    mk = marks{1+mod(j-1,numel(marks))};
    plot(pI_list, Sp(:,j), mk, 'LineWidth', 1.2, 'DisplayName', sprintf('drive\\_amp = %.0f', drive_amp_list(j)));
end
xlabel('Inhibition proportion p_I');
ylabel('Total excitatory spikes');
title('Excitatory activity vs p_I (multi-drive)');
legend('Location','best');

% 3) Energía vs pI (diagnóstico)
figure('Color','w'); hold on; grid on;
for j = 1:nD
    mk = marks{1+mod(j-1,numel(marks))};
    plot(pI_list, En(:,j), mk, 'LineWidth', 1.2, 'DisplayName', sprintf('drive\\_amp = %.0f', drive_amp_list(j)));
end
xlabel('Inhibition proportion p_I');
if strcmpi(energy_mode,'simple')
    ylabel('Energy proxy \Sigma_E |I| \cdot dt');
else
    ylabel('Energy proxy \Sigma_E (|I_{exc}||V-E_{exc}| + |I_{inh}||V-E_{inh}|) \cdot dt');
end
title('Energy consumption vs p_I (multi-drive)');
legend('Location','best');

%% --- Guardar CSV con todas las condiciones ---
if save_csv
    % Unir resultados en tabla larga: columnas (drive_amp, pI, Ne, Ni, SpikesE, EnergyE, Efficiency)
    rows = {};
    for j = 1:nD
        for k = 1:nP
            rows(end+1,:) = {drive_amp_list(j), pI_list(k), NeM(k,j), NiM(k,j), Sp(k,j), En(k,j), Ef(k,j)}; %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', {'drive_amp','pI','Ne','Ni','spikesE','energyE','efficiency'});
    writetable(T, out_csv);
    fprintf('Tabla combinada guardada en: %s\n', out_csv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           FUNCIÓN LOCAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = run_izh_efficiency_once( ...
    pI, N_total, T_sec, dt_ms, seed, ...
    p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
    wI_mu, wI_sd, wI_min, wI_max, ...
    tau_syn_E, tau_syn_I, ...
    drive_amp, tonic_E, tonic_I, ...
    energy_mode, Eexc_mV, Einh_mV)
% Simula red E/I tipo Izhikevich y calcula métricas (sin plasticidad)

    rng(seed);

    % --- tamaños ---
    Ni = round(pI * N_total);
    Ne = N_total - Ni;
    N  = N_total;

    % --- Izhikevich RS/FS ---
    a = [0.02*ones(Ne,1); 0.1*ones(Ni,1)];
    b = [0.20*ones(Ne,1); 0.20*ones(Ni,1)];
    c = [-65*ones(Ne,1); -65*ones(Ni,1)];
    d = [ 8.0*ones(Ne,1);  2.0*ones(Ni,1)];
    v = -65 * ones(N,1);
    u = b .* v;

    % --- Conectividad dispersa ---
    % WE: N×Ne (positiva)
    mask_E = rand(N, Ne) < p_conn;
    [rowE, colE] = find(mask_E);
    wE = wE_mu + wE_sd .* randn(numel(rowE), 1);
    wE = max(wE_min, min(wE_max, wE));
    WE = sparse(rowE, colE, wE, N, Ne);

    % WI: N×Ni (negativa)
    mask_I = rand(N, Ni) < p_conn;
    [rowI, colI] = find(mask_I);
    wI = -(wI_mu + wI_sd .* randn(numel(rowI), 1));
    wI = -max(wI_min, min(wI_max, -wI));
    WI = sparse(rowI, colI, wI, N, Ni);

    % W completa (postsyn × presyn)
    W = [WE, WI];
    % Quitar autoconexiones
    W = W - spdiags(diag(W), 0, N, N);

    % --- Sinapsis: traza exponencial ---
    alpha_E = exp(-dt_ms / tau_syn_E);
    alpha_I = exp(-dt_ms / tau_syn_I);
    s = zeros(N,1);
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

        % Total
        I = I_bg + I_syn;

        % Energía en E
        if strcmpi(energy_mode,'simple')
            IE = I(1:Ne);
            energyE = energyE + sum(abs(IE)) * dt_ms;
        else
            % Aproximación driving force usando signo de I_syn en E
            I_postE = I_syn(1:Ne);
            Iexc = max(I_postE, 0);
            Iinh = min(I_postE, 0);
            V_E  = v(1:Ne);
            P = sum( abs(Iexc).*abs(V_E - Eexc_mV) + abs(Iinh).*abs(V_E - Einh_mV) );
            energyE = energyE + P * dt_ms;
        end

        % Izhikevich (dos semipasos)
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
            s(fired) = s(fired) + 1.0; % % (MATLAB permite s(fired)=..., Octave usa s(fired)=...)
        end
    end

    out.Ne = Ne; out.Ni = Ni;
    out.spikesE = spikesE_total;
    out.energyE = energyE;
    out.efficiency = spikesE_total / max(energyE, eps);
end
