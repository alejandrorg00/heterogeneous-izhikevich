% Alejandro Rodriguez-Garcia
% 03/09/25
% Network of RS (E) and FS (I) neurons with STDP. Sweep pI and measure
% spikes, energy, efficiency. For pI=0.20, record and plot a raster and
% class-separated synaptic current traces: E->E, I->E, E->I, I->I.

clear; clc; close all;
seed0       = 2024;

%% Simulation parameters
pI_list = 0.00:0.05:0.50;   % inhibitory sweep
N_total = 1000;             % total neurons
T_sec   = 1;                % duration (s)
dt_ms   = 1;                % time step (ms)

% Thalamic drive
drive_amp = 13.0;           % ~U(-amp/2, amp/2)
tonic_E   = 1.0;            % keeps E from being quickly shut down
tonic_I   = 0.0;

% Connectivity initialization
p_conn = 0.10;  % connection probability

% E weight stats (magnitudes)
wE_mu = 0.6;  wE_sd = 0.2;  wE_min = 0.00; wE_max = 2.5;

% I weight stats (magnitudes)
wI_mu = 1.0;  wI_sd = 0.2;  wI_min = 0.00; wI_max = 4.0;

% Synaptic time-constants (for postsynaptic traces)
tau_syn_E = 8.0;  % ms
tau_syn_I = 5.0;  % ms

%% STDP config
learning_rule = 'stdp';   % 'none' | 'stdp'
stdp.tau_plus   = 20.0;   % ms (pre trace)
stdp.tau_minus  = 20.0;   % ms (post trace)
stdp.A_plus     = 0.01;   % LTP scale
stdp.A_minus    = 0.012;  % LTD scale
stdp.eta        = 0.10;   % learning rate

%% Sweep results
Sp      = zeros(size(pI_list));   % total E spikes
En      = zeros(size(pI_list));   % energy proxy
Ef      = zeros(size(pI_list));   % efficiency = Sp/En
Ne_list = zeros(size(pI_list));
Ni_list = zeros(size(pI_list));

fprintf('--- Running sweep ---\n');
for k = 1:numel(pI_list)
    pI   = pI_list(k);
    seed = seed0 + k;

    out = run_izh_efficiency_once( ...
        pI, N_total, T_sec, dt_ms, seed, ...
        p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
        wI_mu, wI_sd, wI_min, wI_max, ...
        tau_syn_E, tau_syn_I, ...
        drive_amp, tonic_E, tonic_I, ...
        learning_rule, stdp, false);  % no recording during sweep

    Sp(k)      = out.spikesE;
    En(k)      = out.energyE;
    Ef(k)      = out.efficiency;
    Ne_list(k) = out.Ne;
    Ni_list(k) = out.Ni;

    fprintf('pI=%.2f  Ne=%4d Ni=%4d  spikesE=%7d  energyE=%10.1f  eff=%.4g\n', ...
        pI, out.Ne, out.Ni, out.spikesE, out.energyE, out.efficiency);
end

%% --- Sweep plots ---
figure('Color','w'); 
plot(pI_list, Sp, 'o-','LineWidth',1.2); grid on;
xlabel('Inhibition proportion p_I'); ylabel('Total excitatory spikes');
title('Excitatory activity vs inhibition');

figure('Color','w'); 
plot(pI_list, En, 's-','LineWidth',1.2); grid on;
xlabel('Inhibition proportion p_I'); 
ylabel('Energy proxy \Sigma_E |I| \cdot dt');
title('Energy consumption vs inhibition');

figure('Color','w'); 
plot(pI_list, Ef, '^-','LineWidth',1.2); grid on;
xlabel('Inhibition proportion p_I'); ylabel('Efficiency (spikes / energy)');
title('Energy efficiency vs inhibition'); 

%% Energy efficiency vs inhibition
fig = figure('Color','w', 'Position', [100, 100, 560*1.25, 420*1.25]); hold on;

% Plot with coral/orangered line and filled markers
plot(pI_list, Ef, 'o-', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [1.0, 0.5, 0.31], ...   % coral fill
    'Color', [1.0, 0.27, 0.0]);                % orangered line

% Labels and title
xlabel('Inhibition proportion p_I', 'FontSize', 14, 'FontWeight','bold');
ylabel('Efficiency (spikes / energy)', 'FontSize', 14, 'FontWeight','bold');

% Axis improvements
grid off; 
box off;
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 13;
ax.GridAlpha = 0.3;   % softer grid
ax.GridLineStyle = '--';

% Adjust limits if needed
xlim([min(pI_list) max(pI_list)]);
ylim([0.0015 max(Ef)*1.1]);

% Export at 600 dpi (PNG raster + EPS/SVG vectorial)
% print(fig, 'efficiency_plot', '-dpng', '-r600');
% print(fig, 'efficiency_plot', '-depsc');
% print(fig, 'efficiency_plot', '-dsvg');



%% --- Detailed recording only for pI = 0.20 ---
pI_target = 0.20;
seed_rec  = 999;  % separate seed for the recorded run (optional)
rec = run_izh_efficiency_once( ...
    pI_target, N_total, T_sec, dt_ms, seed_rec, ...
    p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
    wI_mu, wI_sd, wI_min, wI_max, ...
    tau_syn_E, tau_syn_I, ...
    drive_amp, tonic_E, tonic_I, ...
    learning_rule, stdp, true);    % <-- record

%% --- Raster plot (E blue, I red) ---
figure('Color','w', 'Position', [100, 100, 900, 400]); hold on;
if ~isempty(rec.raster.t)
    isE = rec.raster.i <= rec.Ne;        % E neuron indices
    isI = ~isE;
    plot(rec.raster.t(isE), rec.raster.i(isE), '.r', 'MarkerSize', 6); % default color
    plot(rec.raster.t(isI), rec.raster.i(isI), '.b', 'MarkerSize', 6); % will be a second color
end
xlabel('Time (ms)', 'FontSize', 12); ylabel('Neuron index', 'FontSize', 12);
title(sprintf('Raster plot (p_I = %.2f, N = %d, Ne = %d, Ni = %d)', ...
      pI_target, rec.Ne+rec.Ni, rec.Ne, rec.Ni), 'FontSize', 14);
legend({'Excitatory spikes','Inhibitory spikes'}, 'Location','southeast', 'Orientation','vertical', 'FontSize', 12, 'Box', 'on');

grid off; box on;

%% --- Class-separated synaptic current traces ---
% These are population means of incoming synaptic currents per target group.
t_ms = (0:rec.steps-1) * dt_ms;

figure('Color','w'); 
plot(t_ms, rec.traces.mean_I_E_fromE, 'LineWidth',1.3); hold on;
plot(t_ms, rec.traces.mean_I_E_fromI, 'LineWidth',1.3);
plot(t_ms, rec.traces.mean_I_I_fromE, 'LineWidth',1.3);
plot(t_ms, rec.traces.mean_I_I_fromI, 'LineWidth',1.3);
grid on; box on;
xlabel('Time (ms)');
ylabel('Mean incoming synaptic current');
title(sprintf('Synaptic input by connection class (p_I = %.2f)', pI_target));
legend({'E \leftarrow E','E \leftarrow I','I \leftarrow E','I \leftarrow I'}, 'Location','best');

%% --- Combined figure: Raster + Synaptic currents ---
fig = figure('Color','w', 'Position', [100, 100, 1100, 400]);

% --- Left subplot: Raster plot ---
subplot(1,2,1); hold on;

if ~isempty(rec.raster.t)
    isE = rec.raster.i <= rec.Ne;   % E neuron indices
    isI = ~isE;
    plot(rec.raster.t(isE), rec.raster.i(isE), '.', ...
        'Color', [1.0 0.4 0.4], 'MarkerSize', 10); % Excitatory (red pastel)
    plot(rec.raster.t(isI), rec.raster.i(isI), '.', ...
        'Color', [0.2 0.6 1.0], 'MarkerSize', 10); % Inhibitory (blue pastel)
end

xlabel('Time (ms)', 'FontSize', 16); 
ylabel('Neuron index', 'FontSize', 16);

legend({'Excitatory spikes','Inhibitory spikes'}, ...
       'Location','northoutside', 'Orientation','horizontal', ...
       'FontSize', 13, 'Box','off');
box on; grid off;

% --- Right subplot: Synaptic current traces ---
subplot(1,2,2); hold on;
t_ms = (0:rec.steps-1) * dt_ms;

plot(t_ms, rec.traces.mean_I_E_fromE, 'LineWidth',1.5, 'Color',[0.1 0.5 0.1]); % dark green
plot(t_ms, rec.traces.mean_I_E_fromI, 'LineWidth',1.5, 'Color',[0.7 0.2 0.2]); % dark red
plot(t_ms, rec.traces.mean_I_I_fromE, 'LineWidth',1.5, 'Color',[0.2 0.4 0.8]); % dark blue
plot(t_ms, rec.traces.mean_I_I_fromI, 'LineWidth',1.5, 'Color',[0.8 0.4 0.8]); % purple

xlabel('Time (ms)', 'FontSize', 16);
ylabel('Mean synaptic current', 'FontSize', 16);

legend({'E \leftarrow E','E \leftarrow I','I \leftarrow E','I \leftarrow I'}, ...
       'Location','northoutside', 'Orientation','horizontal', ...
       'FontSize', 13, 'Box','off');
box on; grid on;

% --- Export in multiple formats ---
% print(fig, 'cortical_case_pi02', '-dpng', '-r600');
% print(fig, 'cortical_case_pi02', '-depsc');
% print(fig, 'cortical_case_pi02', '-dsvg');

%% ------------------------------------------------------------------------
function out = run_izh_efficiency_once( ...
    pI, N_total, T_sec, dt_ms, seed, ...
    p_conn, wE_mu, wE_sd, wE_min, wE_max, ...
    wI_mu, wI_sd, wI_min, wI_max, ...
    tau_syn_E, tau_syn_I, ...
    drive_amp, tonic_E, tonic_I, ...
    learning_rule, stdp, do_record)

    if nargin < 22, do_record = false; end
    rng(seed);

    % --- population sizes ---
    Ni = round(pI * N_total);
    Ne = N_total - Ni;
    N  = N_total;

    % --- Izhikevich RS (E) / FS (I) params ---
    a = [0.02*ones(Ne,1); 0.10*ones(Ni,1)];
    b = [0.20*ones(Ne,1); 0.20*ones(Ni,1)];
    c = [-65 *ones(Ne,1); -65 *ones(Ni,1)];
    d = [ 8.0*ones(Ne,1);  2.0*ones(Ni,1)];

    v = -65 * ones(N,1);
    u = b .* v;

    % --- Connectivity masks (store magnitudes only) ---
    % E -> *
    mask_E = rand(N, Ne) < p_conn;                 % postsyn x presyn(E)
    [rowE, colE] = find(mask_E);
    wE = wE_mu + wE_sd .* randn(numel(rowE), 1);   % magnitudes
    wE = min(max(wE, wE_min), wE_max);
    WE_abs = sparse(rowE, colE, wE, N, Ne);

    % I -> *
    mask_I = rand(N, Ni) < p_conn;                 % postsyn x presyn(I)
    [rowI, colI] = find(mask_I);
    wI = wI_mu + wI_sd .* randn(numel(rowI), 1);   % magnitudes
    wI = min(max(wI, wI_min), wI_max);
    WI_abs = sparse(rowI, colI, wI, N, Ni);

    % --- Sign matrix for full W: columns 1..Ne are +1 (E), Ne+1..N are -1 (I)
    W_sign = [ones(N, Ne), -ones(N, Ni)];
    build_W = @() zero_diag(W_sign .* [WE_abs, WI_abs]);

    % --- Presynaptic traces ---
    alpha_E = exp(-dt_ms / tau_syn_E);
    alpha_I = exp(-dt_ms / tau_syn_I);
    s       = zeros(N,1);
    alpha   = [alpha_E*ones(Ne,1); alpha_I*ones(Ni,1)];

    % --- Time ---
    steps = round(1000*T_sec / dt_ms);

    % --- Metrics ---
    energyE        = 0;
    spikesE_total  = 0;

    % --- Traces to record (only if requested) ---
    if do_record
        ras_t = []; ras_i = [];              % raster times/indices
        mean_I_E_fromE = zeros(steps,1);
        mean_I_E_fromI = zeros(steps,1);
        mean_I_I_fromE = zeros(steps,1);
        mean_I_I_fromI = zeros(steps,1);
    end

    % --- STDP bookkeeping ---
    do_plastic = ~strcmpi(learning_rule,'none');
    if do_plastic
        tau_plus   = stdp.tau_plus;
        tau_minus  = stdp.tau_minus;
        A_plus     = stdp.A_plus;
        A_minus    = stdp.A_minus;
        eta        = stdp.eta;

        decay_plus  = exp(-dt_ms / tau_plus);
        decay_minus = exp(-dt_ms / tau_minus);

        x_preE = zeros(Ne,1);
        x_preI = zeros(Ni,1);
        x_post = zeros(N,1);
    end

    % --- SIMULATION LOOP ---
    W = build_W(); % initial signed W

    for t = 1:steps
        % Background / tonic drive
        I_bg = drive_amp * (rand(N,1) - 0.5);
        I_bg(1:Ne)     = I_bg(1:Ne)     + tonic_E;
        I_bg(Ne+1:end) = I_bg(Ne+1:end) + tonic_I;

        % Synaptic current
        I_syn = W * s;

        % Total input
        I = I_bg + I_syn;

        % Energy proxy over E pop
        energyE = energyE + sum(abs(I(1:Ne))) * dt_ms;

        % Izhikevich updates
        v = v + 0.5*((0.04*v + 5).*v + 140 - u + I);
        v = v + 0.5*((0.04*v + 5).*v + 140 - u + I);
        u = u + a .* (0.2*v - u);

        % Spikes
        fired = find(v >= 30);
        if ~isempty(fired)
            v(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
            spikesE_total = spikesE_total + sum(fired <= Ne);
            if do_record
                ras_t = [ras_t; (t-1)*dt_ms * ones(numel(fired),1)]; %#ok<AGROW>
                ras_i = [ras_i; fired];                              %#ok<AGROW>
            end
        end

        % Presyn traces
        s = alpha .* s;
        if ~isempty(fired), s(fired) = s(fired) + 1.0; end

        % --- Class-separated syn currents (recording) ---
        if do_record
            sE = s(1:Ne);
            sI = s(Ne+1:end);

            % Build signed blocks on-the-fly (cheap vs readability)
            WEE = W(1:Ne,         1:Ne);
            WEI = W(1:Ne,   Ne+1:Ne+Ni);
            WIE = W(Ne+1:Ne+Ni,    1:Ne);
            WII = W(Ne+1:Ne+Ni, Ne+1:Ne+Ni);

            I_E_fromE = WEE * sE;  % input to E from E
            I_E_fromI = WEI * sI;  % input to E from I
            I_I_fromE = WIE * sE;  % input to I from E
            I_I_fromI = WII * sI;  % input to I from I

            mean_I_E_fromE(t) = mean(I_E_fromE);
            mean_I_E_fromI(t) = mean(I_E_fromI);
            mean_I_I_fromE(t) = mean(I_I_fromE);
            mean_I_I_fromI(t) = mean(I_I_fromI);
        end

        % --- STDP ---
        if do_plastic
            x_preE = x_preE * decay_plus;
            x_preI = x_preI * decay_plus;
            x_post = x_post * decay_minus;

            if ~isempty(fired)
                x_post(fired) = x_post(fired) + 1;

                firedE = fired(fired <= Ne);
                firedI = fired(fired >  Ne) - Ne;
                if ~isempty(firedE), x_preE(firedE) = x_preE(firedE) + 1; end
                if ~isempty(firedI), x_preI(firedI) = x_preI(firedI) + 1; end

                % --- Excitatory magnitude updates (columns 1..Ne) ---
                for jj = firedE.'
                    rows_i = find(WE_abs(:, jj));
                    if ~isempty(rows_i)
                        dWij = eta * A_plus * x_post(rows_i);
                        new_w = full(WE_abs(rows_i, jj)) + dWij;
                        WE_abs(rows_i, jj) = min(max(new_w, wE_min), wE_max);
                    end
                end
                for ii = fired.'
                    cols_j = find(WE_abs(ii, :));
                    if ~isempty(cols_j)
                        dWij = - eta * A_minus * x_preE(cols_j).';
                        new_w = full(WE_abs(ii, cols_j)) + dWij;
                        WE_abs(ii, cols_j) = min(max(new_w, wE_min), wE_max);
                    end
                end

                % --- Inhibitory magnitude updates (columns Ne+1..N) ---
                for jjI = firedI.'
                    rows_i = find(WI_abs(:, jjI));
                    if ~isempty(rows_i)
                        dWij = eta * A_plus * x_post(rows_i);
                        new_w = full(WI_abs(rows_i, jjI)) + dWij;
                        WI_abs(rows_i, jjI) = min(max(new_w, wI_min), wI_max);
                    end
                end
                for ii = fired.'
                    cols_jI = find(WI_abs(ii, :));
                    if ~isempty(cols_jI)
                        dWij = - eta * A_minus * x_preI(cols_jI).';
                        new_w = full(WI_abs(ii, cols_jI)) + dWij;
                        WI_abs(ii, cols_jI) = min(max(new_w, wI_min), wI_max);
                    end
                end
            end

            % Rebuild signed W and zero diagonal
            W = build_W();
        end
    end

    % Outputs
    out.Ne         = Ne;
    out.Ni         = Ni;
    out.steps      = steps;
    out.spikesE    = spikesE_total;
    out.energyE    = energyE;
    out.efficiency = spikesE_total / max(energyE, eps);

    if do_record
        out.raster.t = ras_t;
        out.raster.i = ras_i;
        out.traces.mean_I_E_fromE = mean_I_E_fromE;
        out.traces.mean_I_E_fromI = mean_I_E_fromI;
        out.traces.mean_I_I_fromE = mean_I_I_fromE;
        out.traces.mean_I_I_fromI = mean_I_I_fromI;
    else
        out.raster = struct('t',[],'i',[]);
        out.traces = struct();
    end
end

%% --- Helpers -------------------------------------------------------------
function Wz = zero_diag(W)
    n  = size(W,1);
    Wz = W - spdiags(diag(W), 0, n, n);
end
