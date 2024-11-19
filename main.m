% Alejandro Rodriguez-Garcia
% 20/05/24
% Based on Izhikevich, 2007 (daspnet.m, spnet.m).
% This script simulates a spiking neural network with dopamine-modulated 
% spike-timing dependent plasticity (DA-STDP) and regular STDP learning 
% rules.
% The neurons are divided into excitatory and inhibitory populations with
% different spiking types. The network dynamics include random thalamic 
% input and learning is driven by either STDP or DA-STDP. Synaptic 
% plasticity is updated at each timestep, and dopamine levels adjust 
% based on reward signals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETERS
close all; clear; clc;
seed = 2024; 
rng(seed);
c_date_time = datestr(clock,'YYYY-mm-dd_HH-MM');
hFig = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');

%% Network definition

% Simulation Parameters
N = 1000; % number of neurons
M = 100;  % number of synapses per neuron
D = 5;    % maximal conduction delay
% Time
T = 3;        % simulation time (in seconds)
sec_ms = 1000;   % (1s = 1000ms)

%% Neuron populations
% Define excitatory and inhibitory neuron populations
Ne = round(0.8*N); colorE = [1 0.7 0.4];  % 80% of neurons are excitatory (orange)
Ni = N - Ne; colorI = [0.6 0.8 1]; % 20% of neurons are inhibitory (blue)

% Spiking types
neuron_type_map = containers.Map(...
    {'RS', 'IB', 'CH', 'FS', 'TH', 'LTS', 'RZ', 'AC'}, ...
    {[0.02 0.20 -65 2], ...    % REGULAR SPIKING - RS
     [0.02 0.20 -55 4], ...    % INTRINSICALLY BURSTING - IB
     [0.02 0.20 -50 2], ...    % CHATTERING - CH
     [0.10 0.20 -65 2], ...    % FAST SPIKING - FS
     [0.02 0.25 -65 0.05], ... % THALAMO-CORTICAL - TH
     [0.02 0.25 -65 2], ...    % LOW-THRESHOLD SPIKING - LTS
     [0.10 0.26 -65 2], ...    % RESONATOR - RZ
     [0.02 1.00 -55 4]});      % ACCOMODATING - AC
% populations
popE1 = 'RS'; popE2 = 'RS'; popI1 = 'FS'; popI2 = 'FS'; 
popE1_params = neuron_type_map(popE1);
popE2_params = neuron_type_map(popE2);
popI1_params = neuron_type_map(popI1);
popI2_params = neuron_type_map(popI2);

% (a, b, c, d) parameters
a = [popE1_params(1)*ones(Ne*0.8,1); popE2_params(1)*ones(Ne*0.2,1); ...
     popI1_params(1)*ones(Ni*0.9,1); popI2_params(1)*ones(Ni*0.1,1)]; % recovery parameter
b = [popE1_params(2)*ones(Ne*0.8,1); popE2_params(2)*ones(Ne*0.2,1); ...
     popI1_params(2)*ones(Ni*0.9,1); popI2_params(2)*ones(Ni*0.1,1)]; % recovery sensitivity
c = [popE1_params(3)*ones(Ne*0.8,1); popE2_params(3)*ones(Ne*0.2,1); ...
     popI1_params(3)*ones(Ni*0.9,1); popI2_params(3)*ones(Ni*0.1,1)]; % recovery reset after a spike
d = [popE1_params(4)*ones(Ne*0.8,1); popE2_params(4)*ones(Ne*0.2,1); ...
     popI1_params(4)*ones(Ni*0.9,1); popI2_params(4)*ones(Ni*0.1,1)]; % recovery reset after a spike

%% Synapses
W = ones(N, M); % synaptic weights initialization
signs = [ones(Ne, M); -ones(Ni, M)]; % signs of synaptic weights
dW = zeros(N, M); % derivatives of synaptic strengths
Wm = 4; % maximal synaptic strength

post = ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); % post-synaptic neurons
delays = cell(N, D); % conduction delay containers
pre = cell(N, 1); % pre-synaptic indexes
aux = cell(N, 1); % auxiliary variable for delay adjustment
for i=1:N
    if i<=Ne
        for j=1:M
            delays{i, ceil(D*rand)}(end+1) = j; % Assign random exc delays
        end
    else
        delays{i,1}=1:M; % fixed delay for inhibitory synapses
    end
    pre{i}=find(post==i); % both E and I pre-synaptic indexes
    aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end

%% DA-STDP (Dopamine-modulated Spike Timing Dependent Plasticity)
% Synapse to be reinforced
n1 = ceil(Ne*rand);  % presynaptic neuron
syn = ceil(M*rand);  % synapse number to the postsynaptic neuron
n2 = post(n1,syn);   % postsynaptic neuron
W(n1,syn) = 0;       % set initial weight to zero

DA = 0;   % dopamine level above baseline
rew = []; % reward events
STDP = zeros(N,sec_ms+1+D); % STDP traces
firings = [-D 0]; % spike timings
interval = 20;    % the coincidence interval for n1 and n2
n1f = -D;         % the last spike of n1
n2f = [];         % the last spike of n2
shist = zeros(sec_ms*T,2); % synaptic strength
DAv = zeros(sec_ms*T,1); % DA release

%% Initializations
v = -65 * ones(N, 1); % initial membrane potential - IZ
u = b .* v; % membrane recovery variable - IZ

%% Learning rule selection
learning_rule = 'DA-stdp'; % change this to 'stdp' for STDP learning rule

%% SIMULATION

for sec = 1:T % simulation of 1 day
    for t = 1:sec_ms % simulation of 1 sec
        I = 13*(rand(N, 1) - 0.5); % random thalamic input

        fired = find(v >= 30); % indices of fired neurons        
        v(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
        STDP(fired, t+D) = 0.1;
        for k = 1:length(fired)
            dW(pre{fired(k)}) = dW(pre{fired(k)}) + STDP(N*t + aux{fired(k)});
        end
        firings = [firings; t*ones(length(fired), 1), fired];
        k = size(firings, 1);
        while firings(k, 1) > t-D
            del = delays{firings(k, 2), t - firings(k, 1) + 1};
            ind = post(firings(k, 2), del);
            I(ind) = I(ind) + (signs(firings(k, 2), del) .* W(firings(k, 2), del))';
            dW(firings(k, 2), del) = dW(firings(k, 2), del) - 1.5 * STDP(ind, t+D)';
            k = k - 1;
        end
        v = v + 0.5 * ((0.04*v + 5).*v + 140 - u + I); % membrane voltage variable
        v = v + 0.5 * ((0.04*v + 5).*v + 140 - u + I); % stability step
        u = u + a .* (0.2*v - u);                      % recovery variable

        STDP(:, t+D+1) = 0.95 * STDP(:, t+D);          % STDP trace decay (tau = 20 ms)

        % Learning
        if (mod(t, 10) == 0)
            if strcmp(learning_rule, 'DA-stdp')
                W = max(0,min(Wm,W + (0.002+DA)*dW)); % DA-STDP
            else
                W = max(0, min(Wm, W + 0.002 * dW)); % STDP
            end
            dW = 0.99 * dW;
        end

        % Dopamine update
        if any(fired==n1)
            n1f=[n1f,sec*sec_ms+t];
        end
        if any(fired==n2)
            n2f=[n2f,sec*sec_ms+t];
            if (sec*sec_ms+t-n1f(end)<interval) && (n2f(end)>n1f(end))
                rew=[rew,sec*sec_ms+t+sec_ms+ceil(2000*rand)];
            end
        end
        
        if any(rew==sec*sec_ms+t)
            DA=DA+0.5; % dopamine release if reward
        end
        DA=DA*0.99; % dopamine decay

        DAv(sec*sec_ms+t) = DA;
        shist(sec*sec_ms+t,:)=[W(n1,syn),dW(n1,syn)];
    end

    % ---- plot -------
    % Spike rasterplot
    excitatoryFirings = firings(firings(:,2) <= Ne, :);
    inhibitoryFirings = firings(firings(:,2) > Ne, :);
    subplot(3,2,[1 2],'Parent', hFig);
    plot(excitatoryFirings(:,1), excitatoryFirings(:,2), '.', 'MarkerSize', 10, 'Color', colorE, 'DisplayName', 'Excitatory firings');
    hold on;
    plot(inhibitoryFirings(:,1), inhibitoryFirings(:,2), '.', 'MarkerSize', 10, 'Color', colorI, 'DisplayName', 'Inhibitory firings');
    plot(firings(firings(:,2) == n1, 1), firings(firings(:,2) == n1, 2), 'k+', 'MarkerSize', 15, 'DisplayName', 'Pre-synaptic neuron'); % Mark neuron n1
    plot(firings(firings(:,2) == n2, 1), firings(firings(:,2) == n2, 2), 'ks', 'MarkerSize', 15, 'DisplayName', 'Post-synaptic neuron'); % Mark neuron n2
    hold off;
    axis([0 sec_ms 0 N]);
    xlabel('Time (ms)');
    ylabel('Neurons'); 
    legend('show', 'Location', 'northeastoutside');
    % Synaptic strength, STDP trace, DA and rewards
    subplot(3, 2, [5 6],'Parent', hFig);
    plot(0.001*(1:(sec*sec_ms+t)), shist(1:sec*sec_ms+t, 1), 'b', 'DisplayName', 'Synaptic strength');
    hold on;
    plot(0.001*(1:(sec*sec_ms+t)), shist(1:sec*sec_ms+t, 2), 'g', 'DisplayName', 'STDP trace');
    plot(0.001*(1:(sec*sec_ms+t)), DAv(1:sec*sec_ms+t), 'k', 'DisplayName', 'DA release');
    plot(0.001*rew, 0*rew, 'rx', 'DisplayName', 'rewards', 'MarkerSize', 20);
    hold off;
    xlabel('Time (s)');
    ylabel('Value');
    legend('show', 'Location', 'northeastoutside');
    % Synaptic weight distributions
    subplot(3, 2, [3 4],'Parent', hFig);
    histogram(signs(1:Ne, :) .* W(1:Ne, :), 'FaceColor', colorE, 'EdgeColor', colorE, 'DisplayName', 'Excitatory distribution'); % excitatory synapses
    hold on;
    histogram(signs(Ne+1:end, :) .* W(Ne+1:end, :), 'FaceColor', colorI, 'EdgeColor', colorI, 'DisplayName', 'Inhibitory distribution'); % inhibitory synapses
    hold on;
    plot(signs(n1,syn).*W(n1,syn), 0, 'b.', 'MarkerSize', 20, 'DisplayName', 'Reinforced synapse');
    hold off;
    xlabel('Synaptic weight');
    ylabel('Count');
    legend('show', 'Location', 'northeastoutside');
    timeText = text(1.25, -1.5, sprintf('t = %d sec', sec), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 12);    
    drawnow;
    % ---------------


    STDP(:, 1:D+1) = STDP(:, sec_ms+1:sec_ms+1+D);
    ind = find(firings(:, 1) > sec_ms + 1 - D);
    firings = [-D 0; firings(ind, 1) - sec_ms, firings(ind, 2)];
    disp(sec); % display current second of simulation
end

% Save file
filename = sprintf('shist_Epop%s_%s_Ipop%s_%s_%s_%s.mat', popE1, popE2, popI1, popI2, learning_rule, c_date_time);
save(filename, 'shist');
