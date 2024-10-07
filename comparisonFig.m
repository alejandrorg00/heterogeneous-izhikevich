% Alejandro Rodriguez-Garcia
% 20/05/24
% Analysis of results for neural heterogeneity
%
% This script loads synaptic plasticity data from different neural 
% population configurations, calculates the learning rate for each
% configuration by fitting a line to the synaptic strength over time, and 
% plots the synaptic strength for comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;

%% Load files
% Define the files containing synaptic histogram data
files = {'shist_EpopRS_IpopFS_stdp.mat', 'shist_EpopRS_IpopRS_stdp.mat', ...
         'shist_EpopRS_IpopFS_DA-stdp.mat', 'shist_EpopRS_IpopRS_DA-stdp.mat'};
% Define labels for each file to use in the plot
labels = {'E:RS, I:FS (stdp)', 'E:RS, I:RS (stdp)', ...
          'E:RS, I:FS (DA-stdp)', 'E:RS, I:RS (DA-stdp)'};
% Define colors for each dataset (lighter versions of the colors)
colors = {[198, 12, 48]/255, [255, 179, 186]/255, [0, 163, 155]/255, [153, 217, 234]/255};

% Define general parameters
T = 3600;        % Number of sequences (1 hour)
sec_ms = 1000;   % Milliseconds per second

% Load all data
% Initialize a cell to store loaded data
data_structs = cell(length(files), 1);
for i = 1:length(files)
    data_structs{i} = load(files{i}); % Load each file and store it
end

% Initialize variables to store learning rates
learning_rates = zeros(length(files), 1);

% Create a new figure for the comparative plot with a larger size
figure('Position', [100, 100, 800, 600]); % Width 800px, Height 600px
hold on; % Keep the plot open to add multiple traces

% Iterate over each loaded file and calculate the learning rate
for i = 1:length(files)
    % Extract the relevant variable (synaptic histogram)
    shist_data = data_structs{i}.shist;
    
    % Create a time vector for the full duration
    time_vector_full = 0.001 * (1:length(shist_data))'; % Convert to seconds

    % Find the first moment when shist_data(:, 1) reaches Wm = 4
    idx_end = find(shist_data(:, 1) >= 4, 1, 'first');
    if isempty(idx_end)
        idx_end = length(shist_data(:, 1)); % If it does not reach 4, use the last index
    end

    % Calculate synaptic strength for shist(:, 1) up to the found index
    synaptic_strength = shist_data(1:idx_end, 1);
    time_vector = time_vector_full(1:idx_end);
    % Fit a straight line to the data to calculate the learning rate
    p = polyfit(time_vector, synaptic_strength, 1);
    learning_rates(i) = p(1); % Store the slope as the learning rate
    disp(['Learning rate (' labels{i} '): ', num2str(learning_rates(i))]); % Display the learning rate

    % Plot the synaptic strength data
    plot(time_vector_full, shist_data(:, 1), 'Color', colors{i}, 'LineWidth', 2.5, ...
         'DisplayName', [labels{i} ' (LR: ' num2str(learning_rates(i)) ')']);
end

% Customize the plot
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial'); % X-axis label
ylabel('Synaptic strength', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial'); % Y-axis label
title('Synaptic Strength Over Time', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial'); % Plot title
legend('show', 'Location', 'northoutside', 'FontSize', 10, 'FontName', 'Arial'); % Show legend at the top
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'FontName', 'Arial'); % Customize axes

% Remove grid lines from the plot
set(gca, 'XGrid', 'off', 'YGrid', 'off');

hold off; % Release the plot

% Save the figure as a PNG file
saveas(gcf, 'synaptic_strength_comparison.png');