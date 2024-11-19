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
files = {'shist_EpopRS_RS_IpopFS_FS_DA-stdp_2024-10-18_13-24.mat', ...
         'shist_EpopRS_CH_IpopFS_FS_DA-stdp_2024-10-18_14-01.mat',... % RS
         % 'ch-shist_EpopRS_CH_IpopFS_FS_DA-stdp_2024-10-18_15-05.mat',... % BU
         };
% Define labels for each file to use in the plot
labels = {'E:RS, I:FS', ...
          'E:0.8*RS/0.2*BU, I:FS', ... % pre-and post- RS neurons
          % 'E:0.8*RS/0.2*BU, I:FS', ... % pre-and post- BU neurons
          };
% Define colors for each dataset (lighter versions of the colors)
colors = {[102, 153, 255]/255, ... 
          [255, 80, 80]/255, ... % pre-and post- RS neurons
          % [192, 0, 0]/255, ... % pre-and post- BU neurons
          };

% Define general parameters
T = 3600;        % Number of sequences
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
% figure('Position', [100, 100, 1200, 600]); % Width 800px, Height 600px
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
         'DisplayName', [labels{i}]);
end

% Customize the plot
xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial'); % X-axis label
ylabel('Synaptic strength', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial'); % Y-axis label
legend('show', 'Location', 'northoutside', 'FontSize', 14, 'FontName', 'Arial', 'NumColumns', 3, 'Box', 'off'); 
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'FontName', 'Arial'); % Customize axes

% Remove grid lines from the plot
set(gca, 'XGrid', 'off', 'YGrid', 'off');

hold off; % Release the plot

% Save the figure as a PNG file
% saveas(gcf, 'synaptic_strength_comparison.png');