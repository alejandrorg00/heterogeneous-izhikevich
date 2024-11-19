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

%% Define general parameters
T = 3600;        % Number of sequences (1 hour)
sec_ms = 1000;   % Milliseconds per second

%% Plot 2: Compare second set of files
% Define the files containing synaptic histogram data for the second plot
files2 = {'shist_EpopRS_RS_IpopFS_AC_DA-stdp_2024-10-18_13-47.mat', ...
          'new-shist_EpopRS_RS_IpopFS_AC_DA-stdp_2024-10-25_09-32.mat', ...
          'new73-shist_EpopRS_RS_IpopFS_AC_DA-stdp_2024-10-25_09-37.mat', ...
          'new37-shist_EpopRS_RS_IpopFS_AC_DA-stdp_2024-10-25_09-40.mat', ...
          'shist_EpopRS_RS_IpopFS_FS_DA-stdp_2024-10-18_13-24.mat', ...
          '05-shist_EpopRS_RS_IpopAC_AC_DA-stdp_2024-10-18_16-34.mat'};

% Define labels for each file to use in the plot
labels2 = {'E:RS, I:0.9*FS/0.1*AC', ...
           'E:RS, I:0.5*FS/0.5*AC', ...
           'E:RS, I:0.7*FS/0.3*AC', ...
           'E:RS, I:0.3*FS/0.7*AC', ...
           'E:RS, I:FS', ...
           'E:RS, I:AC'};

% Define colors for each dataset
% Define colors for each dataset
colors2 = {[147, 112, 219]/255, ...       
           [128, 0, 128]/255, ...       
           [255, 140, 0]/255, ...       
           [0, 163, 155]/255, ...       
           [102, 153, 255]/255, ...     
           [51, 51, 255]/255};          


% Load all data for the second plot
data_structs2 = cell(length(files2), 1);
for i = 1:length(files2)
    data_structs2{i} = load(files2{i}); % Load each file and store it
end

% Initialize variables to store learning rates
learning_rates2 = zeros(length(files2), 1);

% Create a new figure for the second comparative plot
figure(2); % Width 1200px, Height 600px
hold on; % Keep the plot open to add multiple traces

% Iterate over each loaded file and calculate the learning rate
for i = 1:length(files2)
    % Extract the relevant variable (synaptic histogram)
    shist_data = data_structs2{i}.shist;
    
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
    learning_rates2(i) = p(1); % Store the slope as the learning rate
    disp(['Learning rate (' labels2{i} '): ', num2str(learning_rates2(i))]); % Display the learning rate

    % Plot the synaptic strength data
    plot(time_vector_full, shist_data(:, 1), 'Color', colors2{i}, 'LineWidth', 2.5, ...
         'DisplayName', [labels2{i}]);
end

% Customize the plot
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial'); % X-axis label
ylabel('Synaptic strength', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial'); % Y-axis label
legend('show', 'Location', 'northoutside', 'FontSize', 14, 'FontName', 'Arial', 'NumColumns', 3, 'Box', 'off'); 
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Arial'); % Customize axes
set(gca, 'XGrid', 'off', 'YGrid', 'off');
hold off; % Release the plot
