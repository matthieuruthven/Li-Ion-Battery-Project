% Path to folder containing all simulation results
dirpath = 'C:\Users\matthieu.ruthven\Documents\BatteryProject';

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_1_00_T_amb_23_85_h_12_0_k_0_9_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
A = surface_t_array(27:end, :, :);

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_1_00_T_amb_23_85_h_12_0_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
B = surface_t_array(27:end, :, :);

% Ensure that arrays have the same dimensions
if size(A, 3) ~= size(B, 3)
    A_n = size(A, 3);
    B_n = size(B, 3);
    if A_n < B_n
        B = B(:,:, 1:A_n);
    else
        A = A(:, :, 1:B_n);
    end
end

% Calculate maximum difference between surface temperature values
T_diff = A - B;
T_diff_max = squeeze(max(T_diff, [], [1 2]));

% Create time series
ts1 = timeseries(T_diff_max, 1:size(T_diff_max, 1));
ts1.Name = 'Maximum difference in surface temperature (degrees C)';
ts1.TimeInfo.Units = 'seconds';
plot(ts1)
grid on

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_12_0_k_0_9_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
A = surface_t_array(27:end, :, :);

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
B = surface_t_array(27:end, :, :);

% Ensure that arrays have the same dimensions
if size(A, 3) ~= size(B, 3)
    A_n = size(A, 3);
    B_n = size(B, 3);
    if A_n < B_n
        B = B(:,:, 1:A_n);
    else
        A = A(:, :, 1:B_n);
    end
end

% Calculate maximum difference between surface temperature values
T_diff = A - B;
T_diff_max = squeeze(max(T_diff, [], [1 2]));

% Create time series
hold on
ts2 = timeseries(T_diff_max, 1:size(T_diff_max, 1));
ts2.Name = 'Maximum difference in surface Temperature (degrees C)';
ts2.TimeInfo.Units = 'seconds';
plot(ts2)

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_1_00_T_amb_23_85_h_12_0_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
A = surface_t_array(27:end, :, :);

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_1_00_T_amb_23_85_h_12_0_k_1_3_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
B = surface_t_array(27:end, :, :);

% Ensure that arrays have the same dimensions
if size(A, 3) ~= size(B, 3)
    A_n = size(A, 3);
    B_n = size(B, 3);
    if A_n < B_n
        B = B(:,:, 1:A_n);
    else
        A = A(:, :, 1:B_n);
    end
end

% Calculate maximum difference between surface temperature values
T_diff = A - B;
T_diff_max = squeeze(max(T_diff, [], [1 2]));

% Create time series
hold on
ts3 = timeseries(T_diff_max, 1:size(T_diff_max, 1));
ts3.Name = 'Maximum difference in surface temperature (degrees C)';
ts3.TimeInfo.Units = 'seconds';
plot(ts3)

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
A = surface_t_array(27:end, :, :);

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_12_0_k_1_3_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
B = surface_t_array(27:end, :, :);

% Ensure that arrays have the same dimensions
if size(A, 3) ~= size(B, 3)
    A_n = size(A, 3);
    B_n = size(B, 3);
    if A_n < B_n
        B = B(:,:, 1:A_n);
    else
        A = A(:, :, 1:B_n);
    end
end

% Calculate maximum difference between surface temperature values
T_diff = A - B;
T_diff_max = squeeze(max(T_diff, [], [1 2]));

% Create time series
hold on
ts4 = timeseries(T_diff_max, 1:size(T_diff_max, 1));
ts4.Name = 'Maximum difference in surface temperature (degrees C)';
ts4.TimeInfo.Units = 'seconds';
plot(ts4)

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_9_6_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
A = surface_t_array(27:end, :, :);

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
B = surface_t_array(27:end, :, :);

% Ensure that arrays have the same dimensions
if size(A, 3) ~= size(B, 3)
    A_n = size(A, 3);
    B_n = size(B, 3);
    if A_n < B_n
        B = B(:,:, 1:A_n);
    else
        A = A(:, :, 1:B_n);
    end
end

% Calculate maximum difference between surface temperature values
T_diff = A - B;
T_diff_max = squeeze(max(T_diff, [], [1 2]));

% Create time series
hold on
ts5 = timeseries(T_diff_max, 1:size(T_diff_max, 1));
ts5.Name = 'Maximum difference in surface temperature (degrees C)';
ts5.TimeInfo.Units = 'seconds';
plot(ts5)

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
A = surface_t_array(27:end, :, :);

% Path to MAT file of temperature values
filepath = fullfile(dirpath, 'C_rate_4_00_T_amb_23_85_h_14_4_k_1_1_Discharge', 'Surface_T.mat');

% Read file
load(filepath)

% Extract surface temperature of battery body only (i.e. exclude tabs)
B = surface_t_array(27:end, :, :);

% Ensure that arrays have the same dimensions
if size(A, 3) ~= size(B, 3)
    A_n = size(A, 3);
    B_n = size(B, 3);
    if A_n < B_n
        B = B(:,:, 1:A_n);
    else
        A = A(:, :, 1:B_n);
    end
end

% Calculate maximum difference between surface temperature values
T_diff = A - B;
T_diff_max = squeeze(max(T_diff, [], [1 2]));

% Create time series
hold on
ts6 = timeseries(T_diff_max, 1:size(T_diff_max, 1));
ts6.Name = 'Maximum difference in surface temperature (degrees C)';
ts6.TimeInfo.Units = 'seconds';
plot(ts6)

% Change title
title('')
% xlabel('Time (seconds)')
% ylabel('Maximum difference in surface temperature (degrees C)')

% Add legend
legend('C-rate = 1, k = 0.9 vs 1.1, h = 12.0', ...
    'C-rate = 4, k = 0.9 vs 1.1, h = 12.0', ...
    'C-rate = 1, k = 1.1 vs 1.3, h = 12.0', ...
    'C-rate = 4, k = 1.1 vs 1.3, h = 12.0', ...
    'C-rate = 4, k = 1.1, h = 9.6 vs 12.0', ...
    'C-rate = 4, k = 1.1, h = 12.0 vs 14.4', ...
    'Location', 'northeast')