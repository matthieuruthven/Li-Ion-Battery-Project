function Energetic_Generate_Dataset(sim_type, sim_time, C_rate, k, h, dirpath, sq_wave_period, img_syn_int)
% Energetic_Generate_Dataset.m
% Function to generate dataset of synthetic thermal images using
% computational model developed by Lin et al. (2022)
% https://www.nature.com/articles/s44172-022-00005-8#Abs1
% 
% Function includes code adapted from:
% https://github.com/Battery-Intelligence-Lab/multiscale-coupling
%
% Input arguments:
% sim_type (string): specifies whether the pouch cell is charged, 
% discharged or both in the simulation (options: Charge, Discharge or 
% Square)
% sim_time (integer): the duration of the simulation in seconds
% C_rate (float): the C-rate at which the pouch cell is (dis)charged at
% k (float): the thermal conductivity of the pouch cell in W/m/K
% h (float): the heat transfer coefficient of the pouch cell in W/m^2/K
% dirpath (string): the path to the folder in which the dataset should be
% saved
% sq_wave_period (integer): the duration of the square wave in seconds (only
% relevant when the sim_type argument is 'Square')
% img_syn_int (integer): the time interval between synthesis of
% consecutive images

% Author: Matthieu Ruthven (matthieu.ruthven@uni.lu)
% Last modified: 22nd February 2024

% Start timing
tic

% If sim_type is either Charge or Discharge, sq_wave_period should be 2 * 
% sim_time
if strcmp(sim_type, 'Charge') || strcmp(sim_type, 'Discharge')
    sq_wave_period = 2 * sim_time;
elseif strcmp(sim_type, 'Square')
    if sim_time < sq_wave_period / 2
        error('When the "sim_type" argument is "Square", the "sim_time" argument should be greater than half the "sq_wave_period" argument')
    end
else
    error('"sim_type" argument should be "Charge", "Discharge" or "Square"')
end

% Vector of battery model parameter values
% Values are those from Table S3 in Lin et al. (2022)
x0 = [1.0220    1.1000     12.0000     0.1800     81.0000       4.2000     6.0000     29.5000   4.5000];
x = x0;
% NB
% x0(1) is 2.35 * 10^6 / 2300 [Cp in Table S3] / [electrode density]
% x0(4) is ku 70% SOC in Table 1
% x0(5) is 0.046 * 42^2 [kappa in Table S3] * [number of cell layers]^2
% x0(6) is 2.4 * 10^-3 * 42^2 [alpha in Table S3] * [number of cell layers]^2

% Extract battery model parameter values from vector x
% Code source: https://github.com/Battery-Intelligence-Lab/multiscale-coupling/blob/dba3649f2eeb4228d2849e4d4e748c91ed6ad7a4/ObjFunc.m#L5C1-L17C48
tf          = round(sim_time);       % Total simulation time 
Cp          = x(1)*1e3;   % Specific heat
% k           = x(2);       % Thermal conductivity
% h           = x(3);       % Heat transfer coefficient   
ku          = x(4);       % OCP gradient
kappa       = x(5);       % Ionic conductivity
alfa_ka     = x(6);       % Temperature coefficient of kappa
i0ref       = x(7);       % Exchange current density
Ei0         = x(8);       % Reaction activation energy
Ds          = x(9)*1e-14; % Diffusion coeffcient
DS          = -0.07;      % Entropy change
Vinit       = 3.2381;      % Initial cell voltage
T0          = 23.85;      % Ambient temperature

% Create folder in which dataset will be saved (an error will be raised if
% the folder already exists)
if strcmp(sim_type, 'Charge') || strcmp(sim_type, 'Discharge')
    dirpath = fullfile(dirpath, strrep(sprintf('C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s', C_rate, T0, h, k, sim_type), '.', '_'));
else
    dirpath = fullfile(dirpath, strrep(sprintf('C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s_%d', C_rate, T0, h, k, sim_type, round(sq_wave_period)), '.', '_'));
end
[status, msg] = mkdir(dirpath);
if or(status ~= 1, contains(msg,'already')) 
    error(msg)
end

% Send battery model to COMSOL server
% Code source: https://github.com/Battery-Intelligence-Lab/multiscale-coupling/blob/dba3649f2eeb4228d2849e4d4e748c91ed6ad7a4/ObjFunc.m#L20
out = Battery_model_energetic(sim_type,round(sq_wave_period),round(img_syn_int),C_rate,dirpath,tf,T0,Cp,k,h,ku,kappa,alfa_ka,i0ref,Ei0,Ds,DS,Vinit);

% Finish timing
toc
