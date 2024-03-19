function Energetic_Generate_Dataset(sim_type, sim_time, batt_type, anomaly, C_rate, k_in, h_in, T0_in, dirpath, sq_wave_period, img_syn_int, SOC, def_or_cus)
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
% batt_type (string): specifies the configuration and dimensions of the
% battery (options: Energetic or Lin)
% anomaly (string): specifies if there should be an anomaly in the thermal
% conductivity of the battery (options: None, Circle, Stripe_W, Stripe_H)
% C_rate (float): the C-rate at which the pouch cell is (dis)charged at
% k_in (float): the thermal conductivity of the pouch cell in W/m/K
% h_in (float): the heat transfer coefficient of the pouch cell in W/m^2/K
% T0_in (float): the ambient temperature and the initial temperature of the
% pouch cell
% dirpath (string): the path to the folder in which the dataset should be
% saved
% sq_wave_period (integer): the duration of the square wave in seconds (only
% relevant when the sim_type argument is 'Square')
% img_syn_int (integer): the time interval between synthesis of
% consecutive images
% SOC (float): the initial state of charge (SOC) of the positive electrode
% (options: 0.3, 0.5 or 0.7)
% def_or_cus (string): specifies whether to use custom values of k, h and 
% T0 (i.e. the values of the k, h and T0 arguments input to the function) 
% or default ones (i.e. the k, h and T0 arguments input to the function 
% will be overridden) (options: custom or default)

% Author: Matthieu Ruthven (matthieu.ruthven@uni.lu)
% Last modified: 18th March 2024

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

% Check battery
if strcmp(batt_type, 'Energetic') || strcmp(batt_type, 'Lin')
else
    error('"batt_type" argument should be either "Energetic" or "Lin"')
end

% Check anomaly
if strcmp(anomaly, 'None') || strcmp(anomaly, 'Circle')
elseif strcmp(anomaly, 'Stripe_H') || strcmp(anomaly, 'Stripe_W')
else
    error('"anomaly" argument should be "None", "Circle", "Stripe_H" or "Stripe_W"')
end

% Battery model parameter values
% NB
% 2.35 * 10^6 / 2300 [Cp in Table S3] / [electrode density]
% ku 70% SOC in Table 1
% 0.046 * 42^2 [kappa in Table S3] * [number of cell layers]^2
% 2.4 * 10^-3 * 42^2 [alpha in Table S3] * [number of cell layers]^2


tf          = round(sim_time);       % Total simulation time 
Cp          = 1022;                  % Specific heat (derived from Table S6)
i0ref       = 6;                     % Exchange current density (From Table S6)
Vinit       = 3.2381;                % Initial cell voltage (from original code)
if strcmp(sim_type, 'Square')
    if SOC == 0.3
        k = 1.1;         % Thermal conductivity
        h = 12.0;        % Heat transfer coefficient
        ku = 0.35;       % OCP gradient
        kappa = 81.1;    % Ionic conductivity
        alfa_ka = 4.2;   % Temperature coefficient of kappa
        Ei0 = 29.5;      % Reaction activation energy
        Ds = 4.5*1e-14;  % Diffusion coeffcient
        DS = -0.14;      % Entropy change
        T0 = 23.85;      % Ambient temperature
    elseif SOC == 0.5
        k = 1.2;         % Thermal conductivity
        h = 11.9;        % Heat transfer coefficient
        ku = 0.24;       % OCP gradient
        kappa = 75.9;    % Ionic conductivity
        alfa_ka = 4.8;   % Temperature coefficient of kappa
        Ei0 = 29.2;      % Reaction activation energy
        Ds = 4.2*1e-14;  % Diffusion coeffcient
        DS = 0.08;       % Entropy change
        T0 = 21.95;      % Ambient temperature
    elseif SOC == 0.7
        k = 1.1;         % Thermal conductivity
        h = 11.8;        % Heat transfer coefficient 
        ku = 0.18;       % OCP gradient
        kappa = 79.4;    % Ionic conductivity
        alfa_ka = 4.2;   % Temperature coefficient of kappa
        Ei0 = 29.7;      % Reaction activation energy
        Ds = 4*1e-14;    % Diffusion coeffcient
        DS = 0.10;       % Entropy change
        T0 = 25.35;      % Ambient temperature
    else
        error('"SOC" argument should be 0.3, 0.5 or 0.7')
    end
else
    k = 1.1;         % Thermal conductivity
    h = 12.0;        % Heat transfer coefficient
    ku = 0.35;       % OCP gradient (NB will be overridden)
    kappa = 81.1;    % Ionic conductivity
    alfa_ka = 4.2;   % Temperature coefficient of kappa
    Ei0 = 29.5;      % Reaction activation energy
    Ds = 4.5*1e-14;  % Diffusion coeffcient
    DS = -0.14;      % Entropy change (NB will be overridden)
    T0 = 23.85;      % Ambient temperature
end

% If required, change values of k, h and T0
if strcmp(def_or_cus, 'custom')
    k = k_in;         % Thermal conductivity
    h = h_in;         % Heat transfer coefficient
    T0 = T0_in;       % Ambient temperature
elseif strcmp(def_or_cus, 'default')
else
    error('"def_or_cus" argument should be either "custom" or "default"')
end

% Create folder in which dataset will be saved (an error will be raised if
% the folder already exists)
if strcmp(sim_type, 'Charge') || strcmp(sim_type, 'Discharge')
    if strcmp(anomaly, 'None')
        dirpath = fullfile(dirpath, strrep(sprintf('%s_C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s_SOC_%.01f', batt_type, C_rate, T0, h, k, sim_type, SOC), '.', '_'));
    else
        dirpath = fullfile(dirpath, strrep(sprintf('%s_C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s_SOC_%.01f_Anomaly_%s', batt_type, C_rate, T0, h, k, sim_type, SOC, anomaly), '.', '_'));
    end    
else
    if strcmp(anomaly, 'None')
        dirpath = fullfile(dirpath, strrep(sprintf('%s_C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s_%d_SOC_%.01f', batt_type, C_rate, T0, h, k, sim_type, round(sq_wave_period), SOC), '.', '_'));
    else
        dirpath = fullfile(dirpath, strrep(sprintf('%s_C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s_%d_SOC_%.01f_Anomaly_%s', batt_type, C_rate, T0, h, k, sim_type, round(sq_wave_period), SOC, anomaly), '.', '_'));
    end
end
[status, msg] = mkdir(dirpath);
if or(status ~= 1, contains(msg,'already')) 
    error(msg)
end

% Send battery model to COMSOL server
out = Battery_model_energetic(sim_type,round(sq_wave_period),round(img_syn_int),batt_type,anomaly,C_rate,dirpath,tf,T0,Cp,k,h,ku,kappa,alfa_ka,i0ref,Ei0,Ds,DS,Vinit,SOC);

% Finish timing
toc