function Energetic_Generate_Dataset(sim_type, C_rate, dirpath)
% Energetic_Generate_Dataset.m
% Function to generate dataset of synthetic thermal images using
% computational model developed by Lin et al. (2022)
% https://www.nature.com/articles/s44172-022-00005-8#Abs1
% 
% Function includes code adapted from:
% https://github.com/Battery-Intelligence-Lab/multiscale-coupling
%
% Input arguments:
% C_rate (float): the C-rate at which the pouch cell is (dis)charged at
% dirpath (string): the path to the folder in which the dataset should be
% saved

% Author: Matthieu Ruthven (matthieu.ruthven@uni.lu)
% Last modified: 7th February 2024

% Start timing
tic

% Array of required simulation times (specific to Lin et al. (2022) pouch cell)
if strcmp(sim_type, 'Charge')
    sim_time = [6800, 3400, 2400, 1800]; % [floor(C_rate=1), floor(C_rate=2), ..., floor(C_rate=3)]
elseif strcmp(sim_type, 'Discharge')
    sim_time = [5400, 2700, 1800, 1800]; % [floor(C_rate=1), floor(C_rate=2), ..., floor(C_rate=3)]
else
    error('"sim_type" argument should be either "Charge" or "Discharge"')
end

% Vector of battery model parameter values
% Values are those from Table S3 in Lin et al. (2022)
x0 = [1.0220    1.1000     12.0000     0.3500     81.0000       4.2000     6.0000     29.5000   4.5000];
x = x0;

% Extract battery model parameter values from vector x
% Code source: https://github.com/Battery-Intelligence-Lab/multiscale-coupling/blob/dba3649f2eeb4228d2849e4d4e748c91ed6ad7a4/ObjFunc.m#L5C1-L17C48
tf          = sim_time(floor(C_rate));       % Total simulation time 
Cp          = x(1)*1e3;   % Specific heat
k           = x(2);       % Thermal conductivity
h           = x(3);       % Heat transfer coefficient   
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
dirpath = fullfile(dirpath, strrep(sprintf('C_rate_%.02f_T_amb_%.02f_h_%.01f_k_%.01f_%s', C_rate, T0, h, k, sim_type), '.', '_'));
[status, msg] = mkdir(dirpath);
if or(status ~= 1, contains(msg,'already')) 
    error(msg)
end

% Send battery model to COMSOL server
% Code source: https://github.com/Battery-Intelligence-Lab/multiscale-coupling/blob/dba3649f2eeb4228d2849e4d4e748c91ed6ad7a4/ObjFunc.m#L20
out = Battery_model_energetic(sim_type,C_rate,dirpath,tf,T0,Cp,k,h,ku,kappa,alfa_ka,i0ref,Ei0,Ds,DS,Vinit);

% Finish timing
toc