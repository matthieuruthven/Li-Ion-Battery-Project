# Framework to generate synthetic thermal images of a pouch cell

This repository contains the MATLAB code to generate synthetic thermal images of a pouch cell. The framework uses the computational model developed by [Lin *et al.* (2022)](https://www.nature.com/articles/s44172-022-00005-8) and implemented using [COMSOL Multiphysics](https://www.comsol.com/).

## Requirements

### Software

- COMSOL Multiphysics 6.1 with the Battery Design Module and LiveLink for MATLAB
- MATLAB R2023b with the Image Processing Toolbox

### Code

- Energetic_Generate_Dataset.m
- Battery_model_energetic.m
- txt_to_array.m
- Anomaly_Img_to_Mask_Conversion.m

## Instructions to generate synthetic thermal images

1. Clone the cvi2 branch of this repository.

```
git clone https://github.com/matthieuruthven/Li-Ion-Battery-Project.git
```

2. Start LiveLink for MATLAB. This will open MATLAB.

3. Run the `Energetic_Generate_Dataset` function in MATLAB.

```
Energetic_Generate_Dataset(sim_type, sim_time, batt_type, anomaly, C_rate, k_in, h_in, T0_in, dirpath, sq_wave_period, img_syn_int, SOC, def_or_cus)
```

Where:

- **sim_type** (string) is *Charge*, *Discharge* or *Square* to generate images of the surface temperature of the pouch cell while it is charging, discharging or cycled using a square wave
- **sim_time** (integer) is the duration of the simulation in seconds
- **batt_type** (string) is either *Energetic* or *Lin* to indicate whether simulations using the Energetic or original pouch cell should be performed
- **anomaly** (string) is either *None*, *Circle*, *Stripe_H* or *Stripe_W* to indicate the spatial anomaly in the battery heat transfer coefficient
- **C_rate** (float) is the C-rate at which the pouch cell is (dis)charged
- **k_in** (float) is the thermal conductivity, *k*, of the pouch cell in W/m/K
- **h_in** (float) is the heat transfer coefficient, *h*, of the pouch cell in W/m<sup>2</sup>/K
- **T0_in** (float) is the initial temperature, *T<sub>0</sub>*, of the pouch cell in &deg;C
- **dirpath** (string) is the path to the folder where the images should be saved
- **sq_wave_period** (integer) is the period of the square wave in seconds (NB this argument is only relevant when the **sim_type** argument is *Square* and its value is overridden when the **sim_type** argument is *Charge* or *Discharge*)
- **img_syn_int** (integer) is the time interval in seconds between the synthesis of consecutive images
- **SOC** (float) is the initial state of charge (SOC) of the pouch cell (0 $\leq$ SOC $\leq$ 1)
- **def_or_cus** (string) is either *custom* or *default* to indicate if custom values of *k*, *h* and *T<sub>0</sub>* should be used or the default values (NB in the latter case, the values of the **k_in**, **h_in** and **T0_in** arguments will be overriden)

For example:

```
Energetic_Generate_Dataset('Square', 2500, 'Lin', 'None', 4, 1, 1, 1, 'C:\Users\matthieu.ruthven\Documents\Li-Ion-Battery-Project', 100, 1, 0.3, 'default')
```

would simulate the square wave cycling of the Lin *et al.* (2022) pouch cell at a C-rate of 4 for 2500 seconds and then generate sythetic images of the surface temperature of the pouch cell every second and save these images in a folder with path *C:\Users\matthieu.ruthven\Documents\Li-Ion-Battery-Project\Lin_C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Square_SOC_0_3*. The initial SOC of the pouch cell would be 30% and the default values of *k*, *h* and *T<sub>0</sub>* would be used in the simulation. The *h* of the pouch cell would not include any spatial anomaly.

## Summary of generated data

Running the `Energetic_Generate_Dataset` function will generate the following data:

```bash
Lin_C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Square_SOC_0_3/
┣ Raw_Surface_T/
┃  ┣ Surface_T_Battery_0001.txt      # Surface temperature in degrees Centigrade (at time point 1) at each node in battery mesh
┃  ┣ Surface_T_Neg_Tab_0001.txt      # Surface temperature in degrees Centigrade (at time point 1) at each node in negative tab mesh
┃  ┣ Surface_T_Pos_Tab_0001.txt      # Surface temperature in degrees Centigrade (at time point 1) at each node in positive tab mesh
┃  ┣ Surface_T_Battery_0002.txt      # Surface temperature in degrees Centigrade (at time point 2) at each node in battery mesh
┃  ┣ Surface_T_Neg_Tab_0002.txt      # Surface temperature in degrees Centigrade (at time point 2) at each node in negative tab mesh
┃  ┣ Surface_T_Pos_Tab_0002.txt      # Surface temperature in degrees Centigrade (at time point 2) at each node in positive tab mesh
┃  ┗ ...
┣ Average_Electrode_SOC.txt          # Average state of charge in negative and positive electrode at each time point
┣ Model_Parameter_Values.txt         # Values of parameters of computational model of battery, and parameters of exported images
┣ Surface_T.mat                      # Array of all surface temperatures in TXT files in Raw_Surface_T folder
┣ Surface_T_0001.png                 # Surface temperature of battery at time point 1
┣ Surface_T_0002.png                 # Surface temperature of battery at time point 2
┣ ...
┣ Surface_T_Labelled.gif             # Surface temperature distribution of battery at each time point, including labels
┗ Surface_T_Unlabelled.gif           # Surface temperature distribution of battery at each time point, not including labels
```

If **anomaly** is *Circle*, *Stripe_H* or *Stripe_W*, the function will also generate the following data:
```bash
Lin_C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Square_SOC_0_3/
┣ Heat_Transfer_Coefficient_1.png          # Distribution of heat transfer coefficient values in battery
┣ Heat_Transfer_Coefficient_2.png          # Distribution of heat transfer coefficient values in battery
┗ Heat_Transfer_Coefficient_Anomaly.png    # Binary mask showing region where heat transfer coefficient is anomalous
```

## Key limitation of binary mask conversion code

Two key assumptions made by the code to create a binary mask showing the region where *h* is anomalous are:

1. The largest region in the image (in terms of the number of pixels with the same value) is the background
2. The second largest region is the part of the battery where the transfer heat coefficient is **NOT** anomalous

The second assumption is not met when the anomalous region is particularly large. Consequently, the binary mask created is not of the anomalous region. This is a key limitation of the code.