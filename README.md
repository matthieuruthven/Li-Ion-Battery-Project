# Framework to generate synthetic thermal images of a pouch cell

This repository contains the MATLAB code to generate synthetic thermal images of a pouch cell. The framework uses the computational model developed by [Lin *et al.* (2022)](https://www.nature.com/articles/s44172-022-00005-8) and implemented using [COMSOL Multiphysics](https://www.comsol.com/).

## Requirements

### Software

- COMSOL Multiphysics 6.1 with the Battery Design Module and LiveLink for MATLAB
- MATLAB R2023b

### Code

- Energetic_Generate_Dataset.m
- Battery_model_energetic.m
- txt_to_array.m

## Instructions to generate synthetic thermal images

1. Clone the cvi2 branch of this repository.

```
git clone https://github.com/matthieuruthven/Li-Ion-Battery-Project.git
```

2. Start LiveLink for MATLAB. This will open MATLAB.

3. Run the `Energetic_Generate_Dataset` function in MATLAB.

```
Energetic_Generate_Dataset(sim_type, sim_time, C_rate, k, h, dirpath, sq_wave_period, img_syn_int)
```

Where:

- sim_type (string) is *Charge*, *Discharge* or *Square* to generate images of the pouch cell charging, discharging or square wave loading
- sim_time (integer) is the duration of the simulation in seconds
- C_rate (float) is the C-rate at which the pouch cell is (dis)charged
- k (float) is the thermal conductivity of the pouch cell
- h (float) is the heat transfer coefficient of the pouch cell
- dirpath (string) is the path to the folder where the images should be saved
- sq_wave_period (integer): the duration of the square wave in seconds (only relevant when the sim_type argument is *Square*)
- img_syn_int (integer) is the time interval between the synthesis of consecutive images

For example:

```
Energetic_Generate_Dataset('Charge', 2000, 4, 1.1, 12, 'C:\Users\matthieu.ruthven\Documents\Li-Ion-Battery-Project', 4000, 10)
```

would simulate the charging of a pouch cell at a C-rate = 4 for 2000 seconds and generate sythetic images of the surface of the pouch cell every 10 seconds and save these images in a folder with path 'C:\Users\matthieu.ruthven\Documents\Li-Ion-Battery-Project\C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Charge'