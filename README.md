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
git clone -b cvi2 https://github.com/matthieuruthven/multiscale-coupling.git
```

2. Start LiveLink for MATLAB. This will open MATLAB.

3. Run the `Energetic_Generate_Dataset` function in MATLAB.

```
Energetic_Generate_Dataset(sim_type, C_rate, k, h, dirpath)
```

Where:

- sim_type (string) is either 'Charge' or 'Discharge' to generate images of the pouch cell charging or discharging
- C_rate (float) is the C-rate at which the pouch cell is (dis)charged
- k (float) is the thermal conductivity of the pouch cell
- h (float) is the heat transfer coefficient of the pouch cell
- dirpath (string) is the path to the folder where the images should be saved

For example:

```
Energetic_Generate_Dataset('Charge', 4, 1.1, 12, 'C:\Users\matthieu.ruthven\Documents\multiscale-coupling')
```

would generate sythetic images of the surface of a pouch cell charged at a C-rate=4 and save these images in a folder with path 'C:\Users\matthieu.ruthven\Documents\multiscale-coupling\C_rate_4_00_T_amb_23_85_h_12_0_k_1_1_Charge'