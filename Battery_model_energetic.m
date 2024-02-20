function out = Battery_model_energetic(sim_type,sq_wave_period,img_syn_int,C_rate,dirpath,tf,T0,Cpel,keff,h,ku,kappa,alfa_ka,i0ref,Ei0,Ds,DS,U0)
%
% Battery_model_energetic.m
% 
% Function to implement pouch cell computational model developed by Lin et
% al. (2022)
% 
% Function adapted version of Battery_model.m file (available here:
% https://github.com/Battery-Intelligence-Lab/multiscale-coupling)
% 
% Function adapted by Matthieu Ruthven (matthieu.ruthven@uni.lu)
% Last modification made on 20th February 2024

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

localpath = pwd();
model.modelPath(localpath);

model.label('Battery model.mph');

model.param.set('L_sep', '20 [um]*N', 'Separator thickness');
model.param.set('L_pos', '70 [um]*N', 'Positive electrode thickness');
model.param.set('L_neg', '40 [um]*N', 'Negative electrode thickness');
model.param.set('L_pos_cc', '25 [um]*N/2', 'Positive current collector thickness');
model.param.set('L_neg_cc', '25 [um]*N/2', 'Negative current collector thickness');
model.param.set('W_cell', '150 [mm]', 'Cell width');
model.param.set('H_cell', '200 [mm]', 'Cell height');
model.param.set('L_tab', '0.25 [mm]', 'Tab thickness');
model.param.set('H_tab', '25 [mm]', 'Negative tab height');
model.param.set('W_tab', '48 [mm]', 'Negative tab width');
model.param.set('H_bar', '6 [in]', 'Copper bar height');
model.param.set('L_bar', '3.25 [mm]', 'Copper bar thickness');
model.param.set('Itab', '10 [mm]', 'Tab distance');
model.param.set('rp_pos', '5 [um]', 'Positive electrode particle radius');
model.param.set('rp_neg', '5 [um]', 'Negative electrode particle radius');
model.param.set('epss_pos', '0.5', 'Possitive electrode solid volume');
model.param.set('epsl_pos', '0.5', 'Positive electrode porosity');
model.param.set('epss_neg', '0.5', 'Negative electrode solid volume');
model.param.set('epsl_neg', '0.5', 'Negative electrode porosity');
model.param.set('Q_cell', '20[Ah]', 'Capacity of simulated cell geometry');
model.param.set('I_1C', 'Q_cell/1[h]/W_tab/L_bar', 'Cell 1C current for this geometry');
model.param.set('C_rate', num2str(C_rate), 'C-rate during simulation');
model.param.set('Iapp', 'I_1C*C_rate', 'Applied current');
model.param.set('cl0', '1.2[mol/l]', 'initial Electrolyte concentration');
model.param.set('T0', [num2str(T0) '[degC]'], 'Ambient temperature');
model.param.set('U0', [num2str(U0) '[V]'], 'Initial cell voltage');
if strcmp(sim_type, 'Discharge')
    model.param.set('x0', '0.199', 'Initial Negative Electrode SOC');
    model.param.set('y0', '0.725', 'Initial Positive Electrode SOC');
else
    model.param.set('x0', '0.92', 'Initial Negative Electrode SOC');
    model.param.set('y0', '0.01', 'Initial Positive Electrode SOC');
end
model.param.set('N', '42', 'Number of cell layers');
model.param.set('L_cell', 'L_sep+L_pos+L_neg+L_neg_cc+L_pos_cc', 'Cell thickness');
model.param.set('Rho', '2300 [kg/m^3]', 'Electrode density');
model.param.set('Cp', [num2str(Cpel) '[J/kg/K]'], 'Electrode specific heat');
model.param.set('k', [num2str(keff) '[W/m/K]'], 'Thermal conductivity of electrode');
model.param.set('ha', [num2str(h) '[W/m^2/K]'], 'Heat transfer coefficient');
model.param.set('ku', [num2str(ku) '[V]'], 'Slope of OCV curve');
model.param.set('kappa', [num2str(kappa) '[S/m]'], 'Reference ionic conductivity');
model.param.set('alfa_ka', [num2str(alfa_ka) '[S/m/K]'], 'Temperature coefficient for kappa');
model.param.set('i0ref', [num2str(i0ref) '[A/m^2]'], 'Exch current dens at T0');
model.param.set('Ei0', [num2str(Ei0) '[kJ/mol]'], 'Temp dependence of i0');
model.param.set('Ds', [num2str(Ds) '[m^2/s]'], 'Diffusion coefficient');
model.param.set('sq_wave_period', [num2str(sq_wave_period) '[s]'], 'Period of square wave');
model.param.set('img_syn_int', [num2str(img_syn_int) '[s]'], 'Time interval between synthesis of consecutive images');
model.param.set('sim_time', [num2str(tf) '[s]'], 'Duration of simulation');

% Save battery model parameters
model.param.saveFile(fullfile(dirpath, 'Model_Parameter_Values.txt'));

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.result.table.create('evl3', 'Table');
model.result.table.create('tbl1', 'Table');
model.result.table.create('evl2', 'Table');

model.func.create('wv1', 'Wave');
model.func('wv1').set('type', 'square');
model.func('wv1').set('freq', '2*pi*.01');
if strcmp(sim_type, 'Charge')
    model.func('wv1').set('amplitude', -1);
end
model.func('wv1').set('period', num2str(sq_wave_period));

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh.create('mesh2');

model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'W_cell' 'H_cell'});
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').set('distance', {'L_pos_cc' 'L_pos_cc+L_pos' 'L_pos_cc+L_pos+L_sep' 'L_pos_cc+L_pos+L_sep+L_neg' 'L_pos_cc+L_pos+L_sep+L_neg+L_neg_cc'});
model.component('comp1').geom('geom1').feature('ext1').set('scale', [1 1; 1 1; 1 1; 1 1; 1 1]);
model.component('comp1').geom('geom1').feature('ext1').set('displ', [0 0; 0 0; 0 0; 0 0; 0 0]);
model.component('comp1').geom('geom1').feature('ext1').set('twist', [0 0 0 0 0]);
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').label('Block - Positive Tab');
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'0' 'Itab' 'H_cell'});
model.component('comp1').geom('geom1').feature('blk1').set('size', {'L_tab' 'W_tab' 'H_tab'});
model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').label('Block - Negative Tab');
model.component('comp1').geom('geom1').feature('blk2').set('pos', {'L_pos_cc+L_pos+L_sep+L_neg+L_neg_cc-L_tab' 'W_cell-Itab-W_tab' 'H_cell'});
model.component('comp1').geom('geom1').feature('blk2').set('size', {'L_tab' 'W_tab' 'H_tab'});
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').active(false);
model.component('comp1').geom('geom1').feature('arr1').set('fullsize', {'1' '1' 'Ncc-1'});
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'0' '0' 'L_neg_cc+2*L_neg+2*L_sep+2*L_pos+L_pos_cc'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'blk1' 'blk2' 'ext1'});
model.component('comp1').geom('geom1').create('blk4', 'Block');
model.component('comp1').geom('geom1').feature('blk4').active(false);
model.component('comp1').geom('geom1').feature('blk4').label('Block - Negative CC');
model.component('comp1').geom('geom1').feature('blk4').set('pos', {'0' '0' '(Ncc-1)*L_neg_cc+2*(Ncc-1)*L_neg+2*(Ncc-1)*L_sep+2*(Ncc-1)*L_pos+(Ncc-1)*L_pos_cc'});
model.component('comp1').geom('geom1').feature('blk4').set('size', {'W_cell/2' 'H_cell' 'L_neg_cc'});
model.component('comp1').geom('geom1').create('blk3', 'Block');
model.component('comp1').geom('geom1').feature('blk3').active(false);
model.component('comp1').geom('geom1').feature('blk3').label('Block - Negative Tab 1');
model.component('comp1').geom('geom1').feature('blk3').set('pos', {'0.04875' '-H_neg_tab' '(Ncc-1)*L_neg_cc+2*(Ncc-1)*L_neg+2*(Ncc-1)*L_sep+2*(Ncc-1)*L_pos+(Ncc-1)*L_pos_cc'});
model.component('comp1').geom('geom1').feature('blk3').set('size', {'W_neg_tab/2' 'H_neg_tab' 'L_neg_cc'});
model.component('comp1').geom('geom1').create('wp2', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp2').active(false);
model.component('comp1').geom('geom1').feature('wp2').set('quickz', 'H_cell');
model.component('comp1').geom('geom1').feature('wp2').set('unite', true);
model.component('comp1').geom('geom1').feature('wp2').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('r1').set('pos', {'0' 'Itab'});
model.component('comp1').geom('geom1').feature('wp2').geom.feature('r1').set('size', {'L_pos_cc/2' 'W_tab'});
model.component('comp1').geom('geom1').feature('wp2').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('r2').set('pos', {'L_pos_cc/2+L_pos+L_sep+L_neg' 'W_cell-Itab-W_tab'});
model.component('comp1').geom('geom1').feature('wp2').geom.feature('r2').set('size', {'L_neg_cc/2' 'W_tab'});
model.component('comp1').geom('geom1').create('wp3', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp3').set('quickz', 'H_cell+H_tab');
model.component('comp1').geom('geom1').feature('wp3').set('unite', true);
model.component('comp1').geom('geom1').feature('wp3').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('r1').set('pos', {'L_tab-L_bar' 'Itab'});
model.component('comp1').geom('geom1').feature('wp3').geom.feature('r1').set('size', {'L_bar' 'W_tab'});
model.component('comp1').geom('geom1').feature('wp3').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('r2').set('pos', {'L_pos_cc+L_pos+L_sep+L_neg+L_neg_cc-L_bar' 'W_cell-Itab-W_tab'});
model.component('comp1').geom('geom1').feature('wp3').geom.feature('r2').set('size', {'L_bar' 'W_tab'});
model.component('comp1').geom('geom1').create('ext2', 'Extrude');
model.component('comp1').geom('geom1').feature('ext2').setIndex('distance', 'H_bar', 0);
model.component('comp1').geom('geom1').feature('ext2').selection('input').set({'wp3'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').selection.create('sel10', 'Explicit');
model.component('comp1').selection('sel10').set([3]);
model.component('comp1').selection.create('sel2', 'Explicit');
model.component('comp1').selection('sel2').set([2]);
model.component('comp1').selection.create('sel3', 'Explicit');
model.component('comp1').selection('sel3').set([4]);
model.component('comp1').selection.create('sel7', 'Explicit');
model.component('comp1').selection('sel7').set([6]);
model.component('comp1').selection.create('sel6', 'Explicit');
model.component('comp1').selection('sel6').set([7]);
model.component('comp1').selection.create('sel5', 'Explicit');
model.component('comp1').selection('sel5').set([8]);
model.component('comp1').selection.create('sel11', 'Explicit');
model.component('comp1').selection('sel11').set([9]);
model.component('comp1').selection.create('sel12', 'Explicit');
model.component('comp1').selection('sel12').set([1 5]);
model.component('comp1').selection.create('sel9', 'Explicit');
model.component('comp1').selection('sel9').geom('geom1', 2);
model.component('comp1').selection('sel9').set([4]);
model.component('comp1').selection.create('sel8', 'Explicit');
model.component('comp1').selection('sel8').geom('geom1', 2);
model.component('comp1').selection('sel8').set([26]);
model.component('comp1').selection.create('uni2', 'Union');
model.component('comp1').selection.create('uni1', 'Union');
model.component('comp1').selection.create('uni3', 'Union');
model.component('comp1').selection('sel10').label('Positive Tab');
model.component('comp1').selection('sel2').label('Positive Current Collector');
model.component('comp1').selection('sel3').label('Positive Electrode');
model.component('comp1').selection('sel7').label('Separator');
model.component('comp1').selection('sel6').label('Negative Electrode');
model.component('comp1').selection('sel5').label('Negative Current Collector');
model.component('comp1').selection('sel11').label('Negative Tab ');
model.component('comp1').selection('sel12').label('Bar');
model.component('comp1').selection('sel9').label('Positive Tab End');
model.component('comp1').selection('sel8').label('Negative Tab End');
model.component('comp1').selection('uni2').label('Positive Current Collector and Tab');
model.component('comp1').selection('uni2').set('input', {'sel2' 'sel10'});
model.component('comp1').selection('uni1').label('Negative Current Collector and Tab');
model.component('comp1').selection('uni1').set('input', {'sel5' 'sel11'});
model.component('comp1').selection('uni3').label('Metal Foil Domains');
model.component('comp1').selection('uni3').set('input', {'uni1' 'uni2'});

model.component('comp1').variable.create('var2');
model.component('comp1').variable('var2').set('cs_neg0', 'x0*mat3.elpot.cEeqref', 'Initial negative electrode concentration');
model.component('comp1').variable('var2').set('cs_pos0', 'y0*mat4.elpot.cEeqref', 'Initial positive electrode concentration');
model.component('comp1').variable('var2').set('Eq_neg0', 'mat3.elpot.Eeq_int3(x0)+mat3.elpot.dEeqdT_int1(x0)*(T-298[K])', 'Initial positive equilibrium potential');
model.component('comp1').variable('var2').set('Eq_pos0', 'mat4.elpot.Eeq_int1(y0)+mat4.elpot.dEeqdT_int1(y0)*(T-298[K])', 'Initial negative equilibrium potential');
model.component('comp1').variable('var2').set('sigmal', 'mat5.ionc.sigmal_iso', 'Effective Ionic conductivity');

model.component('comp1').view('view3').tag('view5');
model.component('comp1').view('view4').tag('view6');
model.view.create('view7', 2);

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material.create('mat4', 'Common');
model.component('comp1').material.create('mat5', 'Common');
model.component('comp1').material.create('mat6', 'Common');
model.component('comp1').material('mat1').selection.set([3 9]);
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').propertyGroup.create('Murnaghan', 'Murnaghan');
model.component('comp1').material('mat1').propertyGroup.create('Lame', ['Lam' native2unicode(hex2dec({'00' 'e9'}), 'unicode') ' parameters']);
model.component('comp1').material('mat2').selection.set([2 8]);
model.component('comp1').material('mat2').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat2').propertyGroup.create('linzRes', 'Linearized resistivity');
model.component('comp1').material('mat3').selection.named('sel6');
model.component('comp1').material('mat3').propertyGroup.create('ElectrodePotential', 'Equilibrium potential');
model.component('comp1').material('mat3').propertyGroup.create('OperationalSOC', 'Operational electrode state-of-charge');
model.component('comp1').material('mat4').selection.named('sel3');
model.component('comp1').material('mat4').propertyGroup.create('ElectrodePotential', 'Equilibrium potential');
model.component('comp1').material('mat4').propertyGroup.create('OperationalSOC', 'Operational electrode state-of-charge');
model.component('comp1').material('mat5').selection.named('sel7');
model.component('comp1').material('mat5').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat5').propertyGroup.create('ElectrolyteConductivity', 'Electrolyte conductivity');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func.create('int1', 'Interpolation');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func.create('an1', 'Analytic');
model.component('comp1').material('mat5').propertyGroup.create('SpeciesProperties', 'Species properties');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').func.create('an1', 'Analytic');
model.component('comp1').material('mat5').propertyGroup.create('ElectrolyteSaltConcentration', 'Electrolyte salt concentration');
model.component('comp1').material('mat6').selection.set([1 5]);
model.component('comp1').material('mat6').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat6').propertyGroup.create('linzRes', 'Linearized resistivity');

model.component('comp1').physics.create('liion', 'LithiumIonBatteryMPH', 'geom1');
model.component('comp1').physics('liion').create('init2', 'init', 3);
model.component('comp1').physics('liion').feature('init2').selection.set([1 2 3 4]);
model.component('comp1').physics('liion').create('pce1', 'PorousElectrode', 3);
model.component('comp1').physics('liion').feature('pce1').selection.named('sel6');
model.component('comp1').physics('liion').create('pce2', 'PorousElectrode', 3);
model.component('comp1').physics('liion').feature('pce2').selection.named('sel3');
model.component('comp1').physics('liion').create('sep1', 'Separator', 3);
model.component('comp1').physics('liion').feature('sep1').selection.named('sel7');
model.component('comp1').physics('liion').create('ece1', 'Electrode', 3);
model.component('comp1').physics('liion').feature('ece1').selection.set([1 3 5 9]);
model.component('comp1').physics('liion').create('ece2', 'Electrode', 3);
model.component('comp1').physics('liion').feature('ece2').selection.set([2 8]);
model.component('comp1').physics('liion').create('ec1', 'ElectrodeCurrent', 2);
model.component('comp1').physics('liion').feature('ec1').selection.named('sel9');
model.component('comp1').physics('liion').create('egnd1', 'ElectricGround', 2);
model.component('comp1').physics('liion').feature('egnd1').selection.named('sel8');
model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').physics('ht').create('solid6', 'SolidHeatTransferModel', 3);
model.component('comp1').physics('ht').feature('solid6').selection.set([4 7]);
model.component('comp1').physics('ht').create('solid4', 'SolidHeatTransferModel', 3);
model.component('comp1').physics('ht').feature('solid4').selection.named('sel7');
model.component('comp1').physics('ht').create('solid5', 'SolidHeatTransferModel', 3);
model.component('comp1').physics('ht').feature('solid5').selection.set([2 8]);
model.component('comp1').physics('ht').create('hs1', 'HeatSource', 3);
model.component('comp1').physics('ht').feature('hs1').selection.all;
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 2);
model.component('comp1').physics('ht').feature('hf1').selection.set([4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 26 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50]);

model.component('comp1').multiphysics.create('ech1', 'ElectrochemicalHeating', -1);

model.component('comp1').mesh('mesh2').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh2').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh2').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh2').feature('ftet1').selection.set([4 6 7]);

model.result.table('evl3').label('Evaluation 3D');
model.result.table('evl3').comments('Interactive 3D values');
model.result.table('tbl1').comments('Surface Average 1 (T)');
model.result.table('evl2').label('Evaluation 2D');
model.result.table('evl2').comments('Interactive 2D values');

model.component('comp1').view('view1').set('scenelight', false);
model.component('comp1').view('view2').axis.set('xmin', -0.09386149793863297);
model.component('comp1').view('view2').axis.set('xmax', 0.30907174944877625);
model.component('comp1').view('view2').axis.set('ymin', -0.00894896686077118);
model.component('comp1').view('view2').axis.set('ymax', 0.204908549785614);
model.component('comp1').view('view5').label('View 5');
model.component('comp1').view('view5').axis.set('xmin', -0.033154115080833435);
model.component('comp1').view('view5').axis.set('xmax', 0.08002246171236038);
model.component('comp1').view('view5').axis.set('ymin', 5.716224550269544E-4);
model.component('comp1').view('view5').axis.set('ymax', 0.05878385901451111);
model.component('comp1').view('view6').label('View 6');
model.component('comp1').view('view6').axis.set('xmin', -40.121212005615234);
model.component('comp1').view('view6').axis.set('xmax', 69.3652572631836);
model.component('comp1').view('view6').axis.set('ymin', 66.00369262695312);
model.component('comp1').view('view6').axis.set('ymax', 143.69981384277344);
model.view('view7').axis.set('xmin', -207.71853637695312);
model.view('view7').axis.set('xmax', 57.718536376953125);
model.view('view7').axis.set('ymin', -5);
model.view('view7').axis.set('ymax', 205);

model.component('comp1').material('mat1').label('Tab');
model.component('comp1').material('mat1').set('family', 'aluminum');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '505[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'300[W/(m*K)]' '0' '0' '0' '300[W/(m*K)]' '0' '0' '0' '300[W/(m*K)]'});
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'3e7[S/m]' '0' '0' '0' '3e7[S/m]' '0' '0' '0' '3e7[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'23e-6[1/K]' '0' '0' '0' '23e-6[1/K]' '0' '0' '0' '23e-6[1/K]'});
model.component('comp1').material('mat1').propertyGroup('def').set('density', '5830[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('youngsmodulus', '70e9[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('poissonsratio', '0.33');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('l', '-2.5e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('m', '-3.3e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Murnaghan').set('n', '-3.5e11[Pa]');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '');
model.component('comp1').material('mat1').propertyGroup('Lame').set('lambLame', '5.1e10[Pa]');
model.component('comp1').material('mat1').propertyGroup('Lame').set('muLame', '2.6e10[Pa]');
model.component('comp1').material('mat2').label('CC');
model.component('comp1').material('mat2').set('family', 'copper');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'N^2*1e7[S/m]' '0' '0' '0' '1e7[S/m]' '0' '0' '0' '1e7[S/m]'});
model.component('comp1').material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'17e-6[1/K]' '0' '0' '0' '17e-6[1/K]' '0' '0' '0' '17e-6[1/K]'});
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', '505[J/kg/K]');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat2').propertyGroup('def').set('density', '5830[kg/m^3]');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'282[W/m/K]' '0' '0' '0' '282*0.18[W/m/K]' '0' '0' '0' '282[W/m/K]'});
model.component('comp1').material('mat2').propertyGroup('Enu').set('youngsmodulus', '110e9[Pa]');
model.component('comp1').material('mat2').propertyGroup('Enu').set('poissonsratio', '0.35');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('rho0', '');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('alpha', '');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('Tref', '');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('rho0', '');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('alpha', '');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('Tref', '');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('rho0', '1.72e-8[ohm*m]');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('alpha', '0.0039[1/K]');
model.component('comp1').material('mat2').propertyGroup('linzRes').set('Tref', '298[K]');
model.component('comp1').material('mat2').propertyGroup('linzRes').addInput('temperature');
model.component('comp1').material('mat3').label('Graphite Electrode, LixC6 MCMB (Negative, Li-ion Battery)');
model.component('comp1').material('mat3').comments(['Eeq for fully lithiated at 0.98\nEeq for delithiated at 0\n\nReferences\nV. Srinivasan, and J. Newman, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Design and Optimization of a Natural Graphite/Iron Phosphate Lithium Ion Cell,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 151, p. 1530, 2004.\n\nK. Kumaresan, G. Sikha, and R. E. White, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Thermal Model for a Li-Ion Cell,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 155, p. A164, 2008.\n\nR. E. Gerver and J. P. Meyers, "Three-Dimensional Modeling of Electrochemical Performance and Heat Generation of Lithium-Ion Betteries in Tabbed Planar Configurations", J. Electrochemical Soc., vol. 158, p. A835, 2011\n\nK. E. Thomas, and J. Newman, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Heats of mixing and of entropy in porous insertion electrodes,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Power Sources., vol. 119-121, p. 844, 2003.']);
model.component('comp1').material('mat3').set('groups', {'electrodes' 'Electrodes'});
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', {'N^2*50[S/m]' '0' '0' '0' '50[S/m]' '0' '0' '0' '50[S/m]'});
model.component('comp1').material('mat3').propertyGroup('def').set('diffusion', {'Ds' '0' '0' '0' 'Ds' '0' '0' '0' 'Ds'});
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'10' '0' '0' '0' '10' '0' '0' '0' '10'});
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '1000');
model.component('comp1').material('mat3').propertyGroup('def').set('density', '3200');
model.component('comp1').material('mat3').propertyGroup('def').set('T_ref', '318[K]');
model.component('comp1').material('mat3').propertyGroup('def').set('T2', 'min(393.15,max(T,223.15))');
model.component('comp1').material('mat3').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('Eeq', '');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('dEeqdT', '');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('cEeqref', '');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('Eeq', '');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('dEeqdT', '');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('cEeqref', '');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('Eeq', '-ku*(soc-x0)-0.02*wv1(root.t[1/s])');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('dEeqdT', [num2str(-DS) '[mV/K]']);
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('cEeqref', '29612[mol/m^3]');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').set('soc', 'c/cEeqref');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').addInput('concentration');
model.component('comp1').material('mat3').propertyGroup('ElectrodePotential').addInput('temperature');
model.component('comp1').material('mat3').propertyGroup('OperationalSOC').set('socmax', '');
model.component('comp1').material('mat3').propertyGroup('OperationalSOC').set('socmin', '');
model.component('comp1').material('mat3').propertyGroup('OperationalSOC').set('socmax', '0.98');
model.component('comp1').material('mat3').propertyGroup('OperationalSOC').set('socmin', '0');
model.component('comp1').material('mat4').label('LFP Electrode, LiFePO4 (Positive, Li-ion Battery)');
model.component('comp1').material('mat4').comments(['vs Li/Li+, T=25 C\nEeq for fully lithiated at 0.90\nEeq for delithiated at 0.01\n\nReferences\nU. S. Kasavajjula, C. Wang, and P. E. Arce, "Discharge Model for LiFePO4 Accounting for the Solid Solution Range", J. Electrochemical Soc., vol. 155, p. A866, 2008\n\nR. E. Gerver and J. P. Meyers, "Three-Dimensional Modeling of Electrochemical Performance and Heat Generation of Lithium-Ion Betteries in Tabbed Planar Configurations", J. Electrochemical Soc., vol. 158, p. A835, 2011']);
model.component('comp1').material('mat4').set('groups', {'electrodes' 'Electrodes'});
model.component('comp1').material('mat4').propertyGroup('def').set('diffusion', {'Ds' '0' '0' '0' 'Ds' '0' '0' '0' 'Ds'});
model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'N^2*50[S/m]' '0' '0' '0' '50[S/m]' '0' '0' '0' '50[S/m]'});
model.component('comp1').material('mat4').propertyGroup('def').set('thermalconductivity', {'10' '0' '0' '0' '10' '0' '0' '0' '10'});
model.component('comp1').material('mat4').propertyGroup('def').set('heatcapacity', '1000');
model.component('comp1').material('mat4').propertyGroup('def').set('density', '3200');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('Eeq', '');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('dEeqdT', '');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('cEeqref', '');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('Eeq', '');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('dEeqdT', '');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('cEeqref', '');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('Eeq', 'U0-ku*(soc-y0)+0.02*wv1(root.t[1/s])'); 
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('dEeqdT', [num2str(DS) '[mV/K]']);
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('cEeqref', '16921[mol/m^3]');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').set('soc', 'c/cEeqref');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').addInput('concentration');
model.component('comp1').material('mat4').propertyGroup('ElectrodePotential').addInput('temperature');
model.component('comp1').material('mat4').propertyGroup('OperationalSOC').set('socmax', '');
model.component('comp1').material('mat4').propertyGroup('OperationalSOC').set('socmin', '');
model.component('comp1').material('mat4').propertyGroup('OperationalSOC').set('socmax', '0.90');
model.component('comp1').material('mat4').propertyGroup('OperationalSOC').set('socmin', '0.01');
model.component('comp1').material('mat5').label('LiPF6 in 1:2 EC:DMC and p(VdF-HFP) (Polymer electrolyte, Li-ion Battery)');
model.component('comp1').material('mat5').comments(['T=25 C\n\n1:2 EC/DMC by volume\n\nReference\nM. Doyle, J. Newman, A.S. Gozdz, C.N. Schmutz, and J.M. Tarascon, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 143, p. 1890, 1996.']);
model.component('comp1').material('mat5').set('groups', {'electrolytes' 'Electrlytes'});
model.component('comp1').material('mat5').propertyGroup('def').func('an1').set('funcname', 'D_iso');
model.component('comp1').material('mat5').propertyGroup('def').func('an1').set('expr', '1.47e3*exp(1.33*c)*exp(-1.69e3/T)*exp(-5.63e2/T*c)*1E-10');
model.component('comp1').material('mat5').propertyGroup('def').func('an1').set('args', {'T' 'c'});
model.component('comp1').material('mat5').propertyGroup('def').func('an1').set('argunit', 'K,1');
model.component('comp1').material('mat5').propertyGroup('def').func('an1').set('fununit', 'm^2/s');
model.component('comp1').material('mat5').propertyGroup('def').func('an1').set('plotargs', {'T' '263.15' '323.15'; 'c' '0' '3'});
model.component('comp1').material('mat5').propertyGroup('def').set('diffusion', {'N^2*7.5e-11' '0' '0' '0' '7.5e-11' '0' '0' '0' '7.5e-11'});
model.component('comp1').material('mat5').propertyGroup('def').set('c_ref', '1000[mol/m^3]');
model.component('comp1').material('mat5').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat5').propertyGroup('def').addInput('concentration');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('int1').set('funcname', 'sigmal_int1');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('int1').set('table', {'0' '0.0108';  ...
'0.2000' '0.1259';  ...
'0.4000' '0.2055';  ...
'0.6000' '0.2553';  ...
'0.8000' '0.2810';  ...
'1.0000' '0.2873';  ...
'1.2000' '0.2788';  ...
'1.4000' '0.2595';  ...
'1.6000' '0.2331';  ...
'1.8000' '0.2027';  ...
'2.0000' '0.1710';  ...
'2.200' '0.1403';  ...
'2.4000' '0.1123';  ...
'2.6000' '0.0885';  ...
'2.8000' '0.0696';  ...
'3.0000' '0.0563'});
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('an1').set('funcname', 'sigmal_int2');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('an1').set('expr', '0.1*0.798*(1+(T-228))*c*(1-1.22*sqrt(c)+0.509*(1-4e-3*exp(1000/T))*c)/(1+c^4*(3.79e-3*exp(1000/T)))');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('an1').set('args', {'T' 'c'});
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('an1').set('argunit', 'K,1');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('an1').set('fununit', 'S/m');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').func('an1').set('plotargs', {'T' '263.15' '323.15'; 'c' '0' '3'});
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').set('sigmal', {'kappa/N^2+(T-T_ref)*alfa_ka/N^2' '0' '0' '0' 'kappa/N^2+(T-T_ref)*alfa_ka/N^2' '0' '0' '0' 'kappa/N^2+(T-T_ref)*alfa_ka/N^2'});
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').set('c_ref', '1000[mol/m^3]');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').set('T_ref', 'T0');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').addInput('temperature');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteConductivity').addInput('concentration');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').func('an1').set('funcname', 'transpNm_int1');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').func('an1').set('expr', '-7.91+2.45e-1*c+5.28e-2*T+6.98e-1*c^2-1.08e-2*c*T-8.21e-5*T^2+7.43e-4*c^3-2.22e-3*c^2*T+3.07e-5*c*T^2');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').func('an1').set('args', {'T' 'c'});
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').func('an1').set('argunit', 'K,1');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').func('an1').set('plotargs', {'T' '263.15' '323.15'; 'c' '0' '3'});
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('transpNum', '');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('fcl', '');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('transpNum', '');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('fcl', '');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('transpNum', '0.363');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('fcl', '0');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').set('c_ref', '1000[mol/m^3]');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').addInput('temperature');
model.component('comp1').material('mat5').propertyGroup('SpeciesProperties').addInput('concentration');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteSaltConcentration').identifier('cElsalt');
model.component('comp1').material('mat5').propertyGroup('ElectrolyteSaltConcentration').set('cElsalt', '1000[mol/m^3]');
model.component('comp1').material('mat6').label('Copper bar');
model.component('comp1').material('mat6').set('family', 'copper');
model.component('comp1').material('mat6').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat6').propertyGroup('def').set('electricconductivity', {'5.998e7[S/m]' '0' '0' '0' '5.998e7[S/m]' '0' '0' '0' '5.998e7[S/m]'});
model.component('comp1').material('mat6').propertyGroup('def').set('thermalexpansioncoefficient', {'17e-6[1/K]' '0' '0' '0' '17e-6[1/K]' '0' '0' '0' '17e-6[1/K]'});
model.component('comp1').material('mat6').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.component('comp1').material('mat6').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat6').propertyGroup('def').set('density', '8960[kg/m^3]');
model.component('comp1').material('mat6').propertyGroup('def').set('thermalconductivity', {'400[W/(m*K)]' '0' '0' '0' '400[W/(m*K)]' '0' '0' '0' '400[W/(m*K)]'});
model.component('comp1').material('mat6').propertyGroup('Enu').set('youngsmodulus', '110e9[Pa]');
model.component('comp1').material('mat6').propertyGroup('Enu').set('poissonsratio', '0.35');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('rho0', '');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('alpha', '');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('Tref', '');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('rho0', '');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('alpha', '');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('Tref', '');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('rho0', '1.72e-8[ohm*m]');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('alpha', '0.0039[1/K]');
model.component('comp1').material('mat6').propertyGroup('linzRes').set('Tref', '298[K]');
model.component('comp1').material('mat6').propertyGroup('linzRes').addInput('temperature');

model.component('comp1').physics('liion').feature('ice1').set('ElectrolyteMaterial', 'mat5');
model.component('comp1').physics('liion').feature('ice1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('init1').set('cl', 'cl0');
model.component('comp1').physics('liion').feature('init2').set('cl', 'cl0');
model.component('comp1').physics('liion').feature('init2').set('phis', 'U0');
model.component('comp1').physics('liion').feature('pce1').set('ElectrolyteMaterial', 'mat5');
model.component('comp1').physics('liion').feature('pce1').set('sigmal', {'sigmal*N^2'; '0'; '0'; '0'; 'sigmal'; '0'; '0'; '0'; 'sigmal'});
model.component('comp1').physics('liion').feature('pce1').set('ElectrodeMaterial', 'mat3');
model.component('comp1').physics('liion').feature('pce1').set('sigma', {'sigmas*N^2'; '0'; '0'; '0'; 'sigmas'; '0'; '0'; '0'; 'sigmas'});
model.component('comp1').physics('liion').feature('pce1').set('epsl', 'epsl_neg');
model.component('comp1').physics('liion').feature('pce1').set('epss', 'epss_neg');
model.component('comp1').physics('liion').feature('pce1').set('IonicCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('pce1').set('fl', 'epsl_neg^brug');
model.component('comp1').physics('liion').feature('pce1').set('ElectricCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('pce1').set('fs', 'epss_neg^brug');
model.component('comp1').physics('liion').feature('pce1').set('DiffusionCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('pce1').set('fDl', 'epsl_neg^brug');
model.component('comp1').physics('liion').feature('pce1').set('Rfilm', '0.01[ohm*m^2]');
model.component('comp1').physics('liion').feature('pce1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('pce1').label('Porous Electrode - Negative');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('csinit', 'cs_neg0');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('Ds', '3.9e-14[m^2/s]');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('ParticleMaterial', 'mat3');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('cEeqref', '31507[mol/m^3]');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('socmax', 0.98);
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('rp', 'rp_neg');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('pce1').feature('per1').set('MaterialOption', 'mat3');
model.component('comp1').physics('liion').feature('pce1').feature('per1').set('ElectrodeKinetics', 'LinButlerVolmer');
model.component('comp1').physics('liion').feature('pce1').feature('per1').set('i0', 'i0ref*exp(-Ei0/R_const*(1/T-1/T0))');
% model.component('comp1').physics('liion').feature('pce1').feature('per1').set('k_a', '6e-11*exp(-Ei0/R_const*(1/T-1/T0))');
% model.component('comp1').physics('liion').feature('pce1').feature('per1').set('k_c', '6e-11*exp(-Ei0/R_const*(1/T-1/T0))');
model.component('comp1').physics('liion').feature('pce1').feature('per1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('pce2').set('ElectrolyteMaterial', 'mat5');
model.component('comp1').physics('liion').feature('pce2').set('sigmal', {'sigmal*N^2'; '0'; '0'; '0'; 'sigmal'; '0'; '0'; '0'; 'sigmal'});
model.component('comp1').physics('liion').feature('pce2').set('ElectrodeMaterial', 'mat4');
model.component('comp1').physics('liion').feature('pce2').set('sigma', {'sigmas*N^2'; '0'; '0'; '0'; 'sigmas'; '0'; '0'; '0'; 'sigmas'});
model.component('comp1').physics('liion').feature('pce2').set('epsl', 'epsl_pos');
model.component('comp1').physics('liion').feature('pce2').set('epss', 'epss_pos');
model.component('comp1').physics('liion').feature('pce2').set('IonicCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('pce2').set('fl', 'epsl_pos^brug');
model.component('comp1').physics('liion').feature('pce2').set('ElectricCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('pce2').set('fs', 'epss_pos^brug');
model.component('comp1').physics('liion').feature('pce2').set('DiffusionCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('pce2').set('fDl', 'epsl_pos^brug');
model.component('comp1').physics('liion').feature('pce2').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('pce2').label('Porous Electrode - Positive');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('csinit', 'cs_pos0');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('Ds', '3.2e-13[m^2/s]');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('ParticleMaterial', 'mat4');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('cEeqref', '21190[mol/m^3]');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('socmin', 0.01);
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('socmax', 0.9);
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('rp', 'rp_pos');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('pce2').feature('per1').set('MaterialOption', 'mat4');
model.component('comp1').physics('liion').feature('pce2').feature('per1').set('Eeq', 'mat4.elpot.Eeq_int1(y0)-mat3.elpot.Eeq_int1(x0)');
model.component('comp1').physics('liion').feature('pce2').feature('per1').set('ElectrodeKinetics', 'LinButlerVolmer');
model.component('comp1').physics('liion').feature('pce2').feature('per1').set('i0', 'i0ref*exp(-Ei0/R_const*(1/T-1/T0))');
% model.component('comp1').physics('liion').feature('pce2').feature('per1').set('k_a', '5e-10*exp(-Ei0/R_const*(1/T-1/T0))');
% model.component('comp1').physics('liion').feature('pce2').feature('per1').set('k_c', '5e-10*exp(-Ei0/R_const*(1/T-1/T0))');
model.component('comp1').physics('liion').feature('pce2').feature('per1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('sep1').set('ElectrolyteMaterial', 'mat5');
model.component('comp1').physics('liion').feature('sep1').set('sigmal', {'sigmal*N^2'; '0'; '0'; '0'; 'sigmal'; '0'; '0'; '0'; 'sigmal'});
model.component('comp1').physics('liion').feature('sep1').set('IonicCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('sep1').set('fl', '(0.4)^brug');
model.component('comp1').physics('liion').feature('sep1').set('DiffusionCorrModel', 'NoCorr');
model.component('comp1').physics('liion').feature('sep1').set('fDl', '(0.4)^brug');
model.component('comp1').physics('liion').feature('sep1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('ece1').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('ece2').set('sigma', {'sigmas_cc*N^2'; '0'; '0'; '0'; 'sigmas_cc'; '0'; '0'; '0'; 'sigmas_cc'});
model.component('comp1').physics('liion').feature('ece2').set('minput_temperature', 'T0');
model.component('comp1').physics('liion').feature('ec1').set('ElectronicCurrentType', 'AverageCurrentDensity');
model.component('comp1').physics('liion').feature('ec1').set('Its', '-I_1C*C_rate');
model.component('comp1').physics('liion').feature('ec1').set('Ias', 'wv1(root.t[1/s])*Iapp');
model.component('comp1').physics('liion').feature('ec1').set('phis0init', 0);
model.component('comp1').physics('ht').prop('PhysicalModelProperty').set('Tref', 'T0');
model.component('comp1').physics('ht').feature('solid1').set('minput_strainreferencetemperature', 'T0');
model.component('comp1').physics('ht').feature('init1').set('Tinit', 'T0');
model.component('comp1').physics('ht').feature('solid6').set('k', {'k'; '0'; '0'; '0'; 'k'; '0'; '0'; '0'; 'k'});
model.component('comp1').physics('ht').feature('solid6').set('rho', 'Rho');
model.component('comp1').physics('ht').feature('solid6').set('Cp', 'Cp');
model.component('comp1').physics('ht').feature('solid6').set('minput_strainreferencetemperature', 'T0');
model.component('comp1').physics('ht').feature('solid6').label('Electrode');
model.component('comp1').physics('ht').feature('solid4').set('k', {'k'; '0'; '0'; '0'; 'k'; '0'; '0'; '0'; 'k'});
model.component('comp1').physics('ht').feature('solid4').set('rho', 'Rho');
model.component('comp1').physics('ht').feature('solid4').set('Cp', 'Cp');
model.component('comp1').physics('ht').feature('solid4').label('Separator');
model.component('comp1').physics('ht').feature('solid5').set('rho', 'Rho');
model.component('comp1').physics('ht').feature('solid5').set('Cp', 'Cp');
model.component('comp1').physics('ht').feature('solid5').label('Current collector');
model.component('comp1').physics('ht').feature('hs1').set('Q0', 'liion.Qh');
model.component('comp1').physics('ht').feature('hs1').set('materialType', 'from_mat');
model.component('comp1').physics('ht').feature('hf1').set('HeatFluxType', 'ConvectiveHeatFlux');
model.component('comp1').physics('ht').feature('hf1').set('h', 'ha');
model.component('comp1').physics('ht').feature('hf1').set('Text', 'T0');
model.component('comp1').physics('ht').feature('hf1').set('materialType', 'from_mat');

model.component('comp1').multiphysics('ech1').active(false);

model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh2').feature('size').set('hauto', 8);
model.component('comp1').mesh('mesh2').feature('ftet1').set('xscale', 2);
model.component('comp1').mesh('mesh2').feature('ftet2').set('xscale', 6);
model.component('comp1').mesh('mesh2').run;

model.component('comp1').physics('liion').feature('ice1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('pce1').set('sigmal_mat', 'userdef');
model.component('comp1').physics('liion').feature('pce1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('socmin_mat', 'userdef');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('socmax_mat', 'userdef');
model.component('comp1').physics('liion').feature('pce1').feature('pin1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('pce1').feature('per1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('pce2').set('sigmal_mat', 'userdef');
model.component('comp1').physics('liion').feature('pce2').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('socmin_mat', 'userdef');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('socmax_mat', 'userdef');
model.component('comp1').physics('liion').feature('pce2').feature('pin1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('pce2').feature('per1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('sep1').set('sigmal_mat', 'userdef');
model.component('comp1').physics('liion').feature('sep1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('ece1').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('liion').feature('ece2').set('minput_temperature_src', 'root.comp1.T');
model.component('comp1').physics('ht').feature('solid1').set('minput_strainreferencetemperature_src', 'userdef');
model.component('comp1').physics('ht').feature('solid6').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid6').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid6').set('Cp_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid6').set('minput_strainreferencetemperature_src', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('rho_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid4').set('Cp_mat', 'userdef');

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('se1', 'Segregated');
model.sol('sol1').feature('t1').create('d1', 'Direct');
model.sol('sol1').feature('t1').create('d2', 'Direct');
model.sol('sol1').feature('t1').create('d3', 'Direct');
model.sol('sol1').feature('t1').create('i1', 'Iterative');
model.sol('sol1').feature('t1').create('i2', 'Iterative');
model.sol('sol1').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ss3', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ss4', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').create('ll1', 'LowerLimit');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ssDef');
model.sol('sol1').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.dataset.create('surf1', 'Surface');
model.result.dataset.create('cpt1', 'CutPoint3D');
model.result.dataset.create('cpt2', 'CutPoint3D');
model.result.dataset.create('cpt3', 'CutPoint3D');
model.result.dataset.create('cpt4', 'CutPoint3D');
model.result.dataset.create('cln1', 'CutLine3D');
model.result.dataset('surf1').selection.set([6]);
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical.create('pev2', 'EvalPoint');
model.result.numerical.create('pev3', 'EvalPoint');
model.result.numerical.create('pev4', 'EvalPoint');
model.result.numerical.create('av1', 'AvSurface');
model.result.numerical('gev1').set('probetag', 'none');
model.result.numerical('pev1').set('probetag', 'none');
model.result.numerical('pev2').set('probetag', 'none');
model.result.numerical('pev3').set('probetag', 'none');
model.result.numerical('pev4').set('probetag', 'none');
model.result.numerical('av1').selection.set([6]);
model.result.numerical('av1').set('probetag', 'none');
model.result.create('pg16', 'PlotGroup1D');
model.result.create('pg24', 'PlotGroup1D');
model.result.create('pg17', 'PlotGroup1D');
model.result.create('pg18', 'PlotGroup3D');
model.result.create('pg19', 'PlotGroup3D');
model.result.create('pg22', 'PlotGroup3D');
model.result.create('pg20', 'PlotGroup3D');
model.result.create('pg21', 'PlotGroup3D');
model.result.create('pg26', 'PlotGroup3D');
model.result.create('pg23', 'PlotGroup3D');
model.result.create('pg25', 'PlotGroup2D');
model.result('pg16').create('glob1', 'Global');
model.result('pg24').create('ptgr1', 'PointGraph');
model.result('pg24').create('tblp1', 'Table');
model.result('pg24').create('ptgr2', 'PointGraph');
model.result('pg17').create('glob1', 'Global');
model.result('pg18').create('mslc1', 'Multislice');
model.result('pg18').create('arwv1', 'ArrowVolume');
model.result('pg19').create('arwv1', 'ArrowVolume');
model.result('pg19').create('vol1', 'Volume');
model.result('pg19').create('con1', 'Contour');
model.result('pg19').feature('arwv1').create('col1', 'Color');
model.result('pg22').create('mslc1', 'Multislice');
model.result('pg22').create('arwv1', 'ArrowVolume');
model.result('pg22').create('arwv2', 'ArrowVolume');
model.result('pg20').create('mslc1', 'Multislice');
model.result('pg21').create('arwv1', 'ArrowVolume');
model.result('pg21').feature('arwv1').create('col1', 'Color');
model.result('pg26').create('vol1', 'Volume');
model.result('pg23').create('vol1', 'Volume');
model.result('pg23').create('con1', 'Contour');
model.result('pg25').create('surf1', 'Surface');
model.result('pg25').create('con1', 'Contour');
model.result.export.create('tbl1', 'Table');
model.result.export.create('plot1', 'Plot');
model.result.export.create('anim1', 'Animation');

model.study('std1').feature('time').set('tlist', sprintf('range(0,%d,%d)', img_syn_int, tf));

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {sprintf('range(0,%d,%d)', img_syn_int, tf) '2.5[s]'});
model.sol('sol1').feature('v1').feature('comp1_cl').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_cl').set('scaleval', 1000);
model.sol('sol1').feature('v1').feature('comp1_liion_pce1_cs').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_liion_pce1_cs').set('scaleval', 10000);
model.sol('sol1').feature('v1').feature('comp1_liion_pce2_cs').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_liion_pce2_cs').set('scaleval', 10000);
model.sol('sol1').feature('v1').feature('comp1_phil').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_phil').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_phis').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_phis').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_liion_ec1_phis0').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_liion_ec1_phis0').set('scaleval', 1);
model.sol('sol1').feature('t1').set('tlist', sprintf('range(0,%d,%d)', img_syn_int, tf));
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('eventout', true);
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').feature('se1').set('segstabacc', 'segaacc');
model.sol('sol1').feature('t1').feature('se1').set('segaaccdim', 5);
model.sol('sol1').feature('t1').feature('se1').set('segaaccmix', 0.9);
model.sol('sol1').feature('t1').feature('se1').set('segaaccdelay', 1);
model.sol('sol1').feature('t1').feature('se1').feature('ss1').label('Heat transfer T');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('segvar', {'comp1_T'});
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('linsolver', 'd1');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('subdamp', 0.8);
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('subjtech', 'once');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_phil' 'comp1_cl' 'comp1_phis' 'comp1_liion_ec1_phis0'});
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('linsolver', 'd2');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subdtech', 'auto');
model.sol('sol1').feature('t1').feature('se1').feature('ss3').set('segvar', {'comp1_liion_pce1_cs'});
model.sol('sol1').feature('t1').feature('se1').feature('ss3').set('linsolver', 'd3');
model.sol('sol1').feature('t1').feature('se1').feature('ss4').set('segvar', {'comp1_liion_pce2_cs'});
model.sol('sol1').feature('t1').feature('se1').feature('ss4').set('linsolver', 'd3');
model.sol('sol1').feature('t1').feature('se1').feature('ll1').set('lowerlimit', 'comp1.T 0');
model.sol('sol1').feature('t1').feature('d1').label('PARDISO (ht)');
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('d2').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d2').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('d3').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d3').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('i1').label('Algebraic Multigrid (ht)');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('so1').set('relax', 0.9);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('i2').label('Geometric Multigrid (ht)');
model.sol('sol1').feature('t1').feature('i2').set('rhob', 20);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

%% Define hot/cold spot locations
model.result.dataset('surf1').label('Front surface');
model.result.dataset('cpt1').label('Hot spot');
model.result.dataset('cpt1').set('pointx', 0);
model.result.dataset('cpt1').set('pointy', 77.791);
model.result.dataset('cpt1').set('pointz', 189.58);
model.result.dataset('cpt2').label('Cold spot');
model.result.dataset('cpt2').set('pointx', 0);
model.result.dataset('cpt2').set('pointy', 'W_cell-20[mm]');
model.result.dataset('cpt2').set('pointz', 'H_cell-175[mm]');
model.result.dataset('cln1').label('Ty');
model.result.dataset('cln1').set('genpoints', {'0' '0' 'H_cell-25[mm]'; '0' 'W_cell' 'H_cell-25[mm]'});
model.result.numerical('gev1').label('Cell voltage');
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').set('expr', {'liion.phis0_ec1'});
model.result.numerical('gev1').set('unit', {'V'});
model.result.numerical('gev1').set('descr', {'Electric potential on boundary'});
model.result.numerical('pev1').label('Hot');
model.result.numerical('pev1').set('data', 'cpt1');
model.result.numerical('pev1').set('table', 'tbl1');
model.result.numerical('pev1').set('unit', {'degC'});
model.result.numerical('pev2').label('Cold');
model.result.numerical('pev2').set('data', 'cpt2');
model.result.numerical('pev2').set('table', 'tbl1');
model.result.numerical('pev2').set('unit', {'degC'});
model.result.numerical('av1').label('Average surface temperature');
model.result.numerical('av1').set('table', 'tbl1');
model.result.numerical('av1').set('unit', {'degC'});
model.result.numerical('gev1').setResult;
model.result.numerical('pev1').appendResult;
model.result.numerical('pev2').appendResult;
model.result.numerical('av1').appendResult;
model.result('pg16').label('Boundary Electrode Potential Vs. Ground (liion)');
model.result('pg16').set('xlabel', 'Capacity (Ah)');
model.result('pg16').set('xlabelactive', true);
model.result('pg16').set('ylabel', 'Electric potential on boundary (V)');
model.result('pg16').set('ylabelactive', false);
model.result('pg16').feature('glob1').set('expr', {'liion.phis0_ec1'});
model.result('pg16').feature('glob1').set('unit', {'V'});
model.result('pg16').feature('glob1').set('descr', {'Electric potential on boundary'});
model.result('pg24').label('Temperature response');
model.result('pg24').set('xlabel', 'Time (s)');
model.result('pg24').set('ylabel', 'Temperature (degC)');
model.result('pg24').set('legendpos', 'lowerright');
model.result('pg24').set('xlabelactive', false);
model.result('pg24').set('ylabelactive', false);
model.result('pg24').feature('ptgr1').set('data', 'cpt1');
model.result('pg24').feature('ptgr1').set('looplevelinput', {'all'});
model.result('pg24').feature('ptgr1').set('expr', 'T');
model.result('pg24').feature('ptgr1').set('unit', 'degC');
model.result('pg24').feature('ptgr1').set('descr', 'Temperature');
model.result('pg24').feature('ptgr1').set('linecolor', 'blue');
model.result('pg24').feature('tblp1').set('table', 'tbl1');
model.result('pg24').feature('tblp1').set('linecolor', 'blue');
model.result('pg24').feature('ptgr2').set('data', 'cpt2');
model.result('pg24').feature('ptgr2').set('looplevelinput', {'all'});
model.result('pg24').feature('ptgr2').set('expr', 'T');
model.result('pg24').feature('ptgr2').set('unit', 'degC');
model.result('pg24').feature('ptgr2').set('descr', 'Temperature');
model.result('pg24').feature('ptgr2').set('linecolor', 'blue');
model.result('pg17').label('Average Electrode State-of-Charge (liion)');
model.result('pg17').set('xlabel', 'Time (s)');
model.result('pg17').set('xlabelactive', false);
model.result('pg17').feature('glob1').set('expr', {'liion.soc_average_pce1' 'liion.soc_average_pce2'});
model.result('pg17').feature('glob1').set('unit', {'1' '1'});
model.result('pg17').feature('glob1').set('descr', {'Average SOC, Porous Electrode - Negative' 'Average SOC, Porous Electrode - Positive'});
model.result('pg18').label('Electrolyte Potential (liion)');
model.result('pg18').set('looplevel', [1]);
model.result('pg18').feature('mslc1').set('resolution', 'normal');
model.result('pg18').feature('arwv1').set('arrowbase', 'center');
model.result('pg18').feature('arwv1').set('scale', 2.6025711135785717E10);
model.result('pg18').feature('arwv1').set('color', 'black');
model.result('pg18').feature('arwv1').set('scaleactive', false);
model.result('pg19').label('Electrolyte Current Density (liion)');
model.result('pg19').set('looplevel', [1]);
model.result('pg19').feature('arwv1').active(false);
model.result('pg19').feature('arwv1').set('arrowbase', 'center');
model.result('pg19').feature('arwv1').set('scale', 5.668630456310266E-6);
model.result('pg19').feature('arwv1').set('color', 'black');
model.result('pg19').feature('arwv1').set('scaleactive', false);
model.result('pg19').feature('arwv1').feature('col1').set('expr', 'root.comp1.liion.IlMag');
model.result('pg19').feature('arwv1').feature('col1').set('unit', 'A/m^2');
model.result('pg19').feature('arwv1').feature('col1').set('descr', 'Electrolyte current density magnitude');
model.result('pg19').feature('vol1').set('expr', 'liion.IlMag');
model.result('pg19').feature('vol1').set('unit', 'A/m^2');
model.result('pg19').feature('vol1').set('descr', 'Electrolyte current density magnitude');
model.result('pg19').feature('vol1').set('colortable', 'Thermal');
model.result('pg19').feature('vol1').set('resolution', 'normal');
model.result('pg19').feature('con1').set('expr', 'liion.IlMag');
model.result('pg19').feature('con1').set('unit', 'A/m^2');
model.result('pg19').feature('con1').set('descr', 'Electrolyte current density magnitude');
model.result('pg19').feature('con1').set('resolution', 'normal');
model.result('pg22').label('Electrolyte Salt Concentration (liion)');
model.result('pg22').set('looplevel', [1]);
model.result('pg22').feature('mslc1').set('expr', 'cl');
model.result('pg22').feature('mslc1').set('unit', 'mol/m^3');
model.result('pg22').feature('mslc1').set('descr', 'Electrolyte salt concentration');
model.result('pg22').feature('mslc1').set('resolution', 'normal');
model.result('pg22').feature('arwv1').set('expr', {'liion.Nposx' 'liion.Nposy' 'liion.Nposz'});
model.result('pg22').feature('arwv1').set('descr', 'Positive ion flux');
model.result('pg22').feature('arwv1').set('scale', 4.1350822164665405E15);
model.result('pg22').feature('arwv1').set('scaleactive', false);
model.result('pg22').feature('arwv2').set('expr', {'liion.Nnegx' 'liion.Nnegy' 'liion.Nnegz'});
model.result('pg22').feature('arwv2').set('descr', 'Negative ion flux');
model.result('pg22').feature('arwv2').set('scale', 2.358230255019835E15);
model.result('pg22').feature('arwv2').set('color', 'black');
model.result('pg22').feature('arwv2').set('inheritplot', 'arwv1');
model.result('pg22').feature('arwv2').set('inheritcolor', false);
model.result('pg22').feature('arwv2').set('scaleactive', false);
model.result('pg20').label('Electrode Potential Vs. Ground (liion)');
model.result('pg20').set('looplevel', [1]);
model.result('pg20').feature('mslc1').set('expr', 'phis');
model.result('pg20').feature('mslc1').set('descr', 'Electric potential');
model.result('pg20').feature('mslc1').set('resolution', 'normal');
model.result('pg21').label('Electrode Current Density (liion)');
model.result('pg21').set('looplevel', [1]);
model.result('pg21').feature('arwv1').set('expr', {'liion.Isx' 'liion.Isy' 'liion.Isz'});
model.result('pg21').feature('arwv1').set('descr', 'Electrode current density vector');
model.result('pg21').feature('arwv1').set('arrowbase', 'center');
model.result('pg21').feature('arwv1').set('scale', 130.1545435962273);
model.result('pg21').feature('arwv1').set('scaleactive', false);
model.result('pg21').feature('arwv1').feature('col1').set('expr', 'root.comp1.liion.IsMag');
model.result('pg21').feature('arwv1').feature('col1').set('unit', 'A/m^2');
model.result('pg21').feature('arwv1').feature('col1').set('descr', 'Electrode current density magnitude');
model.result('pg26').label('SOC');
model.result('pg26').set('looplevel', [1]);
model.result('pg26').feature('vol1').set('expr', 'liion.sodloc_average');
model.result('pg26').feature('vol1').set('unit', '1');
model.result('pg26').feature('vol1').set('descr', 'Local electrode material state-of-discharge, particle average');
model.result('pg26').feature('vol1').set('resolution', 'normal');
model.result('pg23').label('Battery temperature');
model.result('pg23').set('looplevel', [1]);
model.result('pg23').set('showlegendsmaxmin', true);
model.result('pg23').set('showlegendsunit', true);
model.result('pg23').feature('vol1').set('expr', 'T');
model.result('pg23').feature('vol1').set('unit', 'degC');
model.result('pg23').feature('vol1').set('descr', 'Temperature');
model.result('pg23').feature('vol1').set('colortable', 'Thermal');
model.result('pg23').feature('vol1').set('resolution', 'normal');
model.result('pg23').feature('con1').set('expr', 'T');
model.result('pg23').feature('con1').set('unit', 'degC');
model.result('pg23').feature('con1').set('descr', 'Temperature');
model.result('pg23').feature('con1').set('resolution', 'normal');
model.result('pg25').label('Battery surface temperature');
model.result('pg25').feature('surf1').set('expr', 'T');
model.result('pg25').feature('surf1').set('unit', 'degC');
model.result('pg25').feature('surf1').set('descr', 'Temperature');
model.result('pg25').feature('surf1').set('rangecoloractive', true);
model.result('pg25').feature('surf1').set('rangecolormin', 24);
model.result('pg25').feature('surf1').set('rangecolormax', 34);
model.result('pg25').feature('surf1').set('colortable', 'HeatCamera');
model.result('pg25').feature('surf1').set('resolution', 'normal');
model.result('pg25').feature('con1').set('expr', 'T');
model.result('pg25').feature('con1').set('unit', 'degC');
model.result('pg25').feature('con1').set('descr', 'Temperature');
model.result('pg25').feature('con1').set('number', 15);
model.result('pg25').feature('con1').set('resolution', 'normal');

model.label('LinetalModel.mph');
model.label('LinetalModel.mph');

model.result('pg18').run;
model.result('pg19').run;
model.result('pg22').run;
model.result('pg20').run;
model.result('pg21').run;

model.sol('sol1').runAll;

model.result('pg16').run;
model.result('pg16').run;
model.result('pg24').run;
model.result('pg17').run;
model.result('pg18').run;
model.result('pg19').run;
model.result('pg22').run;
model.result('pg20').run;
model.result('pg21').run;
model.result('pg26').run;
model.result('pg23').run;
model.result('pg25').run;
model.result('pg16').run;
model.result('pg24').run;
model.result('pg25').run;
model.result('pg23').run;
model.result('pg26').run;
model.result('pg21').run;

model.label('Lin_et_al_Model_Original.mph');

model.result('pg21').run;

model.label('Lin_et_al_Model_Original.mph');

% Export parameters values
es = struct;
es.img_w_pix = 256;
es.img_h_pix = 256;
es.img_type = 'png';
es.video_type = 'gif';
if strcmp(es.video_type, 'avi')
    es.avi_quality = 1;
end
es.video_fps = 20;
es.colour_table = 'Rainbow';
% es.colour_table_min = 23;
% es.colour_table_max = 30;
es.font_size = 9;

% Save these values
writestruct(es, fullfile(dirpath, 'Export_Parameter_Values.json'))

model.result.dataset.create('surf2', 'Surface');
model.result.dataset('surf2').label('Positive CC');
model.result.dataset('surf2').selection.set([6 10 43]);
model.result.create('pg27', 'PlotGroup2D');
model.result('pg27').run;
model.result('pg27').label('Positive CC Surface T');
model.result('pg27').set('data', 'surf2');
model.result('pg27').set('titletype', 'custom');
model.result('pg27').set('typeintitle', false);
model.result('pg27').set('descriptionintitle', true);
model.result('pg27').set('unitintitle', false);
model.result('pg27').set('edges', false);
model.result('pg27').set('showlegends', true);
model.result('pg27').set('showlegendsmaxmin', true);
model.result('pg27').set('showlegendsunit', true);
model.result('pg27').create('surf1', 'Surface');
model.result('pg27').feature('surf1').set('expr', 'T');
model.result('pg27').feature('surf1').set('unit', 'degC');
model.result('pg27').feature('surf1').set('colortable', es.colour_table);
% model.result('pg27').feature('surf1').set('rangecoloractive', true);
% model.result('pg27').feature('surf1').set('rangecolormin', es.colour_table_min);
% model.result('pg27').feature('surf1').set('rangecolormax', es.colour_table_max);
model.result('pg27').run;

model.view('view10').set('showgrid', false);

% Export average state of charge (SOC) as TXT file
model.result.export.create('plot3', 'Plot');
model.result.export('plot3').label('Average SOC');
model.result.export('plot3').set('plotgroup', 'pg17');
model.result.export('plot3').set('filename', fullfile(dirpath, 'Average_Electrode_SOC.txt'));
model.result.export('plot3').set('header', false);
model.result.export('plot3').run;

% Read TXT file
SOC_vs_t = readmatrix(fullfile(dirpath, 'Average_Electrode_SOC.txt'));

% Extract negative electrode SOCs
neg_SOC_vs_t = SOC_vs_t(1:(size(SOC_vs_t, 1) / 2), :);

% Extract positive electrode SOCs
pos_SOC_vs_t = SOC_vs_t((size(SOC_vs_t, 1) / 2 + 1):end, :);

% Find simulation end point
if strcmp(sim_type, 'Charge')
    
    % Find simulation step ID at which SOC of negative electrode > 0
    neg_step_id = find(neg_SOC_vs_t(:, 2) > 0, 1, 'last');

    % Find simulation step ID at which SOC of positive electrode < 0.90
    pos_step_id = find(pos_SOC_vs_t(:, 2) < 0.90, 1, 'last');

    % Find minimum simulation step ID
    max_step_id = min(neg_step_id, pos_step_id);

elseif strcmp(sim_type, 'Discharge')

    % Find maximum time at which SOC of negative electrode < 0.98
    neg_step_id = find(neg_SOC_vs_t(:, 2) < 0.98, 1, 'last');

    % Find maximum time at which SOC of positive electrode > 0.01
    pos_step_id = find(pos_SOC_vs_t(:, 2) > 0.01, 1, 'last');

    % Find minimum simulation step ID
    max_step_id = min(neg_step_id, pos_step_id);

else

    max_step_id = size(neg_SOC_vs_t, 1);

end

% Export RGB video with labels
model.result.export.create('anim2', 'Animation');
model.result.export('anim2').set('fontsize', num2str(es.font_size));
model.result.export('anim2').set('colortheme', 'globaltheme');
model.result.export('anim2').set('customcolor', [1 1 1]);
model.result.export('anim2').set('background', 'transparent');
model.result.export('anim2').set('gltfincludelines', 'on');
model.result.export('anim2').set('title1d', 'on');
model.result.export('anim2').set('legend1d', 'on');
model.result.export('anim2').set('logo1d', 'on');
model.result.export('anim2').set('options1d', 'on');
model.result.export('anim2').set('title2d', 'on');
model.result.export('anim2').set('legend2d', 'on');
model.result.export('anim2').set('logo2d', 'off');
model.result.export('anim2').set('options2d', 'on');
model.result.export('anim2').set('title3d', 'on');
model.result.export('anim2').set('legend3d', 'on');
model.result.export('anim2').set('logo3d', 'on');
model.result.export('anim2').set('options3d', 'off');
model.result.export('anim2').set('axisorientation', 'on');
model.result.export('anim2').set('grid', 'on');
model.result.export('anim2').set('axes1d', 'on');
model.result.export('anim2').set('axes2d', 'off');
model.result.export('anim2').set('showgrid', 'on');
model.result.export('anim2').label('Positive CC Surface T');
model.result.export('anim2').set('plotgroup', 'pg27');
model.result.export('anim2').set('type', 'movie');
model.result.export('anim2').set('movietype', es.video_type);
if strcmp(es.video_type, 'avi')
    model.result.export('anim2').set('aviqual', es.avi_quality);
end
model.result.export('anim2').set(sprintf('%sfilename', es.video_type), fullfile(dirpath, sprintf('Surface_T_Labelled.%s', es.video_type)));
model.result.export('anim2').set('looplevelinput', 'manualindices');
model.result.export('anim2').set('looplevelindices', sprintf('range(1,1,%d)', max_step_id));
% model.result.export('anim2').set('looplevelindices', '1,101,201,301,401');
model.result.export('anim2').set('framesel', 'all');
model.result.export('anim2').set('width', es.img_w_pix);
model.result.export('anim2').set('height', es.img_h_pix);
model.result.export('anim2').set('fps', es.video_fps);
model.result.export('anim2').set('options2d', true);
model.result.export('anim2').run;

% Export greyscale video with labels
model.result('pg27').feature('surf1').set('colortable', 'GrayScale');
model.result('pg27').run;
model.result.export('anim2').set(sprintf('%sfilename', es.video_type), fullfile(dirpath, sprintf('Surface_T_Labelled_Grayscale.%s', es.video_type)));
model.result.export('anim2').run;

% Export greyscale video without labels
model.result('pg27').set('titletype', 'none');
model.result('pg27').run;
model.result.export('anim2').set('options2d', 'off');
model.result.export('anim2').set(sprintf('%sfilename', es.video_type), fullfile(dirpath, sprintf('Surface_T_Unlabelled_Grayscale.%s', es.video_type)));
model.result.export('anim2').run;

% Export RGB video without labels
model.result('pg27').feature('surf1').set('colortable', es.colour_table);
model.result('pg27').run;
model.result.export('anim2').set(sprintf('%sfilename', es.video_type), fullfile(dirpath, sprintf('Surface_T_Unlabelled.%s', es.video_type)));
model.result.export('anim2').run;

% Export RGB images
model.result.export('anim2').set('type', 'imageseq');
model.result.export('anim2').set('imagefilename', fullfile(dirpath, sprintf('Surface_T_.%s', es.img_type)));
model.result.export('anim2').run;

% Export greyscale images (3-channel)
model.result('pg27').feature('surf1').set('colortable', 'GrayScale');
model.result('pg27').run;
model.result.export('anim2').set('imagefilename', fullfile(dirpath, sprintf('Surface_T_Grayscale_.%s', es.img_type)));
model.result.export('anim2').run;

% Ensure zero padding is consistent in file name
zero_padding = numel(num2str(max_step_id));
assert(zero_padding < 5, 'Zero padding is greater than 4')
if zero_padding ~= 4
    for idx = 1:max_step_id
        if zero_padding == 3
            src_filepath = fullfile(dirpath, sprintf('Surface_T_%03d.%s', idx, es.img_type));
            other_src_filepath = fullfile(dirpath, sprintf('Surface_T_Grayscale_%03d.%s', idx, es.img_type));
        elseif zero_padding == 2
            src_filepath = fullfile(dirpath, sprintf('Surface_T_%02d.%s', idx, es.img_type));
            other_src_filepath = fullfile(dirpath, sprintf('Surface_T_Grayscale_%02d.%s', idx, es.img_type));
        else
            src_filepath = fullfile(dirpath, sprintf('Surface_T_%d.%s', idx, es.img_type));
            other_src_filepath = fullfile(dirpath, sprintf('Surface_T_Grayscale_%d.%s', idx, es.img_type));
        end
        tgt_filepath = fullfile(dirpath, sprintf('Surface_T_%04d.%s', idx, es.img_type));
        movefile(src_filepath, tgt_filepath)
        tgt_filepath = fullfile(dirpath, sprintf('Surface_T_Grayscale_%04d.%s', idx, es.img_type));
        movefile(other_src_filepath, tgt_filepath)
    end
end

% Convert greyscale images to 1-channel
% for idx = 1:max_t
%     tmp_filepath = fullfile(dirpath, sprintf('Surface_T_grayscale_%04d.%s', idx, es.img_type));
%     tmp_img = imread(tmp_filepath);
%     tmp_img = rgb2gray(tmp_img);
%     imwrite(tmp_img, tmp_filepath)
% end

% Create folder in which raw temperature distributions will be saved (an error will be raised if
% the folder already exists)
tmp_dirpath = fullfile(dirpath, 'Raw_Surface_T');
[status, msg] = mkdir(tmp_dirpath);
if or(status ~= 1, contains(msg,'already')) 
    error(msg)
end

% Export temperature distributions as TXT files
surface_t_array = [];
model.result.export.create('plot2', 'Plot');
for idx = 1:max_step_id
    model.result.dataset('surf2').selection.set([6]);
    model.result('pg27').setIndex('looplevel', idx, 0);
    model.result('pg27').run;
    model.result.export('plot2').label('Positive CC Surface T Values');
    model.result.export('plot2').set('plotgroup', 'pg27');
    model.result.export('plot2').set('filename', fullfile(tmp_dirpath, sprintf('Surface_T_Battery_%04d.txt', idx)));
    model.result.export('plot2').set('header', false);
    model.result.export('plot2').run;
    model.result.dataset('surf2').selection.set([10]);
    model.result('pg27').run;
    model.result.export('plot2').set('filename', fullfile(tmp_dirpath, sprintf('Surface_T_Neg_Tab_%04d.txt', idx)));
    model.result.export('plot2').run;
    model.result.dataset('surf2').selection.set([43]);
    model.result('pg27').run;
    model.result.export('plot2').set('filename', fullfile(tmp_dirpath, sprintf('Surface_T_Pos_Tab_%04d.txt', idx)));
    model.result.export('plot2').run;
    % Create single array of temperature values
    tmp_t_vals = txt_to_array(dirpath, idx);
    surface_t_array = cat(3, surface_t_array, tmp_t_vals);
end
% Save surface_t_array
save(fullfile(dirpath, 'Surface_T.mat'), 'surface_t_array')

out = model;