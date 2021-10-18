%%
clear all;clc;
%% add path
%addpath('/Users/xiaoyuwang/Desktop/Research code/VCBM_drug_model-main/TumourModelParticleFiltering')
%addpath('/Users/xiaoyuwang/Desktop/Research code/VCBM_drug_model-main/TumourModelParticleFiltering/Model')

addpath('C:/Users/johnw/Dropbox/ABM/TumourModelParticleFiltering')
addpath('C:/Users/johnw/Dropbox/ABM/TumourModelParticleFiltering/Model')

%%
page = 2; 
% set the time
max_time = 32;

% set true paramter
p0_true = 0.4;
psc_true = 1.05e-05;
dmax_true = 31;
gage_true = 162;

% set the observation value
Control_tumour_growth_data = zeros(1,max_time);
%wn = normrnd(0,1,[1,32]);
%init_vol = exp(linspace(1,8,max_time)) + abs(5*wn);

%%
startingvol = 200;
%for jj = 1:10
pf = clib.Model.SeedAndGrowToStartVolume(p0_true, psc_true, dmax_true, gage_true, page, startingvol); % starts simulation from 1 cell and grows until starting volumes
    for t = 1:max_time
       % 
       pf.SimulateOneDay(1);
       Control_tumour_growth_data(1,t) = pf.TumourVolume(); %calculates tumour volume of one particle's simulation
    end
%end

%%
figure
hold on
plot(linspace(1,max_time,max_time),Control_tumour_growth_data);

%%
save('Control_tumour_growth_data.mat','Control_tumour_growth_data')