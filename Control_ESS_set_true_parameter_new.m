%%
clear all;clc;
%% add path
%addpath('/Users/xiaoyuwang/Desktop/Research code/VCBM_drug_model-main/TumourModelParticleFiltering')
%addpath('/Users/xiaoyuwang/Desktop/Research code/VCBM_drug_model-main/TumourModelParticleFiltering/Model')

addpath('C:\Users\jenner2\OneDrive - Queensland University of Technology\Documents\PF\PF')
addpath('C:\Users\jenner2\OneDrive - Queensland University of Technology\Documents\PF\PF\Model')
%%
load('Control_tumour_growth_data.mat');
Control_tumour_growth = Control_tumour_growth_data;
%% Load parameter dist mean, variance and minimum
% 

%dist means mu for the parameters
p0_mu = 0.3;
psc_mu = 1e-5;
dmax_mu = 30;
gage_mu = 160;

%dist variance for the parameters
p0_sigma = 1;
psc_sigma = 1;
dmax_sigma = 1;
gage_sigma = 1;

%parameter minimum vals
p0_min = 0;
psc_min = 0;
dmax_min = 0;
gage_min = 0;

page = 2; 

%%
% set the time
max_time = 32;

% set true paramter
p0_true = 0.1;
psc_true = 1.05e-05;
dmax_true = 31;
gage_true = 162;

%% generate prior distribution
seed = 42;
rng(seed,'twister') % rng(seed) specifies the seed for the random number generator 

% set the number of particles to simulate for 
% set turning parameters 
n = 20000;
threshold = 0.005;
a = (3*0.96-1)/(2*0.96);
h2 = 1-a^2;

% need to use transformation to get the initial value for theta_vec_normal,
% fisrt need to have the prior, the by use the values of prior we can
% campute the theta_vec_normal via transformation

p0_vec_pd = makedist('Normal', 'mu', 0.0, 'sigma', p0_sigma); %makes a normal dist with mean zero, variance sigma
p0_vec = zeros(max_time,n);
p0_vec_normal = zeros(max_time,n);
p0_vec_normal(1,:) = random(p0_vec_pd, n, 1);
p0_vec(1, :) = exp(log(p0_mu - p0_min) + p0_vec_normal(1,:)) + p0_min;

psc_mu_pd = makedist('Normal', 'mu', 0.0, 'sigma', psc_sigma);
psc_vec = zeros(max_time,n);
psc_vec_normal = zeros(max_time,n);
psc_vec_normal(1,:) = random(psc_mu_pd, n, 1);
psc_vec(1, :) = exp(log(psc_mu - psc_min) + psc_vec_normal(1,:)) + psc_min;

dmax_mu_pd = makedist('Normal', 'mu', 0.0, 'sigma', dmax_sigma);
dmax_vec = zeros(max_time,n);
dmax_vec_normal = zeros(max_time,n);
dmax_vec_normal(1,:) = random(dmax_mu_pd, n, 1);
dmax_vec(1, :) = exp(log(dmax_mu - dmax_min) + dmax_vec_normal(1,:)) + dmax_min;

gage_mu_pd = makedist('Normal', 'mu', 0.0, 'sigma', gage_sigma);
gage_vec = zeros(max_time,n);
gage_vec_normal = zeros(max_time,n);
gage_vec_normal(1,:) = random(gage_mu_pd, n, 1);
gage_vec(1, :) = floor(exp(log(gage_mu - gage_min) + gage_vec_normal(1,:))) + gage_min;

weight = zeros(max_time,n);
weight(1,:) = ones(1,n)/n;

%%
figure
subplot(2,2,1)
[f,xi] = ksdensity(p0_vec(1,:));
plot(xi,f,'b--');
hold on 
plot(p0_true,0,'rx');
subplot(2,2,2)
[f,xi] = ksdensity(psc_vec(1,:));
plot(xi,f,'b--');
hold on 
plot(psc_true,0,'rx');
subplot(2,2,3)
[f,xi] = ksdensity(dmax_vec(1,:));
plot(xi,f,'b--');
hold on 
plot(dmax_true,0,'rx');
subplot(2,2,4)
[f,xi] = ksdensity(gage_vec(1,:));
plot(xi,f,'b--');
hold on 
plot(gage_true,0,'rx');
%% Initialise particle filtering
%grow tumour from 1 cancer cell for n particles
startVolume = Control_tumour_growth(1,1); %volume to grow initial tumours to
PFTumourVolume = zeros(max_time,n);

for i = 1:n
    pf = clib.Model.SeedAndGrowToStartVolume(p0_vec(1,i), psc_vec(1,i), dmax_vec(1,i), gage_vec(1,i), page, startVolume);
    PFTumourVolume(1,i) = pf.TumourVolume(); %calculates tumour volume of one particle's simulation
    
    particle{i} = clib.Model.CreateNewParticle(p0_vec(1,i), psc_vec(1,i), dmax_vec(1,i), gage_vec(1,i), 2, pf); %stores particles start tumour pointer

    %create matrix of residuals of each particle to each mouse
%     for j = 1:n_mice
%         res_mat(1, j, i) =(pf.TumourVolume()- nanmean(Control_tumour_growth(1,:))).^2;
%     end
end

%% track the ESS
ESS_store = zeros(1,max_time);
variance = zeros(1,max_time);
%%
for t = 2:max_time
    t
    for i = 1:n
       p0_vec(t,i) = exp(log(p0_mu - p0_min) + p0_vec_normal(t-1,i)) + p0_min;
       psc_vec(t, i) = exp(log(psc_mu - psc_min) + psc_vec_normal(t-1,i)) + psc_min;
       dmax_vec(t, i) = exp(log(dmax_mu - dmax_min) + dmax_vec_normal(t-1,i)) + dmax_min; 
       gage_vec(t, i) = floor(exp(log(gage_mu - gage_min) + gage_vec_normal(t-1,i))) + gage_min;
       % simulate particle based on previous information
       % calculate the weight based on simulation
       particle{i} = clib.Model.CreateNewParticle(p0_vec(t,i),psc_vec(t, i), dmax_vec(t, i), gage_vec(t, i), page, particle{i});
       PFTumourVolume(t,i) = particle{i}.SimulateOneDay(1);
    end
    
    variance(1,t) = bisection(n,t,threshold,weight,Control_tumour_growth,PFTumourVolume);
    
    for i = 1:n
        weight(t,i) = log(weight(t-1,i)) + log_normpdf(Control_tumour_growth(t),PFTumourVolume(t,i),variance(1,t));
    end
    % Normalize the weight
    weight(t,:) = exp(weight(t,:) - logsumexp(weight(t,:),2));
    
    % get sample mean and covariance
    mean_theta_p0 = sum(weight(t,:) .* p0_vec_normal(t-1, :))
    mean_theta_psc = sum(weight(t,:) .* psc_vec_normal(t-1, :))
    mean_theta_dmax = sum(weight(t,:) .* dmax_vec_normal(t-1, :))
    mean_theta_gage = sum(weight(t,:) .* gage_vec_normal(t-1, :))
    
    V_p0 = weightedcov(p0_vec_normal(t-1, :)', weight(t,:)); 
    V_psc = weightedcov(psc_vec_normal(t-1, :)', weight(t,:));
    V_dmax = weightedcov(dmax_vec_normal(t-1, :)', weight(t,:));
    V_gage = weightedcov(gage_vec_normal(t-1, :)', weight(t,:));
    
    ESS = 1/sum(weight(t,:).^2);
    ESS_store(t) = ESS;
    
    ind = randsample(1:n, n, 'true', weight(t,:));
    p0_vec_normal(t-1,:) = p0_vec_normal(t-1,ind);
    psc_vec_normal(t-1,:) = psc_vec_normal(t-1,ind);
    dmax_vec_normal(t-1,:) = dmax_vec_normal(t-1,ind);
    gage_vec_normal(t-1,:) = gage_vec_normal(t-1,ind);
    
    % resample the simulate value
    PFTumourVolume(t,:) = PFTumourVolume(t,ind);
        
    weight(t,:) = ones(1,n)/n;
    
    for i = 1:n
        m_t_p0 = a*p0_vec_normal(t-1, i) + (1-a) * mean_theta_p0;
        m_t_psc = a*psc_vec_normal(t-1, i) + (1-a) * mean_theta_psc;
        m_t_dmax = a*dmax_vec_normal(t-1, i) + (1-a) * mean_theta_dmax;
        m_t_gage = a*gage_vec_normal(t-1, i) + (1-a) * mean_theta_gage;
        p0_vec_normal(t,i) = mvnrnd(m_t_p0,h2 * V_p0);
        psc_vec_normal(t,i) = mvnrnd(m_t_psc,h2 * V_psc);
        dmax_vec_normal(t,i) = mvnrnd(m_t_dmax,h2 * V_dmax);
        gage_vec_normal(t,i) = mvnrnd(m_t_gage,h2 * V_gage);
    end
    
end

%%
figure
subplot(2,2,1)
[f,xi] = ksdensity(p0_vec(1,:));
plot(xi,f,'b--')
hold on
[f,xi] = ksdensity(p0_vec(max_time,:));
plot(xi,f,'k')
plot(p0_true,0,'rx')
title('p_0')
legend('prior','posterior')

subplot(2,2,2)
[f,xi] = ksdensity(psc_vec(1,:));
plot(xi,f,'b--')
hold on
[f,xi] = ksdensity(psc_vec(max_time,:));
plot(xi,f,'k')
plot(psc_true,0,'rx')
title('psc')
legend('prior','posterior')

subplot(2,2,3)
[f,xi] = ksdensity(dmax_vec(1,:));
plot(xi,f,'b--')
hold on
[f,xi] = ksdensity(dmax_vec(max_time,:));
plot(xi,f,'k')
plot(dmax_true,0,'rx')
title('d_{max}')
legend('prior','posterior')

subplot(2,2,4)
[f,xi] = ksdensity(gage_vec(1,:));
plot(xi,f,'b--')
hold on
[f,xi] = ksdensity(gage_vec(max_time,:));
plot(xi,f,'k')
plot(gage_true,0,'rx')
title('g_{age}')
legend('prior','posterior')

%% plot 
time = linspace(1,max_time,max_time);

p0_vec_mean = zeros(1,max_time);
psc_vec_mean = zeros(1,max_time);
dmax_vec_mean = zeros(1,max_time);
gage_vec_mean = zeros(1,max_time);

for t = 1:max_time
    p0_vec_mean(t) = mean(p0_vec(t,:));
    psc_vec_mean(t) = mean(psc_vec(t,:));
    dmax_vec_mean(t) = mean(dmax_vec(t,:));
    gage_vec_mean(t) = mean(gage_vec(t,:));
end

fig = tiledlayout('flow','TileSpacing','compact');
nexttile%subplot(2,2,1)
quantil025_p0 = zeros(1,max_time);
%quantil95_p0 = zeros(1,max_time);
quantil975_p0 = zeros(1,max_time);
for t = 1:max_time
    quantil025_p0(t) = quantile(p0_vec(t,:),0.025);
    %quantil95_p0(t) = quantile(p0_vec(t,:),0.95);
    quantil975_p0(t) = quantile(p0_vec(t,:),0.975);
end
plot(time,p0_vec_mean);
hold on 
plot(time,quantil025_p0);
%plot(time,quantil95_p0)
plot(time,quantil975_p0);
title('parameter estimate for p_{0}');
xlabel('time');
ylabel('p_0');

nexttile%subplot(2,2,2)
quantil025_psc = zeros(1,max_time);
%quantil95_psc = zeros(1,max_time);
quantil975_psc = zeros(1,max_time);
for t = 1:max_time
    quantil025_psc(t) = quantile(psc_vec(t,:),0.025);
    %quantil95_psc(t) = quantile(psc_vec(t,:),0.95);
    quantil975_psc(t) = quantile(psc_vec(t,:),0.975);
end
plot(time,psc_vec_mean);
hold on 
plot(time,quantil025_psc);
%plot(time,quantil95_psc)
plot(time,quantil975_psc);
title('parameter estimate for psc');
xlabel('time');
ylabel('psc');

nexttile%subplot(2,2,3)
quantil025_dmax = zeros(1,max_time);
%quantil95_dmax = zeros(1,max_time);
quantil975_dmax = zeros(1,max_time);
for t = 1:max_time
    quantil025_dmax(t) = quantile(dmax_vec(t,:),0.025);
    %quantil95_dmax(t) = quantile(dmax_vec(t,:),0.95);
    quantil975_dmax(t) = quantile(dmax_vec(t,:),0.975);
end
plot(time,dmax_vec_mean);
hold on 
plot(time,quantil025_dmax);
%plot(time,quantil95_dmax);
plot(time,quantil975_dmax);
title('parameter estimate for d_{max}');
xlabel('time');
ylabel('d_{max}');

nexttile%subplot(2,2,4)
quantil025_gage = zeros(1,max_time);
%quantil95_gage = zeros(1,max_time);
quantil975_gage = zeros(1,max_time);
for t = 1:max_time
    quantil025_gage(t) = quantile(gage_vec(t,:),0.025);
    %quantil95_gage(t) = quantile(gage_vec(t,:), 0.95);
    quantil975_gage(t) = quantile(gage_vec(t,:),0.975);
end
plot(time,gage_vec_mean);
hold on 
plot(time,quantil025_gage);
%plot(time,quantil95_gage);
plot(time,quantil975_gage);
title('parameter estimate for g_{age}');
xlabel('time');
ylabel('g_{age}');
legend('mean','0.025 quantile','0.975 quantile');
lgd = legend;
lgd.Layout.Tile = 'east';

%% Prior predictive check 
s_n = n;
Simulate_prior_tumour = zeros(max_time,s_n);
startingvol = 200;

for i = 1:s_n
    pf = clib.Model.SeedAndGrowToStartVolume(p0_vec(1, i),psc_vec(1, i), dmax_vec(1, i), gage_vec(1, i), page, startingvol);
    for t = 1:max_time
        pf.SimulateOneDay(1);
        Simulate_prior_tumour(t,i) = pf.TumourVolume(); %calculates tumour volume of one particle's simulation
    end
end

% Posterior predictive check
s_n_pos = n;
Simulate_pos_tumour = zeros(max_time,s_n_pos);

%
for i = 1:s_n
    pf = clib.Model.SeedAndGrowToStartVolume(p0_vec(max_time, i),psc_vec(max_time, i), dmax_vec(max_time, i), gage_vec(max_time, i), page,startingvol);
    for t = 1:max_time 
        pf.SimulateOneDay(1);
        Simulate_pos_tumour(t,i) = pf.TumourVolume(); %calculates tumour volume of one particle's simulation
    end
end


%% Combine prior and posterior
figure
subplot(2,1,1)
plot(linspace(1,max_time,max_time),Control_tumour_growth(1,1:max_time),'LineWidth',2);
hold on
boxplot(Simulate_prior_tumour');
%boxplot(log(Simulate_pos_tumour));
xlabel('time')
ylabel('simulated tumour volumn')
title('Prior predictive check')

% Posterior
subplot(2,1,2)
plot(linspace(1,max_time,max_time),Control_tumour_growth(1,1:max_time),'LineWidth',2);
hold on
%boxplot(log(Simulate_prior_tumour)');
boxplot(Simulate_pos_tumour');
xlabel('time')
ylabel('simulated tumour volumn')
title('Posterior predictive check')

%% Kernel plot
% for t = 1:max_time
% figure 
% ksdensity(Simulate_prior_tumour(t,:));
% hold on
% ksdensity(Simulate_pos_tumour(t,:));
% legend('prior','posterior')
% title('time',t)
% end