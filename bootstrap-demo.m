%{ Bootstrap confidence interval estimation in Matlab

This script demonstrates the application of bootstrap in a nonlinear regression problem

Original code from: Lorena Ruiz del Rio
Modified by: Raibatak Das
Date: Apr 2016
%}

%% Part 1: Bootstrap demo with simulated data

% Load time points
t=xlsread('PBM.xlsx','A2:A104');
% Specify model parameters
Bo = 56.60; %mLCH4/gSV
Kh = 0.007;       %h^1
modelVar = [Bo Kh];
% Compute model prediction
Breal = Hidrolisis(modelVar, t);
%Adicionar error normal o "ruido" a la simualación de los experimentos
% Add noise to model prediction
NoiseStd = 1; % La desviación estándar (sd of noise)
err = NoiseStd*randn(length(Breal), 1);
Bexp = Breal + err; % This is simulated noisy data that we will fit

% Perform nonlinear regression on simulated data
options=statset('Display', 'iter');
% Initialize with reasonable parameter guesses
betaGuess = [max(Bexp) 0.1]; %Valores iniiciales 
% Extract best fit parameters and residuals 
[betaHat, residuals, ~, ~, ~] = nlinfit(t, Bexp, @Hidrolisis, betaGuess, options);
betaHat

% Plot simulated data and fit
figure('Position', [100 100 800 600], 'Color', 'w')
plot(t, Bexp, 'o', 'MarkerSize', 6) % Plot simulated data
hold on
tplot = linspace(0, max(t), 200);
Bplot = Hidrolisis(betaHat, tplot);
plot(tplot, Bplot,'-','LineWidth', 2) % Add best fit curve
plot(t, Breal, '--','LineWidth',1) % Add true model curve
xlabel('t (horas)')
ylabel('B (mLCH4/gSV)')
ylim([0, 60])
grid on
legend('Simulated data','Best fit', 'True model', 'Location', 'SouthEast')
hold off
% Uncomment below to save figure
%print('plots/sim-data-fit.png', '-dpng', '-r100')

% Bootstrap example: Generate a single bootrap dataset and fit 
% Resample residuals
[~, bootIndices] = bootstrp(1, [], residuals);
bootResiduals = residuals(bootIndices);
% Add resampled residuals to model prediction
Bmodel = Hidrolisis(betaHat, t); % Model prediction using estimated parameters
Bbootstrap = Bmodel + bootResiduals;
% Fit bootstrap data  
betaBoot = nlinfit(t, Bbootstrap, @Hidrolisis, betaGuess);
betaBoot

% Plot original data and bootstrap data + fits
clf()
plot(t, Bexp, 'o', 'MarkerSize', 4)
hold on
plot(t, Bmodel, '-', 'LineWidth', 2)
plot(t, Bbootstrap, 'o', 'MarkerSize', 4)
plot(t, Hidrolisis(betaBoot, t), '-', 'LineWidth', 2)
hold off
grid
ylim([0 60])
xlabel('t')
ylabel('B')
legend('Simulated data','Fit to simulated data', 'Bootstrap data', 'Fit to bootstrap data', ...
  'Location', 'SouthEast')
% Uncomment below to save figure
%print('plots/sim-bootstrap-fit.png', '-dpng', '-r100')

% Now, apply bootstrap procedure many times to build up the bootstrap distributions 
nboot= 1000 %Número de replicas para el Bootstrap
[~, bootIndices] = bootstrp(nboot, [], residuals);
bootResiduals = residuals(bootIndices);
Bbootstrap =repmat(Bmodel, 1, nboot) + bootResiduals;
% Fit each bootstrap datset and extract parameter estimates
betaBoot=zeros(nboot,2);
for i=1:nboot
    betaBoot(i,:)=nlinfit(t, Bbootstrap(:,i), @Hidrolisis, betaGuess);
end

% Estimate 95% bootstrap confidence intervals
%Estimación del Intervalo de Confianza del 95%
bootCI = prctile(betaBoot, [2.5 97.5])

% Plot bootstrap parameter distributions
clf()
subplot(211)
[freq loc] = hist(betaBoot(:,1), 30);
stairs(loc, freq, 'LineWidth', 2)
hold on
line([bootCI(1,1), bootCI(1,1)], get(gca, 'YLim'), 'Color', 'r')
line([bootCI(2,1), bootCI(2,1)], get(gca, 'YLim'), 'Color', 'r')
xlabel('Bo')
ylabel('Frequency')
title('Bootstrap distributions')
subplot(212)
[freq loc] = hist(betaBoot(:,2), 30);
stairs(loc, freq, 'LineWidth', 2)
hold on
line([bootCI(1,2), bootCI(1,2)], get(gca, 'YLim'), 'Color', 'r')
line([bootCI(2,2), bootCI(2,2)], get(gca, 'YLim'), 'Color', 'r')
xlabel('Kh')
ylabel('Frequency')
% Uncomment below to save figure
%print('plots/sim-bootstrap-distr.png', '-dpng', '-r100')

%% Part 2: Apply bootstrap to observed data

% Load observed values
Bensayo=xlsread('PBM.xlsx','B2:B104');
% Fit to model
betaGuess = [max(Bensayo), 0.1]
[betaHat, residuals, ~, ~, ~] = nlinfit(t, Bensayo, @Hidrolisis, betaGuess, options);
betaHat

% Plot data and model
clf()
plot(t, Bensayo, 'o', 'MarkerSize', 6) % Plot simulated data
hold on
tplot = linspace(0, max(t), 200);
Bplot = Hidrolisis(betaHat, tplot);
plot(tplot, Bplot,'-','LineWidth', 2) % Add best fit curve
xlabel('t (horas)')
ylabel('B (mLCH4/gSV)')
grid on
legend('Observed values','Fit', 'Location', 'SouthEast')
hold off
% Uncomment below to save figure
%print('plots/obs-data-fit.png', '-dpng', '-r100')

% Plot residuals
clf()
subplot(211)
plot(t, residuals, 'o:')
line(get(gca, 'XLim'), [0 0], 'Color', 'k')
xlabel('t')
ylabel('residuals')
subplot(212)
[freq loc] = hist(residuals, 20);
stairs(loc, freq, 'LineWidth', 2)
line([0 0], get(gca, 'YLim'), 'Color', 'k')
xlabel('residuals')
ylabel('frequency')
% Uncomment below to save figure
%print('plots/obs-residuals.png', '-dpng', '-r100')

% Apply bootstrap procedure to observed data
nboot= 1000 %Número de replicas para el Bootstrap
[~, bootIndices] = bootstrp(nboot, [], residuals);
bootResiduals = residuals(bootIndices);
Bmodel = Hidrolisis(betaHat, t);
Bbootstrap =repmat(Bmodel, 1, nboot) + bootResiduals;
% Fit each bootstrap datset and extract parameter estimates
betaBoot=zeros(nboot,2);
for i=1:nboot
    betaBoot(i,:)=nlinfit(t, Bbootstrap(:,i), @Hidrolisis, betaGuess);
end

% Estimate 95% bootstrap confidence intervals
bootCI = prctile(betaBoot, [2.5 97.5])

% Plot bootstrap parameter distributions
clf()
subplot(211)
[freq loc] = hist(betaBoot(:,1), 30);
stairs(loc, freq, 'LineWidth', 2)
hold on
line([bootCI(1,1), bootCI(1,1)], get(gca, 'YLim'), 'Color', 'r')
line([bootCI(2,1), bootCI(2,1)], get(gca, 'YLim'), 'Color', 'r')
xlabel('Bo')
ylabel('Frequency')
title('Bootstrap distributions')
subplot(212)
[freq loc] = hist(betaBoot(:,2), 30);
stairs(loc, freq, 'LineWidth', 2)
hold on
line([bootCI(1,2), bootCI(1,2)], get(gca, 'YLim'), 'Color', 'r')
line([bootCI(2,2), bootCI(2,2)], get(gca, 'YLim'), 'Color', 'r')
xlabel('Kh')
ylabel('Frequency')
% Uncomment below to save figure
%print('plots/obs-bootstrap-distr.png', '-dpng', '-r100')

% Plot data, model and bootstrap prediction intervals
clf()
plot(t, Bensayo, 'o', 'MarkerSize', 6) % Plot simulated data
hold on
tplot = linspace(0, max(t), 200);
Bplot = Hidrolisis(betaHat, tplot);
plot(tplot, Bplot,'-','LineWidth', 2) % Add best fit curve
% Bootstrap  lower bound
Blo = Hidrolisis(bootCI(1,:), tplot);
plot(tplot, Blo, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
% Bootstrap upper bound
Bhi = Hidrolisis(bootCI(2,:), tplot);
plot(tplot, Bhi, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
xlabel('t (horas)')
ylabel('B (mLCH4/gSV)')
grid on
hold off
% Uncomment below to save figure
%print('plots/bootstrap-prediction.png', '-dpng', '-r100')

