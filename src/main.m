close all
clc
figure('Visible', 'on')

%% Task 1: Construct X[k]
% Set the values of lambda
lambdas = [0.12, 0.3, 0.42, 0.19, 0.17, 0.36]';

% Define the range of k and N
N = 8192;
k = 1:N;

% Calculate the values of omega
omegas = 2*pi*lambdas;

% Generate the random phase offsets
phis = (2*pi*rand(6,1))';
phis(3) = phis(1)+phis(2);
phis(6) = phis(4)+phis(5);

% Calculate the values of cos(omega_i * k + phi_i) for all i and k
cos_values = zeros(6, N);
for i = 1:6
    cos_values(i, :) = cos(omegas(i)*k + phis(i));
end

% Sum the values of cos(omega_i * k + phi_i) for all i to get X[k]
X = sum(cos_values);

% Define a figure counter
figCounter = 1;

% Plot X[k]
figure(figCounter)
figCounter = figCounter + 1;
plot(X)
title("X[k]")
ylabel("X[k]")
xlabel("k")

% Plot a smaller sample of our data for a clearer review.
figure(figCounter)
figCounter = figCounter + 1;
plot(X(1:800))
title("Smaller Sample of X[k]")
ylabel("X[k]")
xlabel("k")

%% Autocorrelation

% Maximum number of shiftings
L2 = 128;  

% Cross-correlation including negative lags
crossCorrFull = xcorr(X, L2);

% Normalizing the cross-correlation values
normCrossCorrFull = crossCorrFull / max(crossCorrFull); 

% Autocorrelation of the signal
autoCorr = autocorr(X, L2);                

% Plot Autocorrelation
autocorr(X, L2);
timeAxis = -L2:1:L2;
figure(figCounter)
figCounter = figCounter + 1;
plot(timeAxis, normCrossCorrFull)
title('Autocorrelation of Signal')
xlabel('Time [sec]')
ylabel('ACF')
grid on

%% Task 2: Power Spectrum Estimation
% Power Spectrum using Autocorrelation
powerSpec = abs(fft(autoCorr));

% Getting the Frequency Axis for Power Spectrum
freqAxis = (0:length(powerSpec)-1) / length(powerSpec);

% Plot Power Spectrum
figure(figCounter)
figCounter = figCounter + 1;
plot(freqAxis, powerSpec)
title("Power Spectrum Estimation")
ylabel("Pxx[f]")
xlabel("f")
% Adding a line at frequency 0.5
line([0.5 0.5], [0 12], 'Color', 'red', 'LineStyle', '--')
hold on
% Finding locations of the 6 highest peaks
[peaks, peakLocations] = findpeaks(powerSpec, freqAxis, 'NPeaks', 6);
% Plotting the peaks on the power spectrum
findpeaks(powerSpec, freqAxis, 'NPeaks', 6);
% Displaying the peak numbers
text(peakLocations(1:6) + 0.02, peaks(1:6) - 0.5, num2str([1; 5; 4; 2; 6; 3]))
% Displaying the peak frequencies
text(peakLocations(1:6), peaks(1:6) + 0.8, num2str(peakLocations(1:6)'))

%% Task 3: Bispectrum Estimation
% N = 8192;
K = 32;
M = 256;
L = 64;

% Split data
subsetsX = reshape(X,[M,K]);

% a1) Indirect Method with Rectangular Window for K=32, M=256, L3=64
figure(figCounter)
figCounter = figCounter + 1;
BispecInRec = indirectBispectrum(subsetsX, L, M, 0, 'unbiased', 128, 1);
hold on;
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Indirect Method with Rectangular Window');
legend('Bispectrum','Primary Area');

% a2) Indirect Method with Parzen Window for K=32, M=256, L3=64
figure(figCounter)
figCounter = figCounter + 1;
BispecInParz = indirectBispectrum(subsetsX, L, M, 0, 'unbiased', 128, 0);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Indirect Method with Parzen Window');
legend('Bispectrum','Primary Area');


% b) Direct Method for K=32, M=256, J=0
J = 0;
D = 2*J + 1;

figure(figCounter)
figCounter = figCounter + 1;
bispecDir = directBispectrum(subsetsX, M, D, M, 0);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Direct Method');
legend('Bispectrum','Primary Area');

%% Task 7a: Change Segment Length
% i) K=16, M=512
K = 16;
M = 512;

% Split Data
subsetsX = reshape(X,[M,K]);

% Indirect Method with Rectangular Window for K=16, M=512, L3=64
figure(figCounter)
figCounter = figCounter + 1;
bispecInRec1 = indirectBispectrum(subsetsX, L, M, 0, 'unbiased', 128, 1);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Indirect Method with Rectangular Window');
legend('Bispectrum','Primary Area');

% Indirect Method with Parzen Window for K=16, M=512, L3=64
figure(figCounter)
figCounter = figCounter + 1;
bispecInParz1 = indirectBispectrum(subsetsX, L, M, 0, 'unbiased', 128, 0);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Indirect Method with Parzen Window');
legend('Bispectrum','Primary Area');

% Direct Method for K=16, M=512, J=0
figure(figCounter)
figCounter = figCounter + 1;
bispecDir1 = directBispectrum(subsetsX, M, D, M, 0);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Direct Method');
legend('Bispectrum','Primary Area');

% ii) K=64, M=128
K = 64;
M = 128;

% Split Data
subsetsX = reshape(X, [M,K]);

% Indirect Method with Rectangular Window for K=64, M=128, L3=64
figure(figCounter)
figCounter = figCounter + 1;
bispecInRec2 = indirectBispectrum(subsetsX, L, M, 0, 'unbiased', 128, 1);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Indirect Method with Rectangular Window');
legend('Bispectrum','Primary Area');

% Indirect Method with Parzen Window for K=64, M=128, L3=64
figure(figCounter)
figCounter = figCounter + 1;
bispecInParz2 = indirectBispectrum(subsetsX, L, M, 0, 'unbiased', 128, 0);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Indirect Method with Parzen Window');
legend('Bispectrum','Primary Area');

% Direct Method for K=64, M=128, J=0
figure(figCounter)
figCounter = figCounter + 1;
bispecDir2 = directBispectrum(subsetsX, M, D, M, 0);
% Plot primary area, where f1=f2, f1+f2=0.5, f2=0
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');
plot([0.25,0.5],[0.25,0],'Color','#D95319');
plot([0,0.5],[0,0],'Color','#D95319');
title('Bispectrum Estimation using the Direct Method');
legend('Bispectrum','Primary Area');

%% Task 7b: Realization and Mean Value Comparison
K = 32;
M = 256;
R = 50;

% Compute the Fast Fourier Transform (FFT) of the autoCorr
fftResult = fft(autoCorr);

% Initialize variables for mean calculations
meanC2 = zeros(length(fftResult), 1);  % Mean of C2
meanC3In1 = zeros(M, M);  % Mean of C3 In1
meanC3In2 = zeros(M, M);  % Mean of C3 In2
meanC3Dir = zeros(M, M);  % Mean of C3 Dir

for i=1:R
    % Generate the random phase offsets
    phis = (2*pi*rand(6,1))';
    phis(3) = phis(1)+phis(2);
    phis(6) = phis(4)+phis(5);
    
    % Calculate the values of cos(omega_j * k + phi_j) for all j and k
    cos_values1 = zeros(6, N);
    for j = 1:6
        cos_values1(j, :) = cos(omegas(j)*k + phis(j));
    end
    
    % Sum the values of cos(omega_j * k + phi_j) for all j to get X[k]
    X = sum(cos_values1);
    
    % Autocorrelation
    L2 = 128;                                      % Maximum number of shiftings
    autoCorr1 = autocorr(X, L2);                % Autocorrelation of the signal
    
    % Power Spectrum using Autocorrelation
    powerSpec1 = abs(fft(autoCorr1));
    
    % Getting the Frequency Axis for Power Spectrum
    freqAxis1 = (0:length(powerSpec1)-1) / length(powerSpec1);
    
    meanC2 = meanC2 + powerSpec1;
    
    % Split data
    subsetsX = reshape(X,[M,K]);
    
    % Indirect Method with Rectangular Window
    bispecInRec3 = indirectBispectrum(subsetsX,L,M,0,'unbiased',128,1);
    meanC3In1 = meanC3In1 + bispecInRec3;
    
    % Indirect Method with Parzen Window
    bispecInParz3 = indirectBispectrum(subsetsX,L,M,0,'unbiased',128,0);
    meanC3In2 = meanC3In2 + bispecInParz3;
    
    % Direct Method
    bispecDir3 = directBispectrum(subsetsX,M,D,M,0);
    meanC3Dir = meanC3Dir + bispecDir3;
end

% Update mean values
meanC2 = meanC2/R;
meanC3In1 = meanC3In1/R;
meanC3In2 = meanC3In2/R;
meanC3Dir = meanC3Dir/R;

nfft = 256;
% Rotate for proper orientation
if rem(nfft,2) == 0
    waxis = (-nfft/2:(nfft/2-1))/nfft;
else
    waxis = (-(nfft-1)/2:(nfft-1)/2)/nfft;
end

% Plot mean power spectrum
figure(figCounter)
figCounter = figCounter + 1;
plot(freqAxis1, meanC2(1,:))
title("Mean Power Spectrum Estimation")
ylabel("Mean Pxx[f]")
xlabel("f")
% Adding a line at frequency 0.5
line([0.5 0.5], [0 12], 'Color', 'red', 'LineStyle', '--')
hold on
% Finding locations of the 6 highest peaks
[peaks1, peakLocations1] = findpeaks(meanC2(1,:), freqAxis1, 'NPeaks', 6);
% Plotting the peaks on the power spectrum
findpeaks(meanC2(1,:), freqAxis1, 'NPeaks', 6);
% Displaying the peak numbers
text(peakLocations1(1:6) + 0.02, peaks1(1:6) - 0.5, num2str([1; 5; 4; 2; 6; 3]))
% Displaying the peak frequencies
text(peakLocations1(1:6), peaks1(1:6) + 0.8, num2str(peakLocations1(1:6)'))

% Plot the mean bispectrum using indirect method with rectangular window
figure(figCounter)
figCounter = figCounter + 1;
contour(waxis, waxis, abs(meanC3In1), 4); 
grid on;
title('Mean Bispectrum Estimation using the Indirect Method with Rectangular Window');
xlabel('f1');
ylabel('f2');

% Plot the bispectrum using indirect method with Parzen window
figure(figCounter)
figCounter = figCounter + 1;
contour(waxis, waxis, abs(meanC3In2), 4); 
grid on;
title('Mean Bispectrum Estimation using the Indirect Method with Parzen Window');
xlabel('f1');
ylabel('f2');

% Plot the bispectrum using direct method
figure(figCounter)
figCounter = figCounter + 1;
contour(waxis, waxis, abs(meanC3Dir), 4); 
grid on;
title('Mean Bispectrum Estimation using the Direct Method');
xlabel('f1');
ylabel('f2');