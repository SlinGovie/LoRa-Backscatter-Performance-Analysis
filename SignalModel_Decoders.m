% =========================================================================
% -- Script to generate LoRa backscatter signals and decoder/demodulator
% =========================================================================

% -- (c) 2023 Ganghui Lin, Ahmed Elzanaty, Mohamed-Slim Alouini

% -- e-mail: ganghui.lin@kaust.edu.sa; a.elzanaty@surrey.ac.uk; slim.alouini@kaust.edu.sa

% =========================================================================

% G. Lin, A. Elzanaty, and M.-S. Alouini, "LoRa Backscatter Communications: Temporal, Spectral, and Error Performance Analysis,"
% in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2023.3268113.

% =========================================================================
%% 
clc
clear
close all
%% Parameters
SF = 6; % spreading factor
M = 2^SF;
N = 3; % 
E_s = 1; % symbol energy

%% SignalModel
% sampled phases of LoRa
phi=zeros(M,M);
for a=1:M
    for k=1:M
        phi(a,k)=(k-1)/M*(2*(a-1)-M+k-1);
        while(phi(a,k)<-1)
            phi(a,k)=phi(a,k)+2;
        end
    end
end

% LoRa signal
x_LoRa = sqrt(E_s/M)*exp(1j*pi*phi);

% Phase quantization
Q=(floor(phi*2^(N-1))+0.5)/2^(N-1);

% Sampled LoRa backscatter 
x_LB=sqrt(E_s/M)*exp(1j*pi*Q); % each row represents the samples of a symbol

%% Decoders
%% ML
% cross-correlation of LoRa symbols
ML_LoRa = abs(x_LoRa*x_LoRa');

% ML decoder output without noise (cross correlation of symbols) for LB
ML_LB=abs(x_LB*x_LB');


%% FFT 
% downchirp
k=0:M-1;
downchirp=(E_s/M)^0.5*exp(-1j*pi*k.^2/M+1j*pi.*k);

% FFT decoder output without noise for LoRa
fft_LoRa=zeros(M,M);
for i=1:M
    fft_LoRa(i,:)=abs(fft(x_LoRa(i,:).*downchirp,M));
end

% FFT decoder output without noise for LB
fft_LB=zeros(M,M);
for i=1:M
    fft_LB(i,:)=abs(fft(x_LB(i,:).*downchirp,M));
end


%% Plots
% ML decoder for LoRa and LB
figure;
sgtitle('ML decoder output without noise')

subplot(1,2,1);
pcolor(ML_LB);
colorbar;
xlabel('a');
ylabel('i');
title('LoRa backscatter');

subplot(1,2,2);
pcolor(ML_LoRa);
colorbar;
xlabel('a');
ylabel('i');
title('LoRa');

% FFT decoder for LoRa and LB
figure;
sgtitle('FFT decoder output without noise')

subplot(1,2,1);
pcolor(fft_LB);
colorbar;
xlabel('a');
ylabel('i');
title('LoRa backscatter');

subplot(1,2,2);
pcolor(fft_LoRa);
colorbar;
xlabel('a');
ylabel('i');
title('LoRa');
