% =========================================================================
% -- Script to compute exact average symbol error rate (SER) for LoRa
% =========================================================================

% -- (c) 2023 Ganghui Lin, Ahmed Elzanaty, Mohamed-Slim Alouini

% -- e-mail: ganghui.lin@kaust.edu.sa; a.elzanaty@surrey.ac.uk; slim.alouini@kaust.edu.sa

% =========================================================================

% G. Lin, A. Elzanaty, and M.-S. Alouini, "LoRa Backscatter Communications: Temporal, Spectral, and Error Performance Analysis,"
% in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2023.3268113.

% =========================================================================
% Example: 
% snr_vec = -20:1:-5; %snr vector in dB
% SF = 7;
% ser = SER_AWGN_LoRa_Num(snr_vec,SF);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function SER = SER_AWGN_LoRa_Num(snr_vec,SF)
M=2^SF;
SER=zeros(1,length(snr_vec));
ricepdf = @(x, v, s) (x ~= 0) .* (x ./ s.^2) .* exp(-(x.^2 + v.^2) ./ (2 .* s.^2)) .* besseli(0, x .* v ./ s.^2);
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    sigma=sqrt(1./snr/(M*2));
    F_rayleigh=@(x) 1-exp(-x.^2/2/sigma^2);
    f_rician=@(x) ricepdf(x,1,sigma);
    intfunc=@(x) (1-F_rayleigh(x).^(M-1)).*f_rician(x);
    SER(j)=integral(intfunc,0,5);
end
end

