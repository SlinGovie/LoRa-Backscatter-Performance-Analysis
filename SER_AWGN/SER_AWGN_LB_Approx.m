% =========================================================================
% -- Script to compute approximated average symbol error rate (SER) for
% LoRa backscatter using Gaussian Hermite quadrature in AWGN channel
% =========================================================================

% -- (c) 2023 Ganghui Lin, Ahmed Elzanaty, Mohamed-Slim Alouini

% -- e-mail: ganghui.lin@kaust.edu.sa; a.elzanaty@surrey.ac.uk; slim.alouini@kaust.edu.sa

% =========================================================================

% G. Lin, A. Elzanaty, and M.-S. Alouini, "LoRa Backscatter Communications: Temporal, Spectral, and Error Performance Analysis,"
% in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2023.3268113.

% =========================================================================
% Example: 
% snr_vec = -20:1:-5; %snr vector in dB
% SF = 7; % spreading factor
% N = 2; % 2^N number of loads
% decoder = 'ML';
% decoder = 'fft';
% N_GH = 20; % # of function samples to approximate the integral
% ser = SER_AWGN_LB_Approx(snr_vec,SF,N,decoder,N_GH);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function SER = SER_AWGN_LB_Approx(snr_vec,SF,N,decoder,N_GH)
%% Params, functions
M=2^SF;
SER = zeros(1,length(snr_vec));
ricecdf = @(x, v, s) (x > 0) .* 1-marcumq(v./s,(x>0).*x./s);
%% roots and weights of hermite polynomial
syms x
p1=hermiteH(N_GH,x);
xt=roots(sym2poly(p1));
wt=(2^(N_GH-1)*factorial(N_GH)*sqrt(pi))/N_GH^2./hermiteH(N_GH-1,xt).^2;
%% SER calculation using GH quadrature
if strcmp(decoder, 'ML') || strcmp(decoder, 'fft')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    s=1./snr/(M*2);%variance
    Fl=zeros(1,length(xt));
    SEP=zeros(1,M);
    ML_LB=abs(decoder_output(SF,N,decoder));
    cdfcell=zeros(M,M);
    for a=1:M
        in=-ML_LB(a,a)^2/2/s;
        mu=(s*pi/2)^0.5*exp(in/2)*((1-in)*besseli(0,-in/2)-in*besseli(1,-in/2));
        sigma2=2*s+ML_LB(a,a)^2-mu^2;
        xin=2^0.5*sigma2^0.5*xt+mu;
        for t=1:length(xin)
            for i=1:M
                cdfcell(a,i)=ricecdf(xin(t),ML_LB(a,i),s^0.5);
            end
            cdfcell(a,a)=1;
            Fl(t)=prod(cdfcell(a,:));
            SEP(a)=1-1/pi^0.5*sum(Fl.*wt.');
        end
    end
    SER(j)=1/M*sum(SEP);
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end

else 
disp('please input correct decoder')
end

end

