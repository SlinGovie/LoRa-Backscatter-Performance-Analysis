% =========================================================================
% -- Script to compute approximated average symbol error rate (SER) for
% LoRa backscatter using Gaussian Hermite quadrature
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
%% Params, signal, and decoder
M=2^SF;
phi=zeros(M,M);
fft_LB=zeros(M,M);
E_s=1;
SER=zeros(length(snr_vec),1);

for a=1:M
    for k=1:M
        phi(a,k)=(k-1)/M*(2*(a-1)-M+k-1);
        while(phi(a,k)<-1)
            phi(a,k)=phi(a,k)+2;
        end
    end
end

k=0:M-1;
downchirp=(E_s/M)^0.5*exp(-1j*pi*k.^2/M+1j*pi.*k);


Q=(floor(phi*2^(N-1))+0.5)/2^(N-1);
x_LB=(E_s/M)^0.5*exp(1j*pi*Q);
x_conj=conj(x_LB).';

for i=1:M
    fft_LB(i,:)=abs(fft(x_LB(i,:).*downchirp,M));
end

ricecdf = @(x, v, s) (x > 0) .* 1-marcumq(v./s,(x>0).*x./s);

%% roots and weights of hermite polynomial
hermite_root=roots(hermite(N_GH));
wi=(2^(N_GH-1)*factorial(N_GH)*pi^0.5)/N_GH^2./hermiteH(N_GH-1,hermite_root).^2;
%% ML decoder
if  strcmp(decoder,'ML')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    s=1./snr/(M*2);%variance
    Fl=zeros(1,length(hermite_root));
    SEP=zeros(1,M);
    ML_LB=abs(x_LB*x_conj);
    cdfcell=zeros(M,M);
    for a=1:M
        in=-ML_LB(a,a)^2/2/s;
        mu=(s*pi/2)^0.5*exp(in/2)*((1-in)*besseli(0,-in/2)-in*besseli(1,-in/2));
        sigma2=2*s+ML_LB(a,a)^2-mu^2;
        xin=2^0.5*sigma2^0.5*hermite_root+mu;
        for t=1:length(xin)
        for i=1:M
            cdfcell(a,i)=ricecdf(xin(t),ML_LB(a,i),s^0.5);
        end
        cdfcell(a,a)=1;
        Fl(t)=prod(cdfcell(a,:));
        SEP(a)=1-1/pi^0.5*sum(Fl.*wi.');
        end
    end
    SER(j)=1/M*sum(SEP);
end
%% FFT decoder
elseif strcmp(decoder,'fft')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    s=1./snr/(M*2); %variance
    fft_LB=zeros(M,M);
    Q=(floor(phi*2^(N-1))+0.5)/2^(N-1);
    x_LB=(E_s/M)^0.5*exp(1j*pi*Q);
    for i=1:M
    fft_LB(i,:)=abs(fft(x_LB(i,:).*downchirp,M));
    end
    Fl=zeros(1,length(hermite_root));
    SEP=zeros(1,M);
    
    cdfcell=zeros(M,M);
    for a=1:M
        in=-fft_LB(a,a)^2/2/s;
        mu=(s*pi/2)^0.5*exp(in/2)*((1-in)*besseli(0,-in/2)-in*besseli(1,-in/2));
        sigma2=2*s+fft_LB(a,a)^2-mu^2;
        xin=2^0.5*sigma2^0.5*hermite_root+mu;
        for t=1:length(xin)
        for i=1:M
            cdfcell(a,i)=ricecdf(xin(t),fft_LB(a,i),s^0.5);
        end
        cdfcell(a,a)=1;
        Fl(t)=prod(cdfcell(a,:));
        SEP(a)=1-1/pi^0.5*sum(Fl.*wi.');
        end
    end
    SER(j)=1/M*sum(SEP);
end

else 
disp('please input correct decoder')
end

end

