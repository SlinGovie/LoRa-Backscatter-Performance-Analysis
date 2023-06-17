% =========================================================================
% -- Script to compute average symbol error rate (SER) for LoRa backscatter
% using Monte Carlo simulation in nakagami-m fading channel
% =========================================================================

% -- (c) 2023 Ganghui Lin, Ahmed Elzanaty, Mohamed-Slim Alouini

% -- e-mail: ganghui.lin@kaust.edu.sa; a.elzanaty@surrey.ac.uk; slim.alouini@kaust.edu.sa

% =========================================================================

% G. Lin, A. Elzanaty, and M.-S. Alouini, "LoRa Backscatter Communications: Temporal, Spectral, and Error Performance Analysis,"
% in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2023.3268113.

% =========================================================================
% Example: 
% snr_vec = 5:5:25; %snr vector in dB
% SF = 7; % spreading factor
% N = 2; % 2^N number of loads
% N_s = 1e5; % # of random symbols for MC simulation
% decoder = 'ML';
% decoder = 'fft';
% d = 10; % distance between Tx & Rx, does not represent the real distance
% ratio = 1; % ratio = d1/d2
% m = [10,2]; % m = [m1,m2], shape parameters of nakagami-m distribution
% ser = SER_Fading_LB_MC(snr_vec,N_s,SF,N,decoder,d,ratio,m);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function ser=SER_Fading_LB_MC(snr_vec,N_s,SF,N,decoder,d,ratio,m)
%% Params
Es=1;
cnt=0;
ser=zeros(1,length(snr_vec));
M=2^SF;
d1=d/(1+ratio);
d2=d-d1;
Omega1=Es/d1^2;
Omega2=Es/d2^2;
m1=m(1);
m2=m(2);
pd1 = makedist('Nakagami','mu',m1,'omega',Omega1);
pd2 = makedist('Nakagami','mu',m2,'omega',Omega2);
%% ML decoder
if strcmp(decoder,'ML')
phi=zeros(M,M);
for a=1:M
    for k=1:M
        phi(a,k)=(k-1)/M*(2*(a-1)-M+k-1);
        while(phi(a,k)<-1)
            phi(a,k)=phi(a,k)+2;
        end
    end
end
Q=(floor(phi*2^(N-1))+0.5)/2^(N-1);
x_LB=sqrt(1/M)*exp(1j*pi*Q);
x_conj=conj(x_LB).';
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    sigma_n=sqrt(1./snr/(M*2)); 
    parfor i=1:N_s
        a=randi([1,M]);
        h1=random(pd1)
        h2=random(pd2)
        ML_LB=x_LB * x_conj
        noiseML=(randn(1,M)+1j*randn(1,M))*sigma_n;
        L_ML=h1*h2*ML_LB(a,:)+(x_conj*noiseML.').';
        [~,ahat]=max(abs(L_ML));
        cnt=cnt+(a~=ahat);
    end
    ser(j)=cnt/N_s;
    cnt=0;

    if ser(j)==0
        break
    end
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end
%% FFT decoder
elseif strcmp(decoder,'fft')
    for j=1:length(snr_vec)
        snr=10.^(snr_vec(j)/10);
        sigma_n=(1./snr/(M*2))^0.5; 
        parfor i=1:N_s
            a=randi([1,M]);
            h1=random(pd1)
            h2=random(pd2)
            fft_LB=decoder_output(SF,N,'fft');
            noisey=(randn(1,M)+1j*randn(1,M))*sigma_n;
            L_y=h1*h2*fft_LB(a,:)+noisey;
            [~,ahat]=max(abs(L_y));
            cnt=cnt+(a~=ahat);
        end
    ser(j)=cnt/N_s;
    cnt=0;
    
    if ser(j)==0
        break
    end
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
    end  

else
    disp('please input correct decoder')
end

end