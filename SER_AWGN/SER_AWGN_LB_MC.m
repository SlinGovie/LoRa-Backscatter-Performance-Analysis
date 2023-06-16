% =========================================================================
% -- Script to compute approximated average symbol error rate (SER) for
% LoRa backscatter using Monte Carlo simulation
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
% N_s = 1e4; % # of random symbols for MC simulation
% decoder = 'ML';
% decoder = 'fft';
% ser = SER_AWGN_LB_MC(snr_vec,N_s,SF,N,decoder);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function SER=SER_AWGN_LB_MC(snr_vec,N_s,SF,N,decoder)
%% Params
cnt=0;
SER=zeros(1,length(snr_vec));
M=2^SF;
phi=zeros(M,M);
E_s=1;
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

%% ML decoder
if strcmp(decoder,'ML')
for j=1:length(snr_vec)
    snr=snr_vec(j);
    parfor i=1:N_s
        a=randi([1,M]);
        x_LB=(E_s/M)^0.5*exp(1j*pi*Q);
        x_conj=x_LB';
        x_LB(a,:)=awgn(x_LB(a,:),snr,'measured');
        ML=abs(x_LB*x_conj);
        [~,ahat]=max(abs(ML(a,:)));
        cnt=cnt+(a~=ahat);
    end
    SER(j)=cnt/N_s;
    cnt=0;
    
    if SER(j)==0
        break
    end
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end
%% FFT decoder
elseif strcmp(decoder,'fft')
for j=1:length(snr_vec)
    snr=snr_vec(j);
    parfor i=1:N_s
        a=randi([1,M]);
        x_LB=(E_s/M)^0.5*exp(1j*pi*Q);
        x_LB(a,:)=awgn(x_LB(a,:),snr,'measured');
        Y=fft(x_LB(a,:).*downchirp,M);
        [~,ahat]=max(abs(Y));
        cnt=cnt+(a~=ahat);
    end
    SER(j)=cnt/N_s;
    cnt=0;
    if SER(j)==0
        break
    end
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end  

else
    disp('please input correct decoder')
end

end