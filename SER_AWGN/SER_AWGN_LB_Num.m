% =========================================================================
% -- Script to compute exact average symbol error rate (SER) for LoRa
% backscatter
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
% ser = SER_AWGN_LB_Num(snr_vec,SF,N,decoder);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function SER = SER_AWGN_LB_Num(snr_vec,SF,N,decoder)
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
x_conj=x_LB';

for i=1:M
    fft_LB(i,:)=abs(fft(x_LB(i,:).*downchirp,M));
end
ricepdf = @(x, v, s) (x > 0) .* (x ./ s.^2) .* exp(-(x.^2 + v.^2) ./ (2 .* s.^2)) .* besseli(0, x .* v ./ s.^2);
ricecdf = @(x, v, s) (x > 0) .* 1-marcumq(v./s,x./s);
%% ML decoder
if  strcmp(decoder,'ML')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    s2=1./snr/(M*2); % variance
    SEP_a=zeros(M,1); % error probability of symbol a
    parfor a=1:M
        ML_LB=abs(x_LB*x_conj);
        pdf=@(beta) ricepdf(beta,ML_LB(a,a),s2^0.5);
        pcell=@(beta) 1;
        S=sort(uniquetol(ML_LB(a,:)));
        V=S(1:end-1);
        Nv=zeros(length(V),1);
        for i=1:length(V)
            Nv(i)=length(find(abs(ML_LB(a,:)-V(i))<1e-5));
            cdfi=@(beta) ricecdf(beta,V(i),s2^0.5);
            cdfi=@(beta) cdfi(beta).^Nv(i);
            pcell=@(beta) pcell(beta).*cdfi(beta);
        end
        pcell=@(beta) pcell(beta).*pdf(beta);
        SEP_a(a)=1-integral(pcell,0,5);
    end
    SER(j)=1/M*sum(SEP_a);
     if isnan(SER(j))
        SER(j)=0;
        break;
     end
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end
%% FFT decoder
elseif strcmp(decoder,'fft')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    s2=1./snr/(M*2); % variance
    SEP_a=zeros(M,1);
    parfor a=1:M
        fft_LB=zeros(M,M);
        x_LB=(E_s/M)^0.5*exp(1j*pi*Q);
        for i=1:M
            fft_LB(i,:)=abs(fft(x_LB(i,:).*downchirp,M));
        end
        pdf=@(beta) ricepdf(beta,fft_LB(a,a),s2^0.5);
        pcell=@(beta) 1;
        S=sort(uniquetol(fft_LB(a,:)));
        V=S(1:end-1);
        Nv=zeros(length(V),1);
        for i=1:length(V)
            Nv(i)=length(find(abs(fft_LB(a,:)-V(i))<1e-5));
            cdfi=@(beta) ricecdf(beta,V(i),s2^0.5);
            cdfi=@(beta) cdfi(beta).^Nv(i);
            pcell=@(beta) pcell(beta).*cdfi(beta);
        end
        pcell=@(beta) pcell(beta).*pdf(beta);
        SEP_a(a)=1-integral(pcell,0,5);
    end
    SER(j)=1/M*sum(SEP_a);
    if isnan(SER(j))
        SER(j)=0;
        break;
    end
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end
else 
    disp('please input correct decoder')
end

end