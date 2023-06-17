% =========================================================================
% -- Script to compute exact average symbol error rate (SER) for LoRa
% backscatter in AWGN channel
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
%% Params, functions
M=2^SF;
SER=zeros(length(snr_vec),1);

ricepdf = @(x, v, s) (x > 0) .* (x ./ s.^2) .* exp(-(x.^2 + v.^2) ./ (2 .* s.^2)) .* besseli(0, x .* v ./ s.^2);
ricecdf = @(x, v, s) (x > 0) .* 1-marcumq(v./s,x./s);
%% SER calculation numerical integration
if strcmp(decoder, 'ML') || strcmp(decoder, 'fft')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    s2=1./snr/(M*2); % variance
    SEP_a=zeros(M,1); % error probability of symbol a
    parfor a=1:M
        ML_LB=abs(decoder_output(SF,N,decoder));
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
else 
    disp('please input correct decoder')
end

end