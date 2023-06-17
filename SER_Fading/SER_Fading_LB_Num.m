% =========================================================================
% -- Script to compute average symbol error rate (SER) for
% LoRa backscatter by numerical integration
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
% decoder = 'ML';
% decoder = 'fft';
% d = 10; % distance between Tx & Rx, does not represent the real distance
% ratio = 1; % ratio = d1/d2
% m = [10,2]; % m = [m1,m2], shape parameters of nakagami-m distribution
% ser = SER_Fading_LB_Num(snr_vec,SF,N,decoder,d,ratio,m);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function ser=SER_Fading_LB_Num(snr_vec,SF,N,decoder,d,ratio,m)
%% Params, functions
Es=1;
M=2^SF;
d1=d/(1+ratio);
d2=d-d1;
Omega1=Es/d1^2;
Omega2=Es/d2^2;
m1=m(1);
m2=m(2);

r1=m1/Omega1;r2=m2/Omega2;v=m1+m2;
fH=@(x,h) 4/(gamma(m1)*gamma(m2))*(r1*r2)^(v/2)*h.^(v-1).*besselk(m1-m2,2*sqrt(r1*r2)*h);

ser=zeros(1,length(snr_vec));
%% ML decoder
if strcmp(decoder, 'ML') || strcmp(decoder, 'fft')
    ML=decoder_output(SF,N,decoder);
    asum=zeros(1,M);
    for a=1:M
        asum(a)=sum(abs(ML(a,:)));
    end
        S=sort(uniquetol(abs(asum),10^(-6)));
        Na=zeros(length(S),1);
        atil=zeros(length(S),1);
    for q=1:length(S)
            fd=find(abs(asum-S(q))<1e-5);
            Na(q)=length(fd);
            atil(q)=fd(1);
    end

    parfor j=1:length(snr_vec)
        ML=decoder_output(SF,N,decoder);
        snr=10.^(snr_vec(j)/10);
        sigma_n=sqrt(1./snr/(M*2)); 
        SEP=zeros(1,length(atil));
        for a=1:length(atil)
            kappa_neg=@(x,h)-(h*abs(ML(atil(a),atil(a)))).^2/(2*sigma_n^2);%-v^2/2sigma_n^2
            mu_g=@(x,h) sigma_n*sqrt(pi/2)*laguerreL(1/2,kappa_neg(x,h));%mu_g
            sigma_g=@(x,h)sqrt(2*sigma_n^2+(h*abs(ML(atil(a),atil(a)))).^2-mu_g(x,h).^2);
    
            flD=@(x,h) pdf('Normal',x,mu_g(x,h),sigma_g(x,h));
            ml=ML(atil(a),:);
            S=sort(uniquetol(abs(ml)));
            V=S(1:end-1);
            Ni=zeros(length(V),1);
            flhat=@(x,h) 1;
            for i=1:length(V)
                Ni(i)=length(find(abs(abs(ml)-V(i))<1e-5));
                cdfi=@(x,h) 1-marcumq(h*V(i)/sigma_n,x/sigma_n);
                cdfi=@(x,h) cdfi(x,h).^Ni(i);
                flhat=@(x,h) flhat(x,h).*cdfi(x,h);
            end
            flhat=@(x,h) 1-flhat(x,h)
            int=@(x,h) flhat(x,h).*flD(x,h).*fH(x,h);
            SEP(a)=integral2(int,0,5,0,5);
        end
        ser(j)=1/M*(SEP*Na);
    end

else
    disp('please input correct decoder')
end

end