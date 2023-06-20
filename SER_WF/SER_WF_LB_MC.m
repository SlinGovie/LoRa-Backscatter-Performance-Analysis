% =========================================================================
% -- Script to compute average symbol error rate (SER) for LoRa backscatter
% using Monte Carlo simulation in nakagami-m fading channel with water
% filling
% =========================================================================

% -- (c) 2023 Ganghui Lin, Ahmed Elzanaty, Mohamed-Slim Alouini

% -- e-mail: ganghui.lin@kaust.edu.sa; a.elzanaty@surrey.ac.uk; slim.alouini@kaust.edu.sa

% =========================================================================

% G. Lin, A. Elzanaty, and M.-S. Alouini, "LoRa Backscatter Communications: Temporal, Spectral, and Error Performance Analysis,"
% in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2023.3268113.

% =========================================================================
% Example: 
% snr_vec = 5:5:25; %snr vector in dB
% SF = 8; % spreading factor
% N = 4; % 2^N number of loads
% N_s = 5e5; % # of random symbols for MC simulation
% decoder = 'ML';
% decoder = 'fft';
% d = 10; % distance between Tx & Rx, does not represent the real distance
% ratio = 16; % ratio = d1/d2
% m = [10,2]; % m = [m1,m2], shape parameters of nakagami-m distribution
% ser = SER_WF_LB_MC(snr_vec,N_s,SF,N,decoder,d,ratio,m);
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function ser=SER_WF_LB_MC(snr_vec,N_s,SF,N,decoder,d,ratio,m)
%% Params
Es=1;
cnt=0;
cnt0=0; %# of time the transmitter does not sent energy
ser=zeros(1,length(snr_vec));
M=2^SF;
d1=d/(1+ratio);
d2=d-d1;
Omega1=Es/d1^2;
Omega2=Es/d2^2;
m1=m(1);
m2=m(2);
r1=m1/Omega1;r2=m2/Omega2;v=m1+m2;Eavg=1;n=m1-m2;
pd1 = makedist('Nakagami','mu',m1,'omega',Omega1);
pd2 = makedist('Nakagami','mu',m2,'omega',Omega2);
%% LB waveforms
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
x_conj=x_LB';
%% ML decoder
if strcmp(decoder,'ML')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    sigma_n=sqrt(Eavg./snr/(M*2)); 
    N0=2*sigma_n^2;
    %% solve nu0 numerically
    p_nu = @(nu) (nu>0).*2*(r1*r2/snr)^(v/2)/gamma(m1)/gamma(m2).*nu.^(v/2-1).*besselk(n,2*sqrt(r1*r2*nu/snr));

    F = @(nu_0) integral(@(nu) ...
    ((1 ./ nu_0) - (1 ./ nu)) .* p_nu(nu), nu_0, inf) - 1;
    % Initial guess. It is important parameter when snr is changed. Set
    % proper value of initial guess will produce the correct results.
    nu_0_guess = 1; 
    nu_0 = fzero(F, nu_0_guess);
    %% Loops for MC simulation
    parfor i=1:N_s
        a=randi([1,M]);
        h1=random(pd1);
        h2=random(pd2);
        h = h1*h2;
        noise=(randn(1,M)+1j*randn(1,M))*sigma_n;
        L_ML = (sqrt(Eavg*h^2/nu_0-M*N0)*x_LB(a,:)+noise)*x_conj;
        if (Eavg*h^2/nu_0-M*N0)<0
            cnt0=cnt0+1;
        else
            [~,ahat]=max(abs(L_ML));
            cnt=cnt+(a~=ahat);
        end
    end
    ser(j)=cnt/(N_s-cnt0);
    cnt=0;
    cnt0=0;
    
    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end
%% FFT decoder
elseif strcmp(decoder,'fft')
for j=1:length(snr_vec)
    snr=10.^(snr_vec(j)/10);
    sigma_n=sqrt(Eavg./snr/(M*2)); 
    N0=2*sigma_n^2;
    %% solve nu0 numerically
    p_nu = @(nu) (nu>0).*2*(r1*r2/snr)^(v/2)/gamma(m1)/gamma(m2).*nu.^(v/2-1).*besselk(n,2*sqrt(r1*r2*nu/snr));

    F = @(nu_0) integral(@(nu) ...
    ((1 ./ nu_0) - (1 ./ nu)) .* p_nu(nu), nu_0, inf) - 1;
    % Initial guess. It is important parameter when snr is changed. Set
    % proper value of initial guess will produce the correct results.
    nu_0_guess = 1; % Initial guess
    nu_0 = fzero(F, nu_0_guess);
    %% Loops for MC simulation
    parfor i=1:N_s
        a=randi([1,M]);
        h1=random(pd1);
        h2=random(pd2);
        h = h1*h2;
        fft_LB=decoder_output(SF,N,'fft');
        noise=(randn(1,M)+1j*randn(1,M))*sigma_n;
        L_fft=sqrt(Eavg*h^2/nu_0-M*N0)*fft_LB(a,:)+noise;
        if Eavg*h^2/nu_0-M*N0<0
            cnt0=cnt0+1;
        else
            [~,ahat]=max(abs(L_fft));
            cnt=cnt+(a~=ahat);
        end
    end
    ser(j)=cnt/(N_s-cnt0);
    cnt=0;
    cnt0=0;

    display=num2str(j/length(snr_vec)*100);
    disp(['process ',display,'%']);
end  
else
    disp('please input correct decoder')
end

end