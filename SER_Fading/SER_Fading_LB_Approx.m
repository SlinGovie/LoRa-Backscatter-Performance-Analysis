% =========================================================================
% -- Script to compute average symbol error rate (SER) for
% LoRa backscatter by GL quadrature. The weights and nodes can be fastly
% computed by executing "gauss_laguerre_weights_nodes.jl" file in the
% directory. A "gl_nodes_weights.mat" will be generated after execution and
% this script will read the .mat file and use the nodes and weights for
% computation.
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
% N_GH = 20;
% N_GL = 100;
% ser = SER_Fading_LB_Approx(snr_vec,SF,N,decoder,d,ratio,N_GH,N_GL)
% figure;
% semilogy(snr_vec,ser)
% grid on
% =========================================================================
function ser=SER_Fading_LB_Approx(snr_vec,SF,N,decoder,d,ratio,N_GH,N_GL)
%% Params
Es=1;
% snr_vec=snrmin:snrstep:snrmax;
M=2^SF;
% d2=(100/d1^2)^0.5;
d1=d/(1+ratio);
d2=d-d1;
Omega1=Es/d1^2;
Omega2=Es/d2^2;
m1=10;
m2=2;

r1=m1/Omega1;
r2=m2/Omega2;
v=m1+m2; % v is the exponent \alpha in generalized Gauss-Laguerre quadrature
n=m1-m2;

%% generate nodes and weights for Gauss Hermite quadrature
% N_GH=20;
syms x k h
p1=hermiteH(N_GH,x);
xt=roots(sym2poly(p1));
wt=(2^(N_GH-1)*factorial(N_GH)*sqrt(pi))/N_GH^2./hermiteH(N_GH-1,xt).^2;%wt

%% load the weights and nodes of GL quadrature. 

% Read the MAT file
% run gauss_laguerre_weights_nodes.jl wto obtain the data
data = load('gl_nodes_weights.mat');

xs = data.x;
ws = data.w;

ser=zeros(1,length(snr_vec));

%% ML decoder
if strcmp(decoder, 'ML') || strcmp(decoder, 'fft')
    ML=decoder_output(SF,N,'ML');
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
            atil(q)=fd(q);
    end

    sumn= symsum(factorial(n+k)/factorial(k)/factorial(n-k)/(2*h).^k,k,0,n);
    sumn_= symsum(factorial(n-1+k)/factorial(k)/factorial(n-1-k)/(2*h).^k,k,0,n-1);
    sumn =matlabFunction(sumn);
    sumn_=matlabFunction(sumn_);

    parfor j=1:length(snr_vec)
        ML=decoder_output(SF,N,'ML');
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
        summation=0;
        SEP_a=zeros(1,length(atil));
    
        %---------reclaim value for parfor
        tmp=snr_vec(j);
        snr1=10.^(tmp/10);
        sigma_n=sqrt(1./snr1/(M*2)); 
            for a=1:length(atil)
                kappa_neg=@(h)-(h*abs(ML(atil(a),atil(a)))).^2/(2*sigma_n^2);%-v^2/2sigma_n^2
                mu_g=@(h) sigma_n*sqrt(pi/2)*laguerreL(1/2,kappa_neg(h));%mu_g
                sigma_g=@(h)sqrt(2*sigma_n^2+(h*abs(ML(atil(a),atil(a)))).^2-mu_g(h).^2);
                for t=1:length(xt)
                    lhatcdff=@(h) 1-lhatcdfawgn(sqrt(2)*sigma_g(h/(2*sqrt(r1*r2)))*xt(t)+mu_g(h/(2*sqrt(r1*r2))),h/(2*sqrt(r1*r2)),ML(atil(a),:),sigma_n);
                    for s=1:length(xs)
                         summation=summation+wt(t)*ws(s)*2^(3/2-v)/(gamma(m1)*gamma(m2))*lhatcdff(xs(s))*sqrt(sumn(xs(s))*sumn_(xs(s))/xs(s))*(sqrt(n-1/2)*gamma(n)/gamma(n+1/2));
        
                    end
                end
                SEP_a(a)=summation
                summation=0;
            end
        ser(j)=1/M*(SEP_a*Na)
    end

else
    disp('please input correct decoder')
end

end

function flhat = lhatcdfawgn(x,h,ml,sigma_n)
    S=sort(uniquetol(abs(ml)));
    V=S(1:end-1);
    Nv=zeros(length(V),1);
    flhat= 1;
    for i=1:length(V)
        Nv(i)=length(find(abs(abs(ml)-V(i))<1e-5));
        cdfi = 1-marcumq(h*abs(V(i))/sigma_n,x.*(x>=0)/sigma_n);
        cdfi= cdfi.^Nv(i);
        flhat= flhat.*cdfi;
    end
end

