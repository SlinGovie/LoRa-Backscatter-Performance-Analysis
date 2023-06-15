% =========================================================================
% -- Script to compute power spectral density of LoRa baskcatter signals,
% with both Monte Carlo simulation and theoretical expressions
% =========================================================================

% -- (c) 2023 Ganghui Lin, Ahmed Elzanaty, Mohamed-Slim Alouini

% -- e-mail: ganghui.lin@kaust.edu.sa; a.elzanaty@surrey.ac.uk; slim.alouini@kaust.edu.sa

% =========================================================================

% G. Lin, A. Elzanaty, and M.-S. Alouini, "LoRa Backscatter Communications: Temporal, Spectral, and Error Performance Analysis,"
% in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2023.3268113.

% =========================================================================
%% LB PSD Monte Carlo vs theoretical
clc
clear
close all
%% Monte Carlo Simulation
%% Modulation parameters
% Modify NQ and SFv to obtain PSD for different SF and N
for NQ=[2]
    Breal=250e3;
    fc=868.3e6;
    B=1;
    oversamp=2^5;
    fs=oversamp*B;
    Tsamp=1/fs;
    SFv=[7];
    Q=@(x)(floor(x*2^(NQ-1))+0.5)/2^(NQ-1);
    for sfindex=1:length(SFv)
        SF=SFv(sfindex);
        disp(['SF= ',num2str(SF)]);
        M=2^SF;
        Ts=M/B;
        time=0:Tsamp:Ts-Tsamp;
        %     N= 2^nextpow2(max(1e3,ceil(10*M*log(M))));%  # of random symbols
        N=1e3; % # of random symbols
        alpha=randi(M,1,N)-1; % Modulating number takes random value between 0 and L-1 (representing SF bits)
        NsampSymbol=fs*Ts; % # of samples in each time symbol
        LoRacomplexenvelop=zeros(1,floor(N*NsampSymbol));
        %% Signal description in time domain
%         phita=@(t,a) 2*pi*B*t.*(a/M-0.5+t.*(0.5/M)-(t>=(M-a)/B)) ; % LoRa
        phita=@(t,a) pi*Q(2*B*t.*(a/M-0.5+t.*(0.5/M)-(t>=(M-a)/B))) ; % LB
        %% Data Generation
        for j=1:N
            LoRacomplexenvelop(1+(j-1)*NsampSymbol:j*NsampSymbol)=exp(1i*phita(time,alpha(j)));
        end
        %% Monte Carlo Spectrum
        IRequiredFreqResolution=B/M/10;
        fftsize=2^nextpow2(1/(IRequiredFreqResolution*Tsamp));
        [psdMC,f_c] = pwelch(LoRacomplexenvelop,hanning(fftsize),[],fftsize,fs,'centered','psd');
        FreqResolution=f_c(2)-f_c(1);
        Bend=8;
        freqindexrandge=find(f_c>=0&f_c<=Bend);
        freqrange=f_c(freqindexrandge)/B;
        psdrange=10*log10(psdMC(freqindexrandge)*B);
        %% Plot
        figure;
        plot(freqrange,psdrange,'o','MarkerIndices',1:5:length(psdrange));
        plot(freqrange,psdrange,'o')
        ylabel('Magnitude (dB)')
        xlabel('Freq/B')
        chartitile=['NQ=',num2str(NQ)];
        title(chartitile)
        title('Monte Carlo PSD')
        hold on
        grid on
        %% Save
%         PSMCData=[freqrange.';psdrange.'].';
%         save(['PSMonteCarloDataSF',num2str(SF),'NQ',num2str(NQ)],'PSMCData')

        %% Theoretical
        df = 0.001; % frequency resolution = 0.001 * Bandwidth
        f_c=0:df:Bend;
        l=0:Bend*M;
        f_d=l*B/M;
        Sa_c=zeros(M,length(f_c));
        Sa_d=zeros(M,length(f_d));

        for a=0:M-1 % traverse all symbols
            %% obtain t_m
            ft=@(t) 2*B*t.*(a/M-0.5+B*t.*(0.5/M)-(t>=(M-a)/B));
            i=ceil(-(M-2*a)^2/(2^(3-NQ)*M)):floor((a*(M-a))/(M*2^(1-NQ)));
            rtpos=((M-2*a)+sqrt((M-2*a)^2+2^(3-NQ)*i*M))/(2*B);
            rtneg=((M-2*a)-sqrt((M-2*a)^2+2^(3-NQ)*i*M))/(2*B);
            rtall=[rtpos rtneg];
            for k=1:length(rtall)
                if rtall(k)<0
                   rtall(k)=rtall(k)+M/B;
                end
            end
            rtsort=[sort(unique(rtall)),Ts];
        
            tm=rtsort(1:end-1);
            tmp=rtsort(2:end);
            tmid=(tm+tmp)/2;
            %% Fourier transform of LB waveforms
            for j=1:length(tm)
                Sa_c(a+1,:)= Sa_c(a+1,:)+1j./(2*pi*f_c) *exp(1j*pi*Q(ft(tmid(j)))).*(exp(-1j*2*pi*f_c*tmp(j))-exp(-1j*2*pi*f_c*tm(j)));
                Sa_d(a+1,:)= Sa_d(a+1,:)+1j./(2*pi*f_d) *exp(1j*pi*Q(ft(tmid(j)))).*(exp(-1j*2*pi*f_d*tmp(j))-exp(-1j*2*pi*f_d*tm(j)));
            end
        
        end
        %% Theoretical PSD Gc->continuous, Gd->discrete
        G_c= 1/(M*Ts)*(sum(abs(Sa_c).^2)-1/M*abs(sum(Sa_c)).^2);
        G_d= abs(sum(Sa_d)).^2/(Ts*M)^2;



        %% Passband Spectrum
        %% Add continous and discrete parts together
        G=zeros(1,length(f_c));
        G(1)=G_c(1);
        for n=2:length(f_c)
            sumGd=0;
            A=find(f_d>f_c(n-1)& f_d<f_c(n));
            for k=1:length(A)
                sumGd=sumGd+G_d(A(k));
            end
            G(n)=(G_c(n)*df+sumGd)/df;
        end
        Gpas1=10*log10(G)+10*log10(0.0502)+10*log10(2); 
        Gpas2=flip(Gpas1);
        Gpas=[Gpas2 Gpas1(2:end)];
        fpas=-Bend:df:Bend;
        fpas=fc+fpas*Breal;
        
        %% Creat spectrum mask
        fmask=8.676e8:0.000004e8:8.69e8;
        mask=zeros(1,length(fmask));
        for i=1:length(fmask)
            if i<=500 && i>=1
                mask(i)=-36;
            end
            if i<=1000 && i>=501
                mask(i)=-30;  
            end
            if i<=2500 && i>=1001
                mask(i)=0; 
            end
            if i<=3000 && i>=2501
                mask(i)=-30;  
            end
            if i<=3501 && i>=3001
                mask(i)=-36;  
            end      
        end

        %% Plot continuous and discrete parts
        G_c = 10*log10(G_c*B);
        G_d = 10*log10(G_d*B);

        figure;
        subplot(2, 1, 1);
        plot(f_c, G_c);
        grid on;
        ylabel('dB');
        xlabel('f/B');
        title('Continuous');
        
        subplot(2, 1, 2);
        plot(f_d, G_d,'.');
        grid on;
        ylabel('dB');
        xlabel('f/B');
        title('Discrete');
        %% Plot passband spectrum with mask
        figure;
        plot(fpas,Gpas)
        xlim([867.6e6 869e6])
        ylim([-60 0])
        grid on
        hold on
        plot(fmask,mask)

        %% Save
        %    PSDConData=[f;Gc].';
        %    save(['PSDConDataSF',num2str(SF),'NQ',num2str(NQ)],'PSDConData')
        %    PSDDisData=[fd;Gd].';
        %    save(['PSDDisDataSF',num2str(SF),'NQ',num2str(NQ)],'PSDDisData')
        %    PSData=[f_all;G].';
        %    save(['PSDataSF',num2str(SF),'NQ',num2str(NQ)],'PSData')
    end
end