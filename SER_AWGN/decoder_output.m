function decoder = decoder_output(SF,N,type)
if SF<6 || SF>12
    disp('please check the value of spreading factor')
end

if N<2 
    disp('please check the bits of quantizer')
end 
E_s=1;
M=2^SF;
phi=zeros(M,M);
y=zeros(M,M);
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
x=(E_s/M)^0.5*exp(1j*pi*Q);
x_conj=conj(x).';
ML=x*x_conj;
for i=1:M
    y(i,:)=fft(x(i,:).*downchirp,M);
end

if  strcmp(type,'ML')
decoder=ML;
elseif  strcmp(type,'fft')
decoder=y;
else 
        disp('please input correct decoder')
end


end
