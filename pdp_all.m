data=xlsread('190524-PHD_LAB-CESA-KONF1-CAL_SlotAnt.xlsx');
f=56:0.1:64;
d=40;

for x=2:40
    m=1;
for j= 4:84
Pr(m)= data(j,x);
m=m+1;
end
for i=1:length(Pr)
    PrAbs(i)=10^(Pr(i)/10);
end
% Amplitude 
for i=1:length(Pr)
    Ar(i)=sqrt(PrAbs(i));
end
%calculation for delay. Delay=2pi*f*d/c...Assuming d=40m (Distance between
%Tx and RX)
for i=1:length(Pr)
    del(i)=2*pi*f(i)*d/(3*10^8);
end
for i=1:length(Pr)
    
    ampi(i)=abs(Ar(i).*cos(del(i)));
end

%Final channel magnitude and phase plot

%k    = hilbert(log(ampi));



Hf_sim     = hilbert((ampi));                                          
Hf_sim_mag = sqrt(real(Hf_sim).^2+imag(Hf_sim).^2);      
H_sim_a    = -angle(Hf_sim);


f= 56:0.1:64;

%Frequency vector
freqlow  = 56;                         %Start frequency
frequp   = 64;                         %Stop frequency
step     = 0.1;                         %Step size
N        = (frequp-freqlow)/step+1;      %Number of points
freq     = freqlow:step:frequp;          %Frequency vector
BW       = frequp-freqlow;               %Overall bandwidth

%PDP
tres=1/BW;                 %Time resolution
timeP=0:tres:(N-1)*tres;   %Time axis for PDP
htm=ifft(Hf_sim,N);            %Measured CIR
PDPm=20*log10(abs(htm));   %Measured PDP

hts      = flipud(ifft(Hf_sim,N));       %Simulated CIR through Hilbert
PDPs     = -20*log10(abs(hts));           %Simulated PDP through Hilbert

figure
subplot(3,1,1)
stem(hts); title('IFFT of H(f)'); xlabel('Time');ylabel('h(t)');
grid on;
subplot(3,1,2)
stem(PDPs);title('Power Delay Profile'); xlabel('Time');ylabel('Power');
grid on;
subplot(3,1,3)
plot(PDPs);title('Power Delay Profile'); xlabel('Time');ylabel('Power');
grid on;
%{
r=num2str(x-1);
Fig_name= strcat('figures\Figure',r);
savefig(Fig_name);
%}

end
