% *************** Simulation BER of BPSK transmission over noise *********%
% with using Monte Carlo simulation
% Mahmoud Aldababsa
% ********************* Initialization ***********************************%
clc;
clear all;
close all;
loop=100; % Monte Carlo 
M=10000;  % Frame length (x_1 x_2 ... x_M)
SNRdB=0:30; % SNR in dB
SNR=10.^(SNRdB./10);
Rate=zeros(1,31);

load('PDP_samp.mat');
%for i=1:39
 
% ********************* Transmitter **************************************%
for dB= 1: length(SNRdB) % start looping by SNR
    dB;
    for lp= 1: loop, % start looping of frame data 
	
% ********************* BPSK signal generation ***************************%    
	x_inp=sign(rand(M,1)- 0.5); % 1 for inphase and -1 for outphase
    N0=1./SNR(dB);
    sigma(dB)=sqrt(N0/2);
    noise=sigma(dB)*randn(M,1);
% ********************* Channel ******************************************%
	  
    %y_channel=awgn(x_inp,SNRdB(dB)); % Additive White Guassiann Noise (AWGN) 
    %load('TDL_taps.mat');
    
    x_chan = conv(x_inp,PDP_samp);
    x_chan= x_chan(1:10000,1);
    y_channel=x_chan+noise;
    
% ********************* Receiver *****************************************% 
    y=y_channel;
    x_out= sign(real(y)); % 
   
% ********************* Bit Error Rate (BER) calulation ******************%    
    
    [err, rate]= symerr(x_inp, x_out);
    Rate(dB)= Rate(dB) + rate;
    
    end % end for loop
   
    Rate(dB)= Rate(dB)/loop; % Average value over Monte Carlo simulation 
                              % loop
   
end % end Monte Carlo
%end
% ********************* Plot the simulation result ***********************%

    f1 = figure(1);
    set(f1,'color',[1 1 1]);
    semilogy(SNRdB,Rate, 'b')
    hold on;
    
    BER_th=(1/2)*erfc(sqrt(SNR)); % theoritical calculation for BER
    semilogy(SNRdB,BER_th,'r');
    
    hold on;
    axis([0 30 0.000001  1.2]);  
    xlabel( 'Signal-to-Noise Ratio (SNR)')
    ylabel( 'Bit Error Rate (BER)')
    title('Simulation BPSK transmission over noise');
    legend('BER measured','BER theoritical')
    grid on;

    % ********************* END ******************************************%