%% This script generates Power Angle Profile (PAP)(disabled) from PDP data.
% This also plots the peaks of the PDP from T_max
%% On data : 190524-PHD_LAB-CESA-KONF1-CAL_SlotAnt.csv
%% Dated : 05.06
%% Clearing all

%clear all;
close all;
%% System settings
% System constants
PlotType = 'IS';  % 'IS' ImageSc
CCFit = 'CIR';      % type of fitting

%% Numerical constants
AZnum = 13;      % [30,35,0,+-5,+-10,+-15,+-20,+-25]
ELnum = 3;        % [+-5,0]
Blim = [56,64];  % lower and upper cutoff frequency
Slices = (1:39);  % AMB - user defined CTFs
Elevation = 3;    %AMB
IntRat = 8;  %Interpolation ratio = 2^IntRat

%% Read file of measured data
[filename, path] = uigetfile('*.csv','Enter file name',' ');
[Freq, magSdBV, EL, AZ] = getSparCSV(strcat(path, filename)); 
ELstep = EL(2) - EL(1);
AZstep = AZ(2) - AZ(1);
dimdBV = size(magSdBV);

%% Frequency span reduction

% lower coordinate for freq. Blim
if (Blim(1)>Freq(1))
    NpL = (Blim(1) - Freq(1)) * 10 + 1;
else
    NpL = 1;
end

% upper coordinate for freq. Blim
if (Blim(2)<Freq(end))
    NpH = numel(Freq) - (Freq(end) - Blim(2))*10;
else
    NpH = numel(Freq);
end

% new frequency and data vectors
FScale = 10^9;
Freq = Freq(NpL:NpH,1)*FScale;
magSdBV = magSdBV(NpL:NpH,:);
freqL = Freq(1);
freqU = Freq(end);
Fs = (freqU-freqL)/(numel(Freq)-1);   %Frequency step in GHz

% FREQUENCY DOMAIN PROCESSING
%% convert measured data matrix to required shape
magSdB = zeros(81,3,13);
SizeMtrx = [13 3];  %[AZ EL]
for i = 1:numel(Freq)
    magSdB(i,:,:) = reshape (magSdBV(i,:), SizeMtrx)'; %magSdB is not used elsewhere
end


Np = NpH - NpL + 1;

%% Hilbert Transform
% phase calculation from antenna distance
H_Rx = 13.5;
H_Tx = 1.6;
D = 107;
%dAnt = sqrt(D^2 + (H_Rx-H_Tx)^2);
dAnt = 107.66;
TauCIR = dAnt/3e8;
magSV = 10.^(magSdBV./20);
SizemagSdBV = size(magSdBV);

% HT over all frequencies for particular AZ,EL
HphaseSRawV = zeros(SizemagSdBV(1),SizemagSdBV(2));
HphaseSV = zeros(SizemagSdBV(1),SizemagSdBV(2));
HSCmpxV = zeros(SizemagSdBV(1),SizemagSdBV(2));
for i = 1:SizemagSdBV(2)
    HphaseSRawV(:,i) = hilbtran(log(magSV(:,i))); % phase with HT only
    HphaseSV(:,i) = HphaseSRawV(:,i) + 2*pi*(Freq - Freq(1))*TauCIR; %phase with HT + Phase Comp.
    %HSCmpxV(:,i) = (magSV(:,i).*exp(1j*HphaseSRawV(:,i))); %Frequency response with HT phase
    HSCmpxV(:,i) = (magSV(:,i).*exp(1j*HphaseSV(:,i))); %Frequency response with HT + Phase Comp.
end


%HSCmpxV = Hf_sim;
%% Plots for phase comparison
 freq_points = 56:0.1:64;
 subplot(2,2,1);
 plot(freq_points,HphaseSRawV(:,19));
 xlabel('Frequency');
 ylabel('Phase');
 title('Phase Recovery using Hilbert Transform');

subplot(2,2,2);
plot(freq_points,angle(HSCmpxV(:,19)));
xlabel('Frequency');
ylabel('Phase');
 title('Phase Recovery using HT + Phase Compensation');
 
subplot(2,2,3);
plot(freq_points,HphaseSRawV(:,28));
 xlabel('Frequency');
 ylabel('Phase');
 title('Phase Recovery using Hilbert Transform');
 
 subplot(2,2,4);
 plot(freq_points,HphaseSV(:,28));
 xlabel('Frequency');
 ylabel('Phase');
 title('Phase Recovery using HT + Phase Compensation');
 


%% IFFT

win(:,1) = ones(Np,1);
win(:,2) = hann (Np, 'periodic');
win(:,3) = hamming(Np, 'periodic');
win(:,4) = blackman(Np, 'periodic');
win(:,5) = flattopwin(Np,'periodic');

wintype = 3; %hamming windowing selected

HTCmpxV = zeros(SizemagSdBV(1),SizemagSdBV(2));
HmagTV = zeros(SizemagSdBV(1),SizemagSdBV(2));
for i=1:39
    HTCmpxV(:,i) = ifft (HSCmpxV(:,i).*win(:,wintype), Np, 1);
    HmagTV(:,i) = abs(HTCmpxV(:,i));
end

HmagdBT = 20*log10(HmagTV);    %magnitude log in time domain


%% RMS Delay Spread

SVar = magSV(:, Slices);    

SVarSize = size(SVar);

Tstop = 70;
BW = freqU - freqL;                 % bandwidth
Tr = abs(1./BW);
time = 0:Tr:Tr*(Np-1);
PDPsel = HmagTV(:, Slices).^2;       % Power Delay profile        
PDPselPks = zeros(SVarSize);      % Prepare array for MPC PDP
AvgTau = zeros(1, SVarSize(2));
RMSDelSpr = zeros(1, SVarSize(2));

 figure;
 title('PDP LOS onwards');
xlabel('Time[ns]');
 ylabel('PDP[dBm]');
 legend_check=1;
 grid on;
 hold on;
for i = 1:SVarSize(2)
    % Find peaks of PDP and their locations
    [PksPDPsel, PksLocs] = findpeaks(PDPsel(:, i));
    
    % Find start time
    [MaxPk, PosPk] = max(PksPDPsel);    % Find max peak and its position
    Tstart = PksLocs(PosPk);
    TimeEval = Tstart:Tstop;            % Tstart:Tstop or 1:SVarSize(1);
    Tau = time(TimeEval);               % TO UNDERSTAND
    %plot(PksLocs(PosPk:size(PksLocs)),10*log(PksPDPsel(PosPk:size(PksPDPsel))));
%     if legend_check==1
%         legend(split(num2str(Slices)));
%         legend_check=0;
%     end
    % Set peaks in MPC PDP
    PDPselPks(PksLocs, i) = PksPDPsel;
    PDPselPks([1:Tstart - 1, Tstop + 1:SVarSize(1)], i) = 0;  % set zeros outside TimeEval
    
    % RMS Delay Spread
    PDP = PDPselPks;                    % PDPselPks or PDPsel
    AvgTau(i) = sum(PDP(TimeEval, i).* Tau')/sum(PDP(TimeEval, i));
    RMSDelSpr(i) = sqrt(sum(PDP(TimeEval, i).*(Tau - AvgTau(i))'.^2)/sum(PDP(TimeEval, i)));
end
PDPseldB = 10*log10(PDPsel);
PDPselPksdB = 10*log10(PDPselPks);


 indx=1;
 for i=1:3
     for j=1:13
         RMSmat(i,j) = RMSDelSpr(indx);
         indx = indx + 1; 
    end
 end
 

for i=1:3
std_azim_5_to5(i)=std(RMSmat(i,:));
 end
 
 for i=1:13
 std_elev_25_to35(i)=std(RMSmat(:,i));
 end


%% Rearranging PDP


%% PDP Plotting -  Multiple
Plotting
 figure;
 plot(Freq, 20*log(abs(HSCmpxV))); %|H(f)| vs f
 xlabel('Frequency[Hz]');
 ylabel('Power [dBm]');
 

 
  figure;
 subplot(2,2,1);
  plot(time*FScale,PDPseldB(:,19));
  xlabel('Propagation Delay (\tau_n)[ns] \rightarrow');
  ylabel('Power[dBm] \rightarrow');
  title('Power Delay Profile-19');
  grid on;
  hold on;
  subplot(2,2,2);
  plot(time*FScale,PDPseldB(:,15));
  xlabel('Propagation Delay (\tau_n)[ns] \rightarrow');
  ylabel('Power[dBm] \rightarrow');
  title('Power Delay Profile-15');
  grid on;
 subplot(2,2,3);
 plot(time*FScale,PDPseldB(:,22));
  xlabel('Propagation Delay (\tau_n)[ns] \rightarrow');
  ylabel('Power[dBm] \rightarrow');
  title('Power Delay Profile-22');
  grid on;
  subplot(2,2,4);
  plot(time*FScale,PDPseldB(:,35));
  xlabel('Propagation Delay (\tau_n)[ns] \rightarrow');
  ylabel('Power[dBm] \rightarrow');
  title('Power Delay Profile-35');
  grid on;
hold on
  stem(PDPselPksdB(:,21), 'filled', 'k',  'MarkerSize', 3, 'LineStyle', 'none', 'Color','red');
  sampled_pks = zeros(81,1);
  for i=28:5:81
      sampled_pks(i)=PDPseldB(i,21);
  end
  
  %stem(sampled_pks,'filled','k', 'MarkerSize',3, 'LineStyle', 'none', 'Color','green');
  %legend('PDP','Actual peaks','Sampled points');
   %xlim([0,85])
  %ylim([min(min(PDPseldB))-5, max(max(PDPseldB))+5])
  
%  %plot(time*FScale, PDPselPksdB)
%  legend(split(num2str(Slices)));
  %ylim([min(min(PDPseldB)), -220])
 
%  hold off

%  clf;
%  plot(time*FScale, PDPseldB)
%  legend(split(num2str(Slices)));
%  xlabel('Time [ns]');
%  ylabel('Magnitude [dB]');
%  grid on
%  hold on
%  stem((time*FScale)', PDPselPksdB, 'filled', 'k',  'MarkerSize', 3, 'LineStyle', 'none');
%  %plot(time*FScale, PDPselPksdB)
%  legend(split(num2str(Slices)));
  
%  hold off

%% PDP Ploting - Single
%{
figure;
lot(time*FScale,PDPseldB(:,19));
xlabel('Propagation Delay (\tau_n) [ns] \rightarrow');
ylabel('Power [dBm] \rightarrow');
title('Power Delay Profile-19');
grid on;
%}
%% Formatting data for Power Angle Profile
% sumPDPsel = sum(PDPsel);
% avgPDPsel = sumPDPsel ./ 81;
% avgPDPseldB = 10 * log(avgPDPsel);
% avgPDPseldB = avgPDPseldB';
% avgPDPseldB_col1 = avgPDPseldB(1:13,1);
% avgPDPseldB_col2 = avgPDPseldB(14:26,1);
% avgPDPseldB_col3 = avgPDPseldB(27:39,1);
% avgPDPseldBmat = [avgPDPseldB_col1 avgPDPseldB_col2 avgPDPseldB_col3];
% 
% thetaDeg = [-25;-20;-15;-10;-5;0;5;10;15;20;25;30;35];
% theta = pi/180 .* thetaDeg;
% 
% figure;
% for i=1:3
%     subplot(2,2,i);
%     polarplot(theta,avgPDPseldBmat(:,i));
%     thetalim([-90 90]);
%     rmin = min(avgPDPseldBmat(:,i));
%     rmax = max(avgPDPseldBmat(:,i));
%     rlim([rmin rmax]);
%     if i==1
%         title('Elevation = +5\circ');
%     elseif i==2
%         title('Elevation = 0\circ');
%     else
%         title('Elevation = -5\circ');
%     end
% end


%% PAP curve fitting ( polar -> rect -> polar)
%[xrec,yrec] = pol2cart(theta, avgPDPseldBmat(:,1));
%plot(xrec,yrec);
%hold on;
%[xrec,yrec] = pol2cart(theta,avgPDPseldBmat(:,2));
%plot(xrec,yrec);
%[xrec,yrec] = pol2cart(theta,avgPDPseldBmat(:,3));
%plot(xrec,yrec);
% hold on;
% coeff = polyfit(xrec,yrec,4);
% y_pre = coeff(1)*(xrec.^4) + coeff(2)*(xrec.^3) + coeff(3)*(xrec.^2) + coeff(4)*(xrec.^1) + coeff(5);
% plot(xrec,y_pre);

%% Channel Parameters using param_dscrt_PDP.m
% RMSdelay = 1.1234 ns (generated value)
% [meanDelay, rmsDelay, symbolRate, coherenceBW] = param_dscrt_PDP(PDPseldB(:,19),time);
% fprintf('Mean Delay: %f (ns)\n',meanDelay*1e9);
% fprintf('RMS Delay: %f (ns)\n',rmsDelay*1e9);
% fprintf('Symbol Rate: %f (Gigabits perSec)\n',symbolRate*1e-9);
% fprintf('Coherence Bandwidth: %f (Hz)\n',coherenceBW*1e-9);
%% TDL Model
 
% N=10; %No of taps
% c=21;  % Combination number (azim,elev)
% intvl= 5; %Regular interval for sampling
% start_pos=28; % LOS components of all combinations are at 28 time units
% end_pos= start_pos+(N*intvl);
% 
% tau=intvl:intvl:(N*intvl); %Start after LOS component
% 
% Sample_pts= transpose(PDPsel(start_pos+intvl:intvl:end_pos,c)); %Sampling PDP at regular interval of 'intvl'
% 
% %Random coeffs (not sure)
% samp= tau.^2;
% random_coeffs=(1/pi)*exp(-(samp*(10^-18))); %Circularly  symmetric gaussian function 
% 
% tap_gain=(Sample_pts.*random_coeffs); %TDL model
% tap_gain_m=tap_gain/PDPsel(start_pos); %TDL normalized with LOS=1;
% 
% tap_gain_m1 = zeros(1,11);
% for i=1:length(tap_gain_m)
%     tap_gain_m1((i+1))= tap_gain_m(i);
% end
% tap_gain_m1(1)=1;

Plot TDL
figure
stem(tau,tap_gain_m);

%% Revised TDL Model

PDP_avg= mean(PDPsel,2);
 TDL_norm = zeros(1,10);
 TDL_norm_dB = zeros(1,10);
 seq = 10;
 sample_pts = 10:4:10+9*4;
 
 tseq=39;
 %for i=1:39

 
 TDL_norm = PDP_avg(sample_pts)./PDP_avg(sample_pts(1));
 TDL_norm_dB = 10*log10(TDL_norm);

 %TDL_norm(i,:) = PDPsel(sample_pts,i)./PDPsel(sample_pts(1),i);
 %end
 
 
 %{
 for i=1:size(sample_pts,2)
     TDL_norm(i) = PDPseldB(sample_pts(1),seq)/PDPseldB(sample_pts(i),seq);
     TDL_norm_dB(i) = abs(PDPseldB(sample_pts(1),seq)) - abs(PDPseldB(sample_pts(i),seq));
 end
  %}

 %{
 figure;
 stem((sample_pts-10)*Tr*1e9,TDL_norm);
 xlabel('Delay (ns)');
 ylabel('Normalised Tap value');
 title('Average TDL');
%}
 
 PDP_samp = TDL_norm;
 save PDP_samp.mat PDP_samp
 
%% Bit Error Rate Calculation

% BER

% num_bit=1e6;%number of bit
% dataIn=randi([0 1],1,num_bit);%random bit generation (1 or 0)data
% dataFrmtIn=2*dataIn-1;%conversion of data for BPSK modulation s
% SNRdB=0:2.5:50; % SNR in dB
% SNR=10.^(SNRdB/10);
% 
% Eb=1;
% 
% R=50e6;
% tp= (1/R)*10e9;
% 
% dataInSamp= dataIn(1:tp:num_bit); %sampleTx dataInSampled
% dataFrmtInSamp= dataFrmtIn(1:tp:num_bit); %sampleddata dataFrmtInSamp
% 
% for k=1:length(SNRdB)%BER (error/bit) calculation for different SNR
% No=Eb/SNR(k); 
%     
%  y1= conv(dataFrmtInSamp,TDL_norm); %Convolution excl. LOS tap
%  
%  N=sqrt(No/2)*randn(1,length(y1));   %Generate AWGN
%          
%  y=y1+N;                           %Received Signal
% 
% error=0;
% ber_result = zeros(1,21);
% for i=1:1:length(dataInSamp)
%     if(y(i)<0)
%         y(i)=0;
%     else
%         y(i)=1;
%     end
% end
% [numErrors,ber2] = biterr(dataInSamp,y(1:size(dataInSamp,2)));
% %error=error/num_bit;%Calculate error/bit
% %m(k)=error;
% ber = numErrors/num_bit;
% ber_result(k)=ber;
% end
% figure 
% %plot start
% plot(SNRdB,ber_result),grid on,hold on;
% 
% %xlabel(' SNR(dB)');
% %ylabel('BER');


numBits = 1e6;
dataIn = randi([0 1],1, numBits);
hModulator = comm.BPSKModulator;
hDemodulator = comm.BPSKDemodulator;

xMod = hModulator(dataIn');
y = conv(TDL_norm,xMod);
yDemod = hDemodulator(y);
yDemod = yDemod(1:numBits,1);
[numErr ber] = biterr(dataIn',yDemod);
