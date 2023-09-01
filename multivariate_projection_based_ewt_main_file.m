clc;
clear all;
close all;

%%%%code released for multivariate projection based EWT filter bank
%%%if you are interested to use this code, you need to cite the following
%%%papers
%1. Tripathy, Rajesh Kumar, et al. "Development of automated sleep stage classification system 
%using multivariate projection-based fixed boundary empirical 
%wavelet transform and entropy features extracted from multichannel EEG signals." Entropy 22.10 (2020): 1141.
%2.Gilles, Jerome. "Empirical wavelet transform." IEEE transactions on signal processing 61.16 (2013): 3999-4010.


load sample01.mat;
val=epo.x(1:2500,:,:);
Fs=epo.fs; % fs=250 Hz
L=size(val,2);

xx=val(:,:,5); %%%%multisensor data: samples x channels
x=xx'; %%%convert it to channel x samples
for i=1:size(x,1)
    x(i,:)=x(i,:)/max(abs(x(i,:))); %%%signal normalization
end

%%%%%%projection of the multichannel data into the direction cosine%%%%
Combined = (sum(x))/sqrt(60); %%%%projected signal
nfft=250;
%%%%%multivariate projection EWT
Z=abs(fftshift(fft(Combined, nfft)));
Z1=Z(Fs/2+1:end); %%%FFT domain signal

Nb=15;
bound = LocalMax(Z1,Nb); %%%evaluation of local maximas

boundaries=(bound*(2*pi))/Fs; %%%evaluation of boundary points
ff=fft(Combined);
% We build the corresponding filter bank
mfb=EWT_Meyer_FilterBank(boundaries,length(ff)); %%filter bank design on projected signal
Bound=1;
xxx=(linspace(0,1,round(length(mfb{1,1}))))*Fs;

for i=1:size(mfb)
plot(xxx,mfb{i,1})
hold on
end
xlim([0 75])
ylim([0 2])
title('Multivariate Projection based EWT filter bank')

for ch=1:size(x,1)
    all_ch_modes(ch,:,:)=mode_eval_each_channel(x(ch,:),mfb);
end

ch1modes=all_ch_modes(1,:,:);
ch1modes=reshape(ch1modes, [Nb, size(x,2)]);
%%%mode plot
figure
plot(ch1modes(1,:))
