clear all
close all
clc

%p=10;

%Mdl = arima('Constant',0.5,'AR',{0.7 0.25},'Variance',.1);

%rng(5)
%Y = simulate(Mdl,50);



%% nouveau

p=4;
N=5000;
nfft=1024;

f= -1/2:1/nfft:1/2-1/nfft;

x=randn(N,1);

R1=0,95;
R2=0.96;
theta1=pi/3;
theta2=pi/7;

p1 = R1*exp(-i*theta1);
p2 = R1*exp(i*theta1);

p3 = R2*exp(-i*theta2);
p4 = R2*exp(i*theta2);

B = [1, -(p1+p2+p3+p4), (p1+p2)*(p3+p4)+p1*p2+p3*p4, -(p1*p2*(p3+p4)+p3*p4*(p1+p2)), p1*p2*p3*p4];
A = [1];
signal_sortie=filter(1,B,x);

fft_signal= fftshift(abs(fft(signal_sortie,nfft)));
spectre_puissance= abs(fft_signal).^2/nfft;
w=2*pi*f;
H=freqz(A,B,w);
dsp=abs(H).^2; %*var^2

figure
plot(signal_sortie);
figure
plot(f, dsp);