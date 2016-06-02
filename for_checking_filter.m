clc, clear all, format compact, close all

dat1= importdata('Data\S1\Trial1\EEGData.mat');
dat=dat1(1:end-1,2:end-1);

fsi=1000;
Wni=[0.4 4];
n_fi=3;
len=300;
data=dat(:,1);
tmax=length(data)/fsi;
t=0:1/fsi:tmax; t1=t(1:end-1);

fnq= 1/(2*(1/fsi));  %nyquist frequency 
Wn = [Wni(1)/fnq Wni(2)/fnq];
[b,a] = butter(n_fi,Wn);
[h,t] = impz(b,a,len,fsi);
plot(t,h)

fdata=fftfilt(h,data);

myfit=fit(t1',data,'poly1');
fmyfit=fit(t1',fdata,'poly1');


% a test signal
drift=0:1/10:(1/10)*length(t1);
yy=0*0.5+sin(2*pi*t1)+cos(2*pi*t1*0.01);
y=yy+0.1*randn(1,length(t1))+cumsum(drift(1:end-1))*4e-7;

fy=fftfilt(h,y);

myf=fit(t1',y','poly1');
fmyf=fit(t1',fy','poly1');



figure,
subplot 211
plot(t1,data,'k'), hold on, plot(t1,fdata,'g'), grid
subplot 212
line([0 t1(end)],[myfit.p2 myfit.p1*t1(end)+myfit.p2],'color',[0 0 0],'linewidth',2), hold on
line([0 t1(end)],[fmyfit.p2 fmyfit.p1*t1(end)+fmyfit.p2],'color',[0 1 0],'linewidth',2), grid

figure,
subplot 211
plot(t1,y,'k'), hold on,plot(t1,yy,'b'), plot(t1,fy,'g'), grid
subplot 212
line([0 t1(end)],[myf.p2 myf.p1*t1(end)+myf.p2],'color',[0 0 0],'linewidth',2), hold on
line([0 t1(end)],[fmyf.p2 fmyf.p1*t1(end)+fmyf.p2],'color',[0 1 0],'linewidth',2), grid

figure
freqz(b,a); xlim([-0.01 0.04])


