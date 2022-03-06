clear all
close all
clc

Ac=0.1;
Rb=1e6;
Tb=1/Rb;
fc=2/Tb;
fs=4*fc;
Tb=1/Rb;
Nbits=1e5;
fill=fs/Rb;

%======== TX ========
bits=round(rand(1,Nbits));
figure (1)
stem(bits(1:10));
title('Slucajno generisani biti');
ylabel('Amplituda (V)');
xlabel('t(s)');

%======== BPSK MODULACIJA ========
NRZ_out=[];
for i=1:length(bits)
    if bits(i)==1
        NRZ_out=[NRZ_out ones(1,fill)*Ac*sqrt(Tb/2)];
    else
        NRZ_out=[NRZ_out ones(1,fill)*-Ac*sqrt(Tb/2)];
    end
end

figure (2)
subplot(3,1,1)
plot(NRZ_out(1:10*fill));
title('Izlaz iz NRZ kodera');
ylabel('Amplituda (V)');
xlabel('t(s)');

t=0:1/fs:length(bits)*Tb-1/fs;
s1=sqrt(2/Tb)*cos(2*pi*fc*t);

subplot(3,1,2)
plot(s1(1:10*fill));
title('Bazna BPSK funkcija');
ylabel('Amplituda (V)');
xlabel('t(s)');

BPSK_out=NRZ_out.*s1;
subplot(3,1,3)
plot(BPSK_out(1:10*fill));
title('BPSK modulisani signal');  
xlabel('t(s)');
ylabel('Amplituda (V)');

%======== KANAL ========
h=[-0.015 0.058 -0.35 1 -0.35 0.058 -0.005];
channel_out=conv(h,BPSK_out);
channel_out=channel_out(4:end-3);
figure (3)
plot(channel_out(1:10*fill))
title('Signal na izlazu iz kanala');
xlabel('t(s)');
ylabel('Amplituda (V)');

%======== AWGN ======== 
SNR=[1 5 10 15 20 25];
awgn_out=awgn(channel_out,SNR(6));
awgn_out(1,:)=awgn(channel_out,SNR(1));
awgn_out(2,:)=awgn(channel_out,SNR(2));
awgn_out(3,:)=awgn(channel_out,SNR(3));
awgn_out(4,:)=awgn(channel_out,SNR(4));
awgn_out(5,:)=awgn(channel_out,SNR(5));
awgn_out(6,:)=awgn(channel_out,SNR(6));

figure (4)
subplot(3,1,1)
plot(awgn_out(1,1:10*fill))
title('SNR 1dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,2)
plot(awgn_out(2,1:10*fill))
title('SNR 5dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,3)
plot(awgn_out(3,1:10*fill))
title('SNR 10dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

figure (5)
subplot(3,1,1)
plot(awgn_out(4,1:10*fill))
title('SNR 15dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,2)
plot(awgn_out(5,1:10*fill))
title('SNR 20dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,3)
plot(awgn_out(6,1:10*fill))
title('SNR 25dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

%======== EKVALIZATOR TRECEG REDA (METOD FORSIRANJA NULA) ========
x=[1 -0.35 0.058;   
   -0.35 1 -0.35;
   0.058 -0.35 1];
z=[0 1 0]';
c=inv(x)*z;

koef=zeros(1,length(h));
x=[0 h 0];
         
for i=1:length(koef)
    koef(i)=x(i+2)*c(1) + x(i+1)*c(2) + x(i)*c(3);   
end
y_zfe=zeros(length(SNR),length(awgn_out));

for i=1:length(SNR)
    equal=conv(koef,awgn_out(i,:));
    y_zfe(i,:)=equal(4:end-3);
end

figure (6)
subplot(3,1,1)
plot(y_zfe(1,1:10*fill))
title('Ekvaliziran signal(ZFE) za SNR 1dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,2)
plot(y_zfe(2,1:10*fill))
title('Ekvaliziran signal(ZFE) za SNR 5dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,3)
plot(y_zfe(3,1:10*fill))
title('Ekvaliziran signal(ZFE) za SNR 10dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

figure (7)
subplot(3,1,1)
plot(y_zfe(4,1:10*fill))
title('Ekvaliziran signal(ZFE) za SNR 15dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,2)
plot(y_zfe(5,1:10*fill))
title('Ekvaliziran signal(ZFE) za SNR 20dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,3)
plot(y_zfe(6,1:10*fill))
title('Ekvaliziran signal(ZFE) za SNR 25dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

%======== EKVALIZATOR TRECEG REDA (METOD MINIMALNE SREDNJE KVADRATNE GRESKE) ========

x=[-0.35 0.058 -0.015;
   1 -0.35 0.058;
   -0.35 1 -0.35;
   0.058 -0.35 1;
   -0.05 0.056 -0.35];

z=[0 0 1 0 0]';
Rxx_inv = inv ( x' * x);
Rxz=x' * z;
c=Rxx_inv * Rxz;

koef=zeros(1,length(h));
x=[0 h 0];
         
for i=1:length(koef)
    koef(i)=x(i+2)*c(1) + x(i+1)*c(2) + x(i)*c(3);   
end

y_mmse=zeros(length(SNR),length(awgn_out));

for i=1:length(SNR)
    equal=conv(koef,awgn_out(i,:));
    y_mmse(i,:)=equal(4:end-3);
end

figure (8)
subplot(3,1,1)
plot(y_mmse(1,1:10*fill))
title('Ekvaliziran signal(MMSE) za SNR 1dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,2)
plot(y_mmse(2,1:10*fill))
title('Ekvaliziran signal(MMSE) za SNR 5dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,3)
plot(y_mmse(3,1:10*fill))
title('Ekvaliziran signal(MMSE) za SNR 10dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

figure (9)
subplot(3,1,1)
plot(y_mmse(4,1:10*fill))
title('Ekvaliziran signal(MMSE) za SNR 15dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,2)
plot(y_mmse(5,1:10*fill))
title('Ekvaliziran signal(MMSE) za SNR 20dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

subplot(3,1,3)
plot(y_mmse(6,1:10*fill))
title('Ekvaliziran signal(MMSE) za SNR 25dB');
xlabel('t(s)');
ylabel('Amplituda (V)');

%======== BPSK DEMODULACIJA ========
bpsk_dem=[];
bpsk_dem_mmse=[];
bpsk_dem_zfe=[];
for i=1:length(SNR)
    bpsk_dem=[bpsk_dem; awgn_out(i,:) .*s1];
    bpsk_dem_mmse=[bpsk_dem_mmse; y_mmse(i,:) .*s1];
    bpsk_dem_zfe=[bpsk_dem_zfe; y_zfe(i,:) .*s1];
end

j=1;
x=zeros(1,length(bits));
x_mmse=zeros(1,length(bits));
x_zfe=zeros(1,length(bits));
for k=1:length(SNR)
    for i=1:fill:length(bpsk_dem)
        if sum(bpsk_dem(k,i:fill*j)) >= 0
            x(k,j)=1;
        else
            x(k,j)=0;
        end
        
        if sum(bpsk_dem_zfe(k,i:fill*j)) >= 0
            x_zfe(k,j)=1;
        else
            x_zfe(k,j)=0;
        end
        
        if sum(bpsk_dem_mmse(k,i:fill*j)) >= 0
            x_mmse(k,j)=1;
        else
            x_mmse(k,j)=0;
        end 
        j=j+1;
    end
    j=1;
end

%======== RX ========  
simb=zeros(1,length(SNR));
simb_zfe=zeros(1,length(SNR));
simb_mmse=zeros(1,length(SNR));
Pe=zeros(1,length(SNR));
Pe_mmse=zeros(1,length(SNR));
Pe_zfe=zeros(1,length(SNR));

for i=1:length(SNR)
    [simb(i), Pe(i)]=biterr(x(i,:),bits);
    [simb_mmse(i), Pe_mmse(i)]=biterr(x_mmse(i,:),bits);
    [simb_zfe(i), Pe_zfe(i)]=biterr(x_zfe(i,:),bits);
end
Pe
Pe_mmse
Pe_zfe

figure (10)
plot(SNR,Pe,'k-o');
title('Vjerovatnoca greske bez ekvalizatora');
xlabel('SNR');
ylabel('BER');
grid on

figure(11)
plot(SNR,Pe_mmse,'k-o');
title('Vjerovatnoca greske sa ekvalizatorom (metod sr. kvadratne greske)');
xlabel('SNR');
ylabel('BER');
grid on

figure(12)
plot(SNR,Pe_zfe,'k-o');
title('Vjerovatnoca greske sa ekvalizatorom (metod forsiranja nula)');
xlabel('SNR');
ylabel('BER');
grid on

figure(13)
plot(SNR,Pe,'k-o');
legend('AWGN','MMSE','ZFE')
hold on
title('Poredjenje vjerovatnoce greske bez ekvalizera,MMSE,ZFE');
xlabel('SNR');
ylabel('BER');

plot(SNR,Pe_mmse,'r-o');
legend('AWGN','MMSE','ZFE')
hold on
xlabel('SNR');
ylabel('BER');

plot(SNR,Pe_zfe,'b-o');
legend('AWGN','MMSE','ZFE')
grid on

