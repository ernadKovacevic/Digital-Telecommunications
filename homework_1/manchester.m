clear all
close all
clc

M=2;
Rb=50000;
Rs=Rb/log2(M);
Tb=1/Rb;
fc=1e6;
fs=2e6;
Ac=5;
fill=fs/Rs;
Nbits=1e6;
bits=round(rand(1,Nbits)); %GENERISANJE BITA
%%%% LINIJSKI KODER %%%%
man_out=[];
for i=1:length(bits)
   if bits(i)==1 
       man_out=[man_out ones(1,fill/2)*Ac ones(1,fill/2)*(-Ac)];
   else
       man_out=[man_out ones(1,fill/2)*(-Ac) ones(1,fill/2)*Ac];
   end
end

figure(1)
subplot(2,1,1)
stem(bits(1:10));
title('Generisani biti');
grid on;
subplot(2,1,2)
plot(man_out(1:400),'r');
title('Izlaz iz Manchester kodera');
grid on;

%%%% KANAL %%%%
A=1.95;
B=59;
delta=0.03;
l=2;
f=0:100:fs;
gama=(1-j*delta/2)*(A*sqrt(f/fc)+j*B*(f/fc));
H=exp(-gama*l);

channel=fir2(10,f/fs,abs(H));
channel_out=conv(channel,man_out);
channel_out=channel_out(6:end-5);

figure(2)
plot(channel_out(1:400))
title('Izlaz iz kanala');
grid on;

%%%% AWGN %%%%

awgn_sig=zeros(5,length(channel_out));
awgn_sig(1,:)=awgn(channel_out,1);
awgn_sig(2,:)=awgn(channel_out,5);
awgn_sig(3,:)=awgn(channel_out,10);
awgn_sig(4,:)=awgn(channel_out,15);
awgn_sig(5,:)=awgn(channel_out,20);

figure(3)
subplot(5,1,1)
plot(awgn_sig(1,1:100))
title('AWGN 1dB');
grid on;
subplot(5,1,2)
plot(awgn_sig(2,1:100))
title('AWGN 5dB');
grid on;
subplot(5,1,3)
plot(awgn_sig(1,1:100))
title('AWGN 10dB');
grid on;
subplot(5,1,4)
plot(awgn_sig(1,1:100))
title('AWGN 15dB');
grid on;
subplot(5,1,5)
plot(awgn_sig(1,1:100))
title('AWGN 20dB');
grid on;

%%%% OPTIMALNI FILTER %%%%

h1=awgn_sig(1,end:-1:1);
h2=awgn_sig(2,end:-1:1);
h3=awgn_sig(3,end:-1:1);
h4=awgn_sig(4,end:-1:1);
h5=awgn_sig(5,end:-1:1);

optimal_filter=zeros(5,length(awgn_sig));
optimal_filter(1,:)=filter(h1,1,awgn_sig(1,:));
optimal_filter(2,:)=filter(h2,1,awgn_sig(2,:));
optimal_filter(3,:)=filter(h3,1,awgn_sig(3,:));
optimal_filter(4,:)=filter(h4,1,awgn_sig(4,:));
optimal_filter(5,:)=filter(h5,1,awgn_sig(5,:));

figure(4)
subplot(5,1,1)
plot(optimal_filter(1,1:400));
title('Izlaz optimalnog filtera za signal sa sumom 1dB');
grid on;
subplot(5,1,2)
plot(optimal_filter(2,1:400));
title('Izlaz optimalnog filtera za signal sa sumom 5dB');
grid on;
subplot(5,1,3)
plot(optimal_filter(3,1:400));
title('Izlaz optimalnog filtera za signal sa sumom 10dB');
grid on;
subplot(5,1,4)
plot(optimal_filter(4,1:400));
title('Izlaz optimalnog filtera za signal sa sumom 15dB');
grid on;
subplot(5,1,5)
plot(optimal_filter(5,1:400));
title('Izlaz optimalnog filtera za signal sa sumom 20dB');
grid on;

%%%% DETEKTOR %%%%
detector=[];
temp=[];
for j=1:size(optimal_filter,1)
    for i=fill/4:fill:length(optimal_filter)
         if optimal_filter(j,i) > 0
             temp=[temp 5];
         else
             temp=[temp -5];
         end
    end
    detector=[detector; temp];
    temp=[];
end
figure(5)
subplot(5,1,1)
stem(detector(1,1:10));
title('Izlaz detektora za signal sa sumom od 1dB');
grid on;
subplot(5,1,2)
stem(detector(2,1:10));
title('Izlaz detektora za signal sa sumom od 5dB');
grid on;
subplot(5,1,3)
stem(detector(3,1:10));
title('Izlaz detektora za signal sa sumom od 10dB');
grid on;
subplot(5,1,4)
stem(detector(4,1:10));
title('Izlaz detektora za signal sa sumom od 15dB');
grid on;
subplot(5,1,5)
stem(detector(5,1:10));
title('Izlaz detektora za signal sa sumom od 20dB');
grid on;

%%%% BER %%%%

temp_detector= detector >0;
N=zeros(1,size(detector,1));
Pe=zeros(1,size(detector,1));
for i=1:size(detector,1)
    [N(i) Pe(i)]=biterr(bits,temp_detector(i,:));
end


