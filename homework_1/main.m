clear all
close all
clc
M=2;
Rb=50000;
Rs=Rb/log2(M);
Tb=1/Rb;
Rs=Rb;
delta=0.03;
l=2000;
fc=1e6;
fs=2e6;
Ac=5;
fill=fs/Rs;
Nbits=10^6;;
Nplot=5;
bits=round(rand(1,Nbits)); %GENERISANJE BITA
%bits=[0 1 1 1 0 1 0 0 1 0]
NRZ_out=[];
man_out=[];

%%% LINIJSKI KODER %%%
for i=1:length(bits)
   if bits(i)==1 
       NRZ_out=[NRZ_out ones(1,fill)*Ac];
       man_out=[man_out ones(1,fill/2)*Ac ones(1,fill/2)*(-Ac)];
   else
       NRZ_out=[NRZ_out ones(1,fill)*(-Ac)];
       man_out=[man_out ones(1,fill/2)*(-Ac) ones(1,fill/2)*(Ac)];
   end
end

figure(1)
subplot(2,1,1)
stem(bits);
title('Generisani biti');
grid on;
subplot(2,1,2)
plot(NRZ_out,'r');
title('NRZ kod');
grid on;

figure(2)
subplot(2,1,1)
stem(bits);
title('Generisani biti');
grid on;
subplot(2,1,2)
plot(man_out,'r');
title('Manchester kod');
grid on;
signal=[];

%%% KANAL %%%
% A=1.95;
% B=59;
% gama=(1-j*delta/2)*(A*sqrt(f/fc)+j*B*f/fc);
% H=exp(-gama*l);
% NRZ_out=fir2(9,0:100:fs,
% man_out
%%AWGN%%%
SNR=[1 5 10 15 20];
nrz_awgn=[];
man_awgn=[];
for i=1:length(SNR)
    nrz_awgn=[nrz_awgn; awgn(NRZ_out,SNR(i))];
    man_awgn=[man_awgn; awgn(man_out,SNR(i))];
end

% j=1;
% for i=3:2+length(SNR)
%     figure(i)
%     subplot(2,1,1)
%     plot(nrz_awgn(j,1:size(nrz_awgn,2)));
%     title(['NRZ kod za SNR = ' num2str(SNR(j))]);
%     grid on;
%     subplot(2,1,2)
%     plot(man_awgn(j,1:size(man_awgn,2)));
%     title(['Manchester kod za SNR = ' num2str(SNR(j))]);
%     grid on;
%     j=j+1;
% end

%%OPTIMALNI FILTER%%%

NRZ_out=[];
man_out=[];
for i=1:size(nrz_awgn,1)
    h1=nrz_awgn(i,end:-1:1);
    h2=man_awgn(i,end:-1:1);
    NRZ_out=[NRZ_out; filter(h1,1,nrz_awgn(i,:))];
    man_out=[man_out; filter(h2,1,man_awgn(i,:))];
end
figure(3)
plot(NRZ_out(3,:));

%%%DETEKTOR%%%
opt_out=[];
optman_out=[];
for i=fill/2:fill:size(NRZ_out,2)
    if NRZ_out(3,i) > 0
        opt_out=[opt_out 1];
    else
        opt_out=[opt_out 0];

    end
end

figure(10)
stem(opt_out);
grid on;
title('DETEKTOR');

N0=1;
Pe=0;
eb=linspace(0,25,50);
for i=1:length(eb)
    Eb=N0 * 10^(eb(i)/10);
    Pe(i)=1/2*erfc(sqrt(Eb/N0));
end
figure(11)
plot(eb,Pe)
[a ber]=biterr(bits,opt_out);
    
    

