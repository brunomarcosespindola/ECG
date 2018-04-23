
%% Tratamento entrada
clc; clear all; close all;
%Name = '100m';
%Name = '121m';
%Name = '101m';
Name = '16265m';


infoName = strcat(Name, '.info');
matName = strcat(Name, '.mat');
Octave = exist('OCTAVE_VERSION');
load(matName);
fid = fopen(infoName, 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[freqint] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = 0.001;
fgetl(fid);

for i = 1:size(val, 1)
    [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
end

fclose(fid);
val(val==-32768) = NaN;
fval = 1./size(val, 1);

for i = 1:size(val, 1)
    val(i, :) = (val(i, :) - base(i)) / gain(i);
end


%% 

load('filtros.mat');

time_signal_s=length(val)/freqint(1,1); %calcula o tempo do sinal em segundos
t=0:time_signal_s/length(val):time_signal_s-(time_signal_s/length(val));% cria vetor de tempo para utilizar nos plots
plot(t,val); % plot sinal original
xlim([0 10])
hold on
y=filter(FIR_HP_300,[val zeros(1,150)]);%filtragem das frequencias menores que 0.5Hz FIR ordem 200
y=y([151:end]);
plot(t,y); % plot sinal filtrado
y2=filter(filtro_IIR,y);%FILTRA 60Hz IIR ORDEM 16
plot(t,y2); % plot sinal filtrado
grid on;


%% DENOISE
%%%apply Wavelet Transform
[C,L]=wavedec(y2,8,'db4');
%[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
% %%%Denoise 
[thr,sorh,keepapp]=ddencmp('den','wv',y2);
cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);
plot(t,cleanecg) %plota sinal sem ruído
legend('Original','Filtro Passa Altas','Filtro 60Hz','Denoised')
plotbrowser('on');

figure()
subplot(1,2,1)
plot(t,val)
xlim([0 4])
title('sinal original')
grid on;
subplot(1,2,2)
plot(t,cleanecg)
xlim([0 4])
title('sinal apos filtros')
grid on;



%% Realiza os 3 niveis de decomposicoes-----------------------------

% [c,l] = wavedec(sumsin,3,'db2');
% approx = appcoef(c,l,'db2');
% [cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
wavelet='db7';
ordem=7; %deve ser colocado o valor correspondente a ordem da wavelet utilizada
[a1,d1] = dwt(cleanecg,wavelet);
a1=a1([ordem:end]); %retira amostras que foram colocadas pela dwt
d1=d1([ordem:end]);
t1=0:time_signal_s/length(a1):time_signal_s-(time_signal_s/length(a1));

[a2,d2] = dwt(a1,wavelet);
a2=a2([ordem:end]);
d2=d2([ordem:end]);
t2=0:time_signal_s/length(a2):time_signal_s-(time_signal_s/length(a2));

[a3,d3] = dwt(a2,wavelet);
a3=a3([ordem:end]);
d3=d3([ordem:end]);
t3=0:time_signal_s/length(a3):time_signal_s-(time_signal_s/length(a3));

% [a4,d4] = dwt(a3,'db6');
% a4=a4([ordem:end]);
% d4=d4([ordem:end]);
% t4=0:time_signal_s/length(a4):time_signal_s-(time_signal_s/length(a4));

%----------------------------------------------------------------


%----------Plots das decomposicoes------------------------------
%%
janela=5; %define o intervalo dos prints
figure()
subplot(3,2,1);
plot(t1,a1);
title('a1')
xlim([0 janela])
subplot(3,2,2);
plot(t1,d1);
title('d1')
xlim([0 janela])
xlabel('Seconds')
ylabel('Amplitude')

subplot(3,2,3);
plot(t2,a2);
title('a2')
xlim([0 janela])
subplot(3,2,4);
plot(t2,d2);
title('d2')
xlim([0 janela])

subplot(3,2,5);
plot(t3,a3);
title('a3')
xlim([0 janela])
subplot(3,2,6);
plot(t3,d3);
title('d3')
xlim([0 janela])


d3_2=d3.^2; %plota o d3^2


%% Encontra os picos
[peak_y1, peak_x1] = findpeaks(d3_2);
[m,n]=max(peak_y1);
[peak_y, peak_x] = findpeaks(d3_2,'minpeakheight',m*0.03,'MinPeakDistance',3);
taxa=(length(d3_2)/time_signal_s);
peak_x = (peak_x./taxa)-1/taxa;



%% analise dos intervalos-------------------------------------
peak_aux = peak_x(2:end);
interval = peak_aux - peak_x(1:end-1);
bpm= (60./interval);

Media_bpm= sum(bpm)/length(bpm)
%figure()
%plot(peak_x,peak_y,'x');
%hold on
%plot(t,y2);



%% encontra picos sinal original

[peak_y2, peak_x2] = findpeaks(cleanecg);
[m2,n2]=max(peak_y2);
[peak_y3, peak_x3] = findpeaks(cleanecg,'minpeakheight',m2*0.5,'MinPeakDistance',3);
taxa2=(length(cleanecg)/time_signal_s);
peak_x3 = (peak_x3./taxa2)-1/taxa2;


%% analise dos intervalos picos sinal original-------------------------------------
peak_aux2 = peak_x3(2:end);
interval2 = peak_aux2 - peak_x3(1:end-1);
bpm2= (60./interval2);

Media_bpm_original= sum(bpm2)/length(bpm2)

%% Plots comparando 
figure()
subplot(2,1,1)
stem(bpm)
xlabel('amostras')
ylabel('BPM')
grid on;
title('Utilizando wavelet')
subplot(2,1,2)
stem(bpm2)
xlabel('amostras')
ylabel('BPM')
title('Direto')
grid on;

figure();
subplot(2,1,1)
plot(peak_x3,peak_y3,'x');
hold on;
plot(t,cleanecg)
title('original')
grid on;
hold off;
subplot(2,1,2)
plot(t3,d3_2);
hold on
title('d3^2')
plot(peak_x,peak_y,'x');
grid on;



