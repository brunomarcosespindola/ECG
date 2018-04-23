
%% Tratamento entrada
clc; clear all; close all;
Name = '100m';
%Name = '121m';
%Name = '101m';
%Name = '16265m';


infoName = strcat(Name, '.info');
matName = strcat(Name, '.mat');
%Octave = exist('OCTAVE_VERSION');
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

%load('filtros.mat')
tam=length(val)/freqint(1,1); %calcula o tempo do sinal em segundos
t=0:tam/length(val):tam-(tam/length(val));% cria vetor de tempo para utilizar nos plots
% plot(t,val); % plot sinal original
% xlim([0 10])
%hold on
%y=filter(filtro,[val zeros(1,100)]);% filtragem das frequencias menores que 0.5Hz FIR ordem 200 
%y=y([101:end]);
%plot(t,y); % plot sinal filtrado
%y2=filter(filtro_IIR,y);%FILTRA 60Hz IIR ORDEM 16
%plot(t,y2); % plot sinal filtrado

%% DENOISE
%%%apply Wavelet Transform
% [C,L]=wavedec(y2,8,'db4');
% %[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
% % %%%Denoise 
% [thr,sorh,keepapp]=ddencmp('den','wv',y2);
% cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);
% plot(t,cleanecg) %plota sinal sem ruído
% figure()
% subplot(1,2,1)
% plot(t,val)
% xlim([0 4])
% title('sinal original')
% subplot(1,2,2)
% plot(t,cleanecg)
% xlim([0 4])
% title('sinal apos filtros')



%% Realiza os 3 niveis de decomposicoes-----------------------------

% [c,l] = wavedec(sumsin,3,'db2');
% approx = appcoef(c,l,'db2');
% [cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);

ordem=6; %deve ser colocado o valor correspondente a ordem da wavelet utilizada
%[a1,d1] = dwt(cleanecg,'db6');
[a1,d1] = dwt(val,'db6');
a1=a1([ordem:end]); %retira amostras que foram colocadas pela dwt
d1=d1([ordem:end]);
t1=0:tam/length(a1):tam-(tam/length(a1));

[a2,d2] = dwt(a1,'db6');
a2=a2([ordem:end]);
d2=d2([ordem:end]);
t2=0:tam/length(a2):tam-(tam/length(a2));

[a3,d3] = dwt(a2,'db6');
a3=a3([ordem:end]);
d3=d3([ordem:end]);
t3=0:tam/length(a3):tam-(tam/length(a3));

% [a4,d4] = dwt(a3,'db6');
% a4=a4([ordem:end]);
% d4=d4([ordem:end]);
% t4=0:tam/length(a4):tam-(tam/length(a4));

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
figure();
plot(t3,d3_2);
hold on
title('d3^2')

%% 



%% Encontra os picos
[peak_y1, peak_x1] = findpeaks(d3_2);
[m,n]=max(peak_y1);
[peak_y, peak_x] = findpeaks(d3_2,'minpeakheight',m*0.04,'MinPeakDistance',3);
taxa=(length(d3_2)/tam);
peak_x = (peak_x./taxa)-1/taxa;
plot(peak_x,peak_y,'x');

%% analise dos intervalos-------------------------------------
peak_aux = peak_x(2:end);
interval = peak_aux - peak_x(1:end-1);
bpm= (60./interval);
figure()
stem(bpm)

Media_bpm= sum(bpm)/length(bpm)
%figure()
%plot(peak_x,peak_y,'x');
%hold on
%plot(t,y2);

%%
%load('Hp100.mat');
load('ordem300.mat');

%y=filter(Hd,[a3 zeros(1,48)]);% filtragem das frequencias menores que 0.5Hz FIR ordem 200 
%y=y([49:end]);

%plot(t3,y);
%hold on

%plot(t3,a3);

%%apply Wavelet Transform

 [C,L]=wavedec(val,8,'db4');
% %[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
% % %%%Denoise 
[thr,sorh,keepapp]=ddencmp('den','wv',val);
cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);


figure()
plot (t,val);
hold on; grid on;
xlim([0 10]);
plot(t,cleanecg) %plota sinal sem ruído

y=filter(teste_ordem30,[cleanecg zeros(1,150)]);% filtragem das frequencias menores que 0.5Hz FIR ordem 200 
y=y([151:end]);
plot(t,y)
legend('original','Denoised','filtrado')
plotbrowser('on')
 %% 
 
 
% figure()
% subplot(1,2,1)
% plot(t,val)
% xlim([0 4])
% title('sinal original')
% subplot(1,2,2)
% plot(t,cleanecg)
% xlim([0 4])
% title('sinal apos filtros')

