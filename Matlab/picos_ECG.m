
%% Tratamento entrada
clc; clear all; close all;
%Name = '100m';
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


%% Filtragem inicial

load('filtros.mat');

time_signal_s=length(val)/freqint(1,1); %calcula o tempo do sinal em segundos
t=0:time_signal_s/length(val):time_signal_s-(time_signal_s/length(val));% cria vetor de tempo para utilizar nos plots
plot(t,val); % plot sinal original
xlim([0 30])
hold on
y=filter(FIR_HP_300,[val zeros(1,150)]);%filtragem das frequencias menores que 0.5Hz FIR ordem 300
y=y([151:end]);
plot(t,y); % plot sinal filtrado
%y2=filter(filtro_IIR,y);%FILTRA 60Hz IIR ORDEM 16
%plot(t,y2); % plot sinal filtrado
grid on;

%%FILTRO PB -- ELIMINA RUÍDO
y2 = filtfilt(SOS_PB,G_PB,y);
plot(t,y2)


%%DENOISE------------------------------
%%%apply Wavelet Transform
[C,L]=wavedec(y2,3,'sym5');
%[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
% %%%Denoise 
[thr,sorh,keepapp]=ddencmp('den','wv',y2);
cleanecg=wdencmp('gbl',C,L,'sym5',3,thr,sorh,keepapp);
%cleanecg=val;
plot(t,cleanecg,'b') %plota sinal sem ruï¿½do
legend('Original','Filtro Passa Altas','PB','Denoised')
plotbrowser('on');

figure()
subplot(2,1,1)
plot(t,y);
grid on
xlim([0 3])
xlabel('Segundos')
ylabel('mV')
title('Entrada do filtro passa-baixas')
subplot(2,1,2)
plot(t,y2)
grid on
xlim([0 3])
title('Saída do filtro passa-baixas')
xlabel('Segundos')
ylabel('mV')

%% Plot comparacao sinal filtrado com original----------------------------------------
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

coe_i = corrcoef(val,cleanecg);
correlacao=abs(coe_i); 
correlacao=correlacao(1,2)

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


d3_2=d3.^2; 


%% Encontra os picos R do d3^2

m=max(d3_2);
[peak_y, peak_x] = findpeaks(d3_2,t3,'minpeakheight',m*0.01,'MinPeakDistance',0.3);

% analise dos intervalos-------------------------------------
peak_aux = peak_x(2:end);
interval = peak_aux - peak_x(1:end-1);
bpm= (60./interval);

Media_bpm= sum(bpm)/length(bpm)



%% encontra picos R sinal PREPROCESSADO
m2=max(cleanecg);
[peak_y3, peak_x3] = findpeaks(cleanecg,t,'minpeakheight',m2*0.5,'MinPeakDistance',0.3);


% analise dos intervalos picos sinal original-------------------------------------
peak_aux2 = peak_x3(2:end);%%CORRIGIR
interval2 = peak_aux2 - peak_x3(1:end-1);
bpm2= (60./interval2);

Media_bpm_original= sum(bpm2)/length(bpm2)

%% Encontra onda Q e S do sinal PREPROCESSADO

t_matrix=repmat(t,length(peak_x3),1); %cria uma matriz de tamanho: length(peak_x3) x length(t) apenas repetindo o t
peak_x3_matrix= repmat(peak_x3',1,length(t)); %cria uma matriz de tamanho: length(peak_x3') x length(t) apenas repetindo o peak_x3'

S=sum(t_matrix>peak_x3_matrix & t_matrix<(peak_x3_matrix+0.3)); %localiza todas as amostras que estão no intervalo de 0,3s após cada onda R
Q=sum(t_matrix<peak_x3_matrix & t_matrix>(peak_x3_matrix-0.3));%localiza todas as amostras que estão no intervalo de 0,3s antes de cada onda R


pos_S=find(S); %pega a posição das amostras
pos_Q=find(Q); %pega a posição das amostras


[Sy,Sx]  = findpeaks(- cleanecg(pos_S), t(pos_S),'MinPeakDistance',0.5); %encontra os picos S
[Qy,Qx]  = findpeaks(- cleanecg(pos_Q), t(pos_Q),'MinPeakDistance',0.5); %encontra os picos Q

Q_S=sum(Sx([1:length(Qx)])-Qx)/length(Qx);

string = ['intervalo QS= ',num2str(Q_S*1000),' ms'];disp(string);


%% Plots comparando 


figure();
subplot(2,1,1)
plot(peak_x3,peak_y3,'x');
xlim([0 10])
hold on;
plot(Sx,-Sy,'vg')%% marca onda S
plot(Qx,-Qy,'vb')%% marca onda Q
plot(t,cleanecg)
title('Complexo QRS no sinal pré-processado')
grid on;
hold off;
subplot(2,1,2)
plot(t3,d3_2);
%xlim([0 30])
hold on
title('d3^2')
plot(peak_x,peak_y,'x');
grid on;





%% encontra picos R sinal a1
m_a1=max(a1);
[peak_a1y, peak_a1x] = findpeaks(a1,t1,'minpeakheight',m_a1*0.5,'MinPeakDistance',0.3);


% analise dos intervalos picos sinal original-------------------------------------
peak_aux_a1 = peak_a1x(2:end);
interval_a1 = peak_aux_a1 - peak_a1x(1:end-1);
bpm_a1= (60./interval_a1);

Media_bpm_a1= sum(bpm_a1)/length(bpm_a1)


figure()
subplot(3,1,2)
stem(peak_x(2:end),bpm)
xlabel('segundos')
ylabel('BPM')
grid on;
title('Método da decomposição wavelet')
subplot(3,1,1)
stem(peak_x3(2:end),bpm2)
xlabel('segundos')
ylabel('BPM')
title('Método Direto')
grid on;
subplot(3,1,3)
stem(peak_a1x(2:end),bpm_a1)
xlabel('segundos')
ylabel('BPM')
title('calculado pelo a1')
grid on;


%% Encontra onda Q e S do sinal a1
t1_matriz=repmat(t1,length(peak_a1x),1);
peak_a1x_matriz= repmat(peak_a1x',1,length(t1));

S_a1=sum(t1_matriz>peak_a1x_matriz & t1_matriz<(peak_a1x_matriz+0.1));
Q_a1=sum(t1_matriz<peak_a1x_matriz & t1_matriz>(peak_a1x_matriz-0.1));


pos_S_a1=find(S_a1); 
pos_Q_a1=find(Q_a1);


[Sy_a1,Sx_a1]  = findpeaks(- a1(pos_S_a1), t1(pos_S_a1),'MinPeakDistance',0.5);
[Qy_a1,Qx_a1]  = findpeaks(- a1(pos_Q_a1), t1(pos_Q_a1),'MinPeakDistance',0.5);

Q_S_a1=sum(Sx_a1([1:length(Qx_a1)])-Qx_a1)/length(Qx_a1);

string = ['intervalo QS pelo a1= ',num2str(Q_S_a1*1000),' ms'];disp(string);


figure()
plot(t1,a1)
hold on
plot(t,cleanecg)
plot(peak_a1x,peak_a1y,'x')
title('a1');
plot(Sx_a1,-Sy_a1,'vg')%% marca onda S
plot(Qx_a1,-Qy_a1,'vb')%% marca onda Q
hold off
legend('a1','Denoised','R','S','Q')
plotbrowser('on');




