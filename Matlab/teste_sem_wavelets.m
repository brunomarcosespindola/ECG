
%% Tratamento entrada
clc; clear all; close all;
Name = '100m'; %arritmia, onda T negativa
%Name = '16265m';% sinusal
%Name ='16272m'; %sinusal
%Name ='16420m'; %sinusal 


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

if (freqint(1,1)<200)  %realiza interpolação caso a frequencia de amostragem do sinal original seja baixa demais para a decomposição  
    r=3;
    desloc=2.4;
    sinal_original= interp(val,r);
    
else
    r=1;
    desloc=7;
    sinal_original= val;
end


time_signal_s=length(sinal_original)/(freqint(1,1)*r); %calcula o tempo do sinal em segundos
t=0:time_signal_s/length(sinal_original):time_signal_s-(time_signal_s/length(sinal_original));% cria vetor de tempo para utilizar nos plots


y=filter(FIR_HP_300,[sinal_original zeros(1,150)]);%filtragem das frequencias menores que 0.5Hz FIR ordem 300
y=y([151:end]);

%%FILTRO PB -- ELIMINA RUÍDO
y2 = filtfilt(SOS_PB,G_PB,y);

%% Plot comparando os sinais iniciais
figure();
plot (t,sinal_original)
hold on
plot(t,y)
plot(t,y2)
grid on;
title('Estágios de Pre-processamento')
plotbrowser('on');
%%
% plot(t,sinal_original); % plot sinal original
% xlim([0 30])
% hold on
% plot(t,y); % plot sinal filtrado
% grid on;
% plotbrowser('on');

% %% Realiza os 3 niveis de decomposicoes-----------------------------
% 
% % [c,l] = wavedec(sumsin,3,'db2');
% % approx = appcoef(c,l,'db2');
% % [cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
% wavelet='sym6';
% ordem=6; %deve ser colocado o valor correspondente a ordem da wavelet utilizada
% [a1,d1] = dwt(y2,wavelet);
% a1=a1([ordem:end]); %retira amostras que foram colocadas pela dwt
% d1=d1([ordem:end]);
% t1=0:time_signal_s/length(a1):time_signal_s-(time_signal_s/length(a1));
% 
% [a2,d2] = dwt(a1,wavelet);
% a2=a2([ordem:end]);
% d2=d2([ordem:end]);
% t2=0:time_signal_s/length(a2):time_signal_s-(time_signal_s/length(a2));
% 
% [a3,d3] = dwt(a2,wavelet);
% a3=a3([ordem:end]);
% d3=d3([ordem:end]);
% t3=0:time_signal_s/length(a3):time_signal_s-(time_signal_s/length(a3));

% [a4,d4] = dwt(a3,'db6');
% a4=a4([ordem:end]);
% d4=d4([ordem:end]);
% t4=0:time_signal_s/length(a4):time_signal_s-(time_signal_s/length(a4));

%----------------------------------------------------------------


% %% ----------Plots das decomposicoes------------------------------
% %
% janela=5; %define o intervalo dos prints
% figure()
% subplot(3,2,1);
% plot(t1,a1);
% title('a1')
% xlim([0 janela])
% grid on;
% subplot(3,2,2);
% plot(t1,d1);
% title('d1')
% xlim([0 janela])
% grid on;
% 
% xlabel('Seconds')
% ylabel('Amplitude')
% 
% subplot(3,2,3);
% plot(t2,a2);
% title('a2')
% xlim([0 janela])
% grid on;
% subplot(3,2,4);
% plot(t2,d2);
% title('d2')
% xlim([0 janela])
% grid on;
% 
% subplot(3,2,5);
% plot(t3,a3);
% title('a3')
% xlim([0 janela])
% grid on;
% subplot(3,2,6);
% plot(t3,d3);
% title('d3')
% xlim([0 janela])
% grid on;



%% encontra picos R sinal y2
m2=max(y2);
[peak_y, peak_x] = findpeaks(y2,t,'minpeakheight',m2*0.3,'MinPeakDistance',0.3);

% analise dos intervalos RR do a3_atenuado-------------------------------------
peak_aux2 = peak_x(2:end);%%CORRIGIR
intervalo = peak_aux2 - peak_x(1:end-1);
figure()
plot(peak_x(2:end),intervalo);
title('intervalos de tempo RR')
xlabel('Segundos')
ylabel('RR (seg)')
grid on;
bpm= (60./intervalo);

Media_bpm= sum(bpm)/length(bpm)

%% Encontra onda Q e S do sinal PREPROCESSADO

t_matriz=repmat(t,length(peak_x),1); %cria uma matriz de tamanho: length(peak_x) x length(t3) apenas repetindo o t3
peak_x_matrix= repmat(peak_x',1,length(t)); %cria uma matriz de tamanho: length(peak_x3') x length(t) apenas repetindo o peak_x3'

S=sum(t_matriz>peak_x_matrix & t_matriz<(peak_x_matrix+0.15)); %localiza todas as amostras que estão no intervalo de 0,15s após cada onda R
Q=sum(t_matriz<peak_x_matrix & t_matriz>(peak_x_matrix-0.15));%localiza todas as amostras que estão no intervalo de 0,15s antes de cada onda R


pos_S=find(S); %pega a posição das amostras
pos_Q=find(Q); %pega a posição das amostras


[Sy,Sx]  = findpeaks(- y2(pos_S), t(pos_S),'MinPeakDistance',0.5); %encontra os picos S
[Qy,Qx]  = findpeaks(- y2(pos_Q), t(pos_Q),'MinPeakDistance',0.5); %encontra os picos Q

Q_S=sum(Sx([1:length(Qx)])-Qx)/length(Qx); %calcula tempo do intervalo QS

string = ['intervalo QS= ',num2str(Q_S*1000),' ms'];disp(string);

%% Encontra ondas P e T
t_matriz_qx=repmat(t,length(Qx),1);
t_matriz_sx=repmat(t,length(Sx),1);
Qx_matriz= repmat(Qx',1,length(t));
Sx_matriz= repmat(Sx',1,length(t));

media_Intervalo_QQ=median(intervalo);

intervalo_P=sum(t_matriz_qx<Qx_matriz & t_matriz_qx>(Qx_matriz - (media_Intervalo_QQ/2)));
intervalo_T=sum(t_matriz_sx>Sx_matriz & t_matriz_sx<(Sx_matriz + (media_Intervalo_QQ/2)));

pos_P=find(intervalo_P);
pos_T=find(intervalo_T);

[Py,Px]  = findpeaks( (y2(pos_P)), t(pos_P),'MinPeakDistance',0.5); %encontra os picos P
[Ty,Tx]  = findpeaks( abs(y2(pos_T)), t(pos_T),'MinPeakDistance',0.5); %encontra os picos T





%% Transportar os pontos do sinal y2 para o sinal original

t_matrix_R=repmat(t,length(peak_x),1);
R_matrix= repmat(peak_x',1,length(t));
pos_R_original=find(sum((R_matrix>t_matrix_R - (time_signal_s/length(y)/2) & R_matrix<= t_matrix_R + (time_signal_s/length(y)/2))));

t_matrix_Q=repmat(t,length(Qx),1);
Q_matrix= repmat(Qx',1,length(t));
pos_Q_original=find(sum((Q_matrix>t_matrix_Q - (time_signal_s/length(y)/2) & Q_matrix<= t_matrix_Q + (time_signal_s/length(y)/2))));

t_matrix_S=repmat(t,length(Sx),1);
S_matrix= repmat(Sx',1,length(t));
pos_S_original=find(sum((S_matrix>t_matrix_S - (time_signal_s/length(y)/2) & S_matrix<= t_matrix_S + (time_signal_s/length(y)/2))));

t_matrix_P=repmat(t,length(Px),1);
P_matrix= repmat(Px',1,length(t));
pos_P_original=find(sum((P_matrix>t_matrix_P - (time_signal_s/length(y)/2) & P_matrix<= t_matrix_P + (time_signal_s/length(y)/2))));

t_matrix_T=repmat(t,length(Tx),1);
T_matrix= repmat(Tx',1,length(t));
pos_T_original=find(sum((T_matrix>t_matrix_T - (time_signal_s/length(y)/2) & T_matrix<= t_matrix_T + (time_signal_s/length(y)/2))));

%% plot
figure()
% PLOT
subplot(2,1,1)
plot(t,y2,'b');
hold on;
grid on;
xlim([0 10])
title('Sinal apos filtro PA e PB')
plot(Px, Py,'xk')% marca onda P
plot(Qx,-Qy,'vb')% marca onda Q
plot(peak_x,peak_y,'x');% marca onda R
plot(Sx,-Sy,'vg')% marca onda S
plot(Tx, Ty,'xr')% marca onda T

subplot(2,1,2)
plot(t,y,'b');
grid on;
hold on;
title('Marcação dos pontos no sinal apos PA')
plot(t(pos_P_original),y(pos_P_original), 'xk');% marca onda P 
plot(t(pos_Q_original),y(pos_Q_original), 'vb');% marca onda Q
plot(t(pos_R_original),y(pos_R_original), 'x');% marca onda R
plot(t(pos_S_original),y(pos_S_original), 'vg');% marca onda S
plot(t(pos_T_original),y(pos_T_original), 'xr');% marca onda T
xlim([0 10])
