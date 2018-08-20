clc; clear all; close all;
Name = '100m'; %arritmia, onda T negativa
%Name = '101m'; %arritmia
% Name = '102m'; %arritmia
Name = '103m'; %arritmia
% Name = '104m'; %arritmia
% Name = '105m'; %arritmia
% Name = '106m'; %arritmia
% Name = '107m'; %arritmia
% Name = '108m'; %arritmia
% Name = '109m'; %arritmia
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
%%
load('filtros.mat');
if (freqint(1,1)<200)  %realiza interpolação caso a frequencia de amostragem do sinal original seja baixa demais para a decomposição  
    r=3;
    sinal_original= interp(val,r);    
else
    r=1;
    sinal_original= val;
end

time_signal_s=length(sinal_original)/(freqint(1,1)*r); %calcula o tempo do sinal em segundos
t=0:time_signal_s/length(sinal_original):time_signal_s-(time_signal_s/length(sinal_original));% cria vetor de tempo para utilizar nos plots

y=filter(FIR_HP_300,[sinal_original zeros(1,150)]);%filtragem das frequencias menores que 0.5Hz FIR ordem 300
y=y([151:end]);

%% Decomposição, reconstrução, filtragem e downsampling
wavelet='sym1';

[a1,d1] = dwt(y,wavelet);
[a2,d2] = dwt(a1,wavelet);
[a3,d3] = dwt(a2,wavelet);


reconst_a2 = idwt(a3,0*d3,wavelet,length(a2));
reconst_a1 = idwt(reconst_a2,0*d2,wavelet,length(a1));
reconst = idwt(reconst_a1,0*d1,wavelet,length(y));

y2 = filtfilt(SOS_PB,G_PB,reconst); %FILTRO PASSA BAIXAS
y2_downsample = downsample(y2,4); % REALIZA O DOWNSAMPLING
t2_downsample=0:time_signal_s/length(y2_downsample):time_signal_s-(time_signal_s/length(y2_downsample));%cria vetor de tempo do downsampling

y3 = filtfilt(SOS_PB,G_PB,y); %FILTRAR COM O PB SEM PASSAR PELA DECOMPOSICAO
y3_downsample = downsample(y3,4);
t3_downsample=0:time_signal_s/length(y3_downsample):time_signal_s-(time_signal_s/length(y3_downsample));%cria vetor de tempo do downsampling




%% DENOISE------------------------------
%%%apply Wavelet Transform
[C,L]=wavedec(y2,3,'sym5');
%[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
% %%%Denoise 
[thr,sorh,keepapp]=ddencmp('den','wv',y2);
cleanecg=wdencmp('gbl',C,L,'sym5',3,thr,sorh,keepapp);
%%
% figure()
% plot(t,y2)
% hold on;
% plot (t3_downsample,y3_downsample)

%% plots
figure()
plot (t,y,'k')
hold on
plot (t,reconst,'b')
plotbrowser('on');
grid on;
plot(t,y2)
plot(t2_downsample,y2_downsample);
plot(t,y3);
plot(t,cleanecg)
legend('original','a3','a3 filtrado altas frequencias','downsampling a3 filtrado','sem wavelet','Denoised')
xlim([0 10])
ylim([min(y) max(y)])



%% TRANSFORMADA DE FOURIER
L=length(sinal_original);
Fs=L/time_signal_s;
T=1/Fs;
Y=fft(sinal_original);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure()
subplot(2,1,1)
plot(f,P1)
grid on
title('Espectro de frequencia do sinal original')
xlim([0 100])
ylim([0 0.08])
%---------------------------------------------------------

L=length(y3_downsample);
Fs=L/time_signal_s;
T=1/Fs;
Y=fft(y3_downsample);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
subplot(2,1,2)
plot(f,P1)
grid on
title('Espectro de frequencia do sinal pre-processado PA+PB')
xlim([0 100])

%% Encontra picos R sinal y3_downsample

m=max(y3_downsample);
[peak_y, peak_x] = findpeaks(y3_downsample,t3_downsample,'minpeakheight',m*0.3,'MinPeakDistance',0.3);

% analise dos intervalos RR do a3_atenuado-------------------------------------
peak_aux2 = peak_x(2:end);
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

t_matriz=repmat(t3_downsample,length(peak_x),1); %cria uma matriz de tamanho: length(peak_x) x length(t3) apenas repetindo o t3
peak_x_matrix= repmat(peak_x',1,length(t3_downsample)); %cria uma matriz de tamanho: length(peak_x3') x length(t) apenas repetindo o peak_x3'

S=sum(t_matriz>peak_x_matrix & t_matriz<(peak_x_matrix+0.15)); %localiza todas as amostras que estão no intervalo de 0,15s após cada onda R
Q=sum(t_matriz<peak_x_matrix & t_matriz>(peak_x_matrix-0.15));%localiza todas as amostras que estão no intervalo de 0,15s antes de cada onda R


pos_S=find(S); %pega a posição das amostras
pos_Q=find(Q); %pega a posição das amostras


[Sy,Sx]  = findpeaks(- y3_downsample(pos_S), t3_downsample(pos_S),'MinPeakDistance',0.5); %encontra os picos S
[Qy,Qx]  = findpeaks(- y3_downsample(pos_Q), t3_downsample(pos_Q),'MinPeakDistance',0.5); %encontra os picos Q


%Q_S=sum(Sx([1:length(Qx)])-Qx)/length(Qx); %calcula tempo do intervalo QS
Q_S=sum(Sx - Qx([1:length(Sx)]))/length(Qx); %calcula tempo do intervalo QS
string = ['intervalo QS= ',num2str(Q_S*1000),' ms'];disp(string);

%% Encontra ondas P e T
t_matriz_qx=repmat(t3_downsample,length(Qx),1);
t_matriz_sx=repmat(t3_downsample,length(Sx),1);
Qx_matriz= repmat(Qx',1,length(t3_downsample));
Sx_matriz= repmat(Sx',1,length(t3_downsample));

media_Intervalo_QQ=median(intervalo);

intervalo_P=sum(t_matriz_qx<Qx_matriz & t_matriz_qx>(Qx_matriz - (media_Intervalo_QQ/3)));
intervalo_T=sum(t_matriz_sx>(Sx_matriz + 0.05)& t_matriz_sx<(Sx_matriz + (media_Intervalo_QQ/3)));

pos_P=find(intervalo_P);
pos_T=find(intervalo_T);

[Py,Px]  = findpeaks( (y3_downsample(pos_P)), t3_downsample(pos_P),'MinPeakDistance',0.5); %encontra os picos P
[Ty,Tx]  = findpeaks( abs(y3_downsample(pos_T)), t3_downsample(pos_T),'MinPeakDistance',0.5); %encontra os picos T

%% LOCALIZAR A POSICAO DOS PONTOS

t3_matrix_P=repmat(t3_downsample,length(Px),1);
P_matrix= repmat(Px',1,length(t3_downsample));
posicao_P=find(sum(t3_matrix_P==P_matrix)) ;
clear t3_matrix_P; clear P_matrix;

t3_matrix_Q=repmat(t3_downsample,length(Qx),1);
Q_matrix= repmat(Qx',1,length(t3_downsample));
posicao_Q=find(sum(t3_matrix_Q==Q_matrix)) ;
clear t3_matrix_Q; clear Q_matrix;

t3_matrix_R=repmat(t3_downsample,length(peak_x),1);
R_matrix= repmat(peak_x',1,length(t3_downsample));
posicao_R=find(sum(t3_matrix_R==R_matrix)) ;
clear t3_matrix_R; clear R_matrix;

t3_matrix_S=repmat(t3_downsample,length(Sx),1);
S_matrix= repmat(Sx',1,length(t3_downsample));
posicao_S=find(sum(t3_matrix_S==S_matrix)) ;
clear t3_matrix_S; clear S_matrix;

t3_matrix_T=repmat(t3_downsample,length(Tx),1);
T_matrix= repmat(Tx',1,length(t3_downsample));
posicao_T=find(sum(t3_matrix_T==T_matrix)) ;
clear t3_matrix_T; clear T_matrix;

%% plot
figure()
% PLOT
subplot(2,1,1)
plot(t3_downsample,y3_downsample,'b');
hold on;
grid on;
%xlim([0 10])
title('Sinal y3downsample')

plot(t3_downsample(posicao_P), y3_downsample(posicao_P),'xk')% marca onda P
plot(t3_downsample(posicao_Q), y3_downsample(posicao_Q),'vb')% marca onda Q
plot(t3_downsample(posicao_R), y3_downsample(posicao_R),'x')% marca onda R
plot(t3_downsample(posicao_S), y3_downsample(posicao_S),'vg')% marca onda S
plot(t3_downsample(posicao_T), y3_downsample(posicao_T),'xr')% marca onda T

subplot(2,1,2)
plot(t,y,'b');
grid on;
hold on;
title('Marcacoes no sinal original')
plot(t(posicao_P*4 -2),y(posicao_P*4 -2), 'xk');% marca onda P 
plot(t(posicao_Q*4 -2),y(posicao_Q*4 -2), 'vb');% marca onda Q
plot(t(posicao_R*4 -2),y(posicao_R*4 -2), 'x');% marca onda R
plot(t(posicao_S*4 -2),y(posicao_S*4 -2), 'vg');% marca onda S
plot(t(posicao_T*4 -2),y(posicao_T*4 -2), 'xr');% marca onda T
%xlim([0 10])
