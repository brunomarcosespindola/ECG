clc; clear all; close all;
%Name ='118e00m'; %NOISE STRESS


%% 
%Name= {'100m','101m','102m','103m','104m','105m','107m','108m','109m','111m'};% arritmia
%Name= {'112m','114m','115m','116m','117m','118m','121m','122m','123m','124m'};% arritmia
%Name= {'200m','201m','202m','203m','205m','209m','210m','212m','214m','219m'};% arritmia
%Name= {'113m','213m','215m','217m'};% arritmia

%% TODOS SINAIS DE ARRITMIA
%Name= {'100m','101m','102m','103m','104m','105m','107m','108m','109m','111m','112m','114m','115m','116m','117m','118m','121m','122m','123m','124m','200m','201m','202m','203m','205m','209m','210m','212m','214m','219m','113m','213m','215m','217m'};
Name= {'100m'}

%% ---- sinusal
%Name= {'16265m','16272m','16420m','16483m','16539m','16273m','16773m','16786m','16795m','17052m','17453m','18177m','18184m'};
%Name={'19090m','19093m','19140m','19830m'}; %sinusal com falha
%% 

%SINAIS MAIORES QUE 1 MIN
%Name={'116m_30min'}
%Name={'103m_30m','100m30m','112m_30m'};


resultados=struct();
%%
for cont=1:length(Name)

    
infoName = strcat(char(Name(cont)), '.info');
matName = strcat(char(Name(cont)), '.mat');

Octave = exist('OCTAVE_VERSION');
load(matName);

resultados(cont).Name=matName;
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

%% filtragem das frequencias menores que 0.5Hz FIR ordem 300
y=filter(FIR_HP_300,[sinal_original zeros(1,150)]);%
y=y([151:end]);

y3 = filtfilt(SOS_PB,G_PB,y); %FILTRAR COM O PB SEM PASSAR PELA DECOMPOSICAO
 y3_downsample = downsample(y3,4);%% downsample
 t3_downsample=0:time_signal_s/length(y3_downsample):time_signal_s-(time_signal_s/length(y3_downsample));%cria vetor de tempo do downsampling

% figure()
% plot(t,y,'b')
% xlabel('segundos')
% ylabel('mV')
% grid on
% hold on
% plot(t,y3,'r')
% legend('original','filtrado')
% xlim([6 8])
%% TRANSFORMADA DE FOURIER

L=length(y3_downsample);
Fs=L/time_signal_s;

T=1/Fs;
Y=fft(y3_downsample(1:60*Fs));
L=length(y3_downsample(1:60*Fs));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

% figure()
% subplot(1,2,2)
% plot(f,(abs(P1)))
% xlabel('Hz')
% ylabel('mV')
% grid on
% title(['Domínio da frequência'])
% xlim([0 50])
% subplot(1,2,1)
% plot(t,y3)
% title(['Domínio do tempo'])
% xlim([0 2])
% xlabel('segundos')
% ylabel('mV')
% grid on



%--------Comparação em frequencia com sinal Sinusal 

load('sinusal_FFT.mat')
 %figure()
 %plot(sinusal_FFT.f(1:length(f)),sinusal_FFT.P(1:length(f)));
 corr=corrcoef(sinusal_FFT.P(1:length(f)),P1);
 resultados(cont).correlacao=corr(1,2);

%% JANELAMENTO

janela= Fs*60;
ii=ceil(length(y3)/janela);

for jj=1:ii
    if(jj==ii && ii ~=1)
       y_janela= y3([((jj-1)*janela)+(jj-1): end]);
       t_janela= t([((jj-1)*janela)+(jj-1): end]);       
    else   
        y_janela= y3([((jj-1)*janela)+(jj):((((jj-1)*janela)+(jj-1))+janela)]);
        t_janela=t([((jj-1)*janela)+(jj):((((jj-1)*janela)+(jj-1))+janela)]);
    end

%% Encontra picos R sinal y3
m=max(y_janela);
[peak_y, peak_x] = findpeaks(y_janela,t_janela,'minpeakheight',m*0.3,'MinPeakDistance',0.3);

%% Encontra onda Q e S do sinal PREPROCESSADO

t_matriz=repmat(t_janela,length(peak_x),1); %cria uma matriz de tamanho: length(peak_x) x length(t3) apenas repetindo o t3
peak_x_matrix= repmat(peak_x',1,length(t_janela)); %cria uma matriz de tamanho: length(peak_x3') x length(t) apenas repetindo o peak_x3'

S=sum(t_matriz>peak_x_matrix & t_matriz<(peak_x_matrix+0.15)); %localiza todas as amostras que estão no intervalo de 0,15s após cada onda R
Q=sum(t_matriz<peak_x_matrix & t_matriz>(peak_x_matrix-0.15));%localiza todas as amostras que estão no intervalo de 0,15s antes de cada onda R


pos_S=find(S); %pega a posição das amostras
pos_Q=find(Q); %pega a posição das amostras


[Sy,Sx]  = findpeaks(- y_janela(pos_S), t_janela(pos_S),'MinPeakDistance',0.5); %encontra os picos S
[Qy,Qx]  = findpeaks(- y_janela(pos_Q), t_janela(pos_Q),'MinPeakDistance',0.5); %encontra os picos Q

%% Encontra ondas P e T
t_matriz_qx=repmat(t_janela,length(Qx),1);
t_matriz_sx=repmat(t_janela,length(Sx),1);
Qx_matriz= repmat(Qx',1,length(t_janela));
Sx_matriz= repmat(Sx',1,length(t_janela));


peak_aux2 = peak_x(2:end);
intervalo = peak_aux2 - peak_x(1:end-1);
media_Intervalo_RR=median(intervalo);

intervalo_P=sum(t_matriz_qx<Qx_matriz & t_matriz_qx>(Qx_matriz - (media_Intervalo_RR/3)));
intervalo_T=sum(t_matriz_sx>(Sx_matriz + 0.10)& t_matriz_sx<(Sx_matriz + (media_Intervalo_RR/2)));

pos_P=find(intervalo_P);
pos_T=find(intervalo_T);


[Py,Px]  = findpeaks( (y_janela(pos_P)), t_janela(pos_P),'MinPeakDistance',0.5); %encontra os picos P
[Ty,Tx]  = findpeaks( abs(y_janela(pos_T)), t_janela(pos_T),'MinPeakDistance',0.5); %encontra os picos T


%% LOCALIZA A POSICAO DOS PONTOS

t3_matrix_P=repmat(t,length(Px),1);
P_matrix= repmat(Px',1,length(t));
posicao_P_parcial=find(sum(t3_matrix_P==P_matrix)) ;
clear t3_matrix_P; clear P_matrix;

t3_matrix_Q=repmat(t,length(Qx),1);
Q_matrix= repmat(Qx',1,length(t));
posicao_Q_parcial=find(sum(t3_matrix_Q==Q_matrix)) ;
clear t3_matrix_Q; clear Q_matrix;

t3_matrix_R=repmat(t,length(peak_x),1);
R_matrix= repmat(peak_x',1,length(t));
posicao_R_parcial=find(sum(t3_matrix_R==R_matrix)) ;
clear t3_matrix_R; clear R_matrix;

t3_matrix_S=repmat(t,length(Sx),1);
S_matrix= repmat(Sx',1,length(t));
posicao_S_parcial=find(sum(t3_matrix_S==S_matrix)) ;
clear t3_matrix_S; clear S_matrix;

t3_matrix_T=repmat(t,length(Tx),1);
T_matrix= repmat(Tx',1,length(t));
posicao_T_parcial=find(sum(t3_matrix_T==T_matrix)) ;
clear t3_matrix_T; clear T_matrix;

%% Agrupar os resultados do janelamento
if jj==1
    posicao_P=posicao_P_parcial;
    posicao_Q=posicao_Q_parcial;
    posicao_R=posicao_R_parcial;
    posicao_S=posicao_S_parcial;
    posicao_T=posicao_T_parcial;
    
    Ry_total=peak_y;
    Rx_total=peak_x;  
    Sy_total= Sy;
    Sx_total= Sx;
    Qy_total=Qy;
    Qx_total=Qx;
    Py_total=Py;
    Px_total=Px;
    Ty_total=Ty;
    Tx_total=Tx;
else
    posicao_P=[posicao_P posicao_P_parcial];
    posicao_Q=[posicao_Q posicao_Q_parcial];
    posicao_R=[posicao_R posicao_R_parcial];
    posicao_S=[posicao_S posicao_S_parcial];
    posicao_T=[posicao_T posicao_T_parcial];
    Ry_total=[Ry_total peak_y];
    Rx_total=[Rx_total peak_x];    
    Sy_total=[Sy_total Sy];
    Sx_total=[Sx_total Sx];    
    Qy_total=[Qy_total Qy];
    Qx_total=[Qx_total Qx];
    Py_total=[Py_total Py];
    Px_total=[Px_total Px];   
    Ty_total=[Ty_total Ty];
    Tx_total=[Tx_total Tx];
end

end
%% 

% analise dos intervalos RR

Rx_aux = Rx_total(2:end);
intervalo_total = Rx_aux - Rx_total(1:end-1);
resultados(cont).media_RR_s=mean(intervalo_total);
bpm= (60./intervalo_total);
resultados(cont).Media_bpm= sum(bpm)/length(bpm);
resultados(cont).variancia_RR=var(intervalo_total);

% figure()
% plot(peak_x(2:end),intervalo);
% title('intervalos de tempo RR')
% xlabel('Segundos')
% ylabel('RR (seg)')
% grid on;

%-------------------- CRIA ALERTA PARA ONDA R FORA DO RITMO
intervalo_aux= intervalo_total(2:end);
dif_intervalo=intervalo_aux - intervalo_total(1:end-1);
t_dif=Rx_total(3:end);

pontos_de_alerta=find(dif_intervalo>0.1);
if length(pontos_de_alerta)>0
    if pontos_de_alerta(1)==1;
        pontos_de_alerta=pontos_de_alerta(2:end);
    end
end
alerta=t_dif(pontos_de_alerta -1);
resultados(cont).n_alertas=length(alerta);

%% ----------------INTERVALO QS
if (length(Sx_total)> length(Qx_total));
    vetorQS=Sx_total([2:length(Qx_total)+1])-Qx_total;
    % Q_S=sum(vetorQS)/length(Qx_total); %calcula tempo do intervalo QS
else (length(Sx_total)< length(Qx_total));   
       vetorQS=Sx_total - Qx_total([1:length(Sx_total)]);
    % Q_S=sum(vetorQS)/length(Qx_total);     
end

resultados(cont).intervalo_QS_ms=median(vetorQS)*1000;
resultados(cont).variancia_QS_ms=var(vetorQS);
%string = ['intervalo QS= ',num2str(Q_S*1000),' ms'];disp(string);
%resultados(cont).intervalo_QS_ms=Q_S*1000;


%% plot
figure()
% PLOT

plot(t,y3,'b');
xlabel('segundos')
ylabel('mV')
hold on;
grid on;

title(['Sinal y3downsample: ',resultados(cont).Name])

plot(t(posicao_P), y3(posicao_P),'xk')% marca onda P
plot(t(posicao_Q), y3(posicao_Q),'vm')% marca onda Q
plot(t(posicao_R), y3(posicao_R),'x')% marca onda R
plot(t(posicao_S), y3(posicao_S),'vg')% marca onda S
plot(t(posicao_T), y3(posicao_T),'xr')% marca onda T
legend('ECG','onda P','onda Q','onda R','onda S','onda T')
%xlim([2 8])
y_alerta=repmat(m/2,1,length(alerta));
text(alerta,y_alerta,'Alerta')


% subplot(2,1,2)
% plot(t,y,'b');
% grid on;
% hold on;
% title('Marcacoes no sinal original')
% plot(t(posicao_P*4 -2),y(posicao_P*4 -2), 'xk');% marca onda P
% plot(t(posicao_Q*4 -2),y(posicao_Q*4 -2), 'vm');% marca onda Q
% plot(t(posicao_R*4 -2),y(posicao_R*4 -2), 'x');% marca onda R
% plot(t(posicao_S*4 -2),y(posicao_S*4 -2), 'vg');% marca onda S
% plot(t(posicao_T*4 -2),y(posicao_T*4 -2), 'xr');% marca onda T
%xlim([0 10]

%% Diagnostico

if resultados(cont).correlacao>0.2
    resultados(cont).sinusal='sinal pode ser sinusal';   
else
    resultados(cont).sinusal='sinal não sinusal';  
end

if resultados(cont).intervalo_QS_ms>100
    resultados(cont).bloqueio_ramo='Possível bloqueio de ramo';
else
     resultados(cont).bloqueio_ramo='não';
end

if resultados(cont).Media_bpm>100
    resultados(cont).Ritmo_cardiaco='Arritmia:Taquicardia';
elseif resultados(cont).Media_bpm<50
    resultados(cont).Ritmo_cardiaco='Arritmia:Bradicardia';
else
    
    if resultados(cont).n_alertas>0
        resultados(cont).Ritmo_cardiaco=['Arritmia: Possui ',num2str(resultados(cont).n_alertas),' alertas'];
    else
        resultados(cont).Ritmo_cardiaco='Batimentos Regulares';
    end
end
    




%% COMPARA RESULTADO PHYSIONET
%faz a leitura do arquivo e converte o tempo localizado em segundos
tabela=importdata('100.txt');
R_physionet=tabela.textdata(:,1);
R_physionet=R_physionet(2:end)';
anotacao=tabela.textdata(:,3);
anotacao=anotacao(2:end)';
clear tabela;
[minutos,segundos]=strtok(R_physionet,':');
segundos = regexprep(segundos,':','');
minutos=str2double(minutos);
segundos=str2double(segundos);
tempo=(minutos*60)+segundos;
%%

tempo=tempo(2:end);

diferenca=Rx_total-tempo;
media_dif=mean(diferenca)

end
%%
for ii=1:length(resultados)
resultados(ii)
end

