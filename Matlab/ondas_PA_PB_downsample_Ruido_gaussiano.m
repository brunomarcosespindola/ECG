clc; clear all; close all;
Name = '100m'; %arritmia, onda T negativa
Name = '101m'; %arritmia
%Name = '102m'; %arritmia
%Name = '103m'; %arritmia
%Name = '104m'; %arritmia
%Name = '105m'; %arritmia
%Name = '106m'; %arritmia TRAVA O ALGORITMO
%Name = '107m'; %arritmia
%Name = '108m'; %arritmia
%Name = '109m'; %arritmia
%Name = '111m'; %arritmia
%Name = '16265m';% sinusal
Name ='16272m'; %sinusal -possui um ruído forte no meio
Name ='16420m'; %sinusal
Name ='16483m'; %sinusal
Name ='16539m'; %sinusal
%Name ='118e00m'; %NOISE STRESS

%snr = 100;
cont=0;
for snr=51:-5:1
    cont=cont + 1;
    
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
    
    
    
    %% ADICIONA RUÍDO GAUSSIANO
    
    
    r = randn(1, length(sinal_original));
    % Adicao da senoide ao ruído ajustado para atender SNR
    % energia do ruido
    er = sum(r.^2);
    % energia do sinal
    es = sum(sinal_original.^2);
    % SNR desejada, dB
    
    k = 10^(snr/10); % deseja-se que es/er seja k, logo:
    r = r * (es / (er*k))^0.5;
    % energia ajustada do ruido
    er = sum(r.^2);
    % SNR
    snr_atual(cont)=10*log10(es/er);
    % sinal resultante
    y_ruido = sinal_original + r;
    
    % figure()
    % plot(t,sinal_original)
    % hold on
    % plot(t,y_ruido)
    % legend('original','ruido')
    
    %% filtragem das frequencias menores que 0.5Hz FIR ordem 300
    y=filter(FIR_HP_300,[y_ruido zeros(1,150)]);%
    y=y([151:end]);
    
    %% Decomposição, reconstrução, filtragem e downsampling
    % wavelet='sym1';
    %
    % [a1,d1] = dwt(y,wavelet);
    % [a2,d2] = dwt(a1,wavelet);
    % [a3,d3] = dwt(a2,wavelet);
    %
    %
    % reconst_a2 = idwt(a3,0*d3,wavelet,length(a2));
    % reconst_a1 = idwt(reconst_a2,0*d2,wavelet,length(a1));
    % reconst = idwt(reconst_a1,0*d1,wavelet,length(y));
    %
    % y2 = filtfilt(SOS_PB,G_PB,reconst); %FILTRO PASSA BAIXAS
    % y2_downsample = downsample(y2,4); % REALIZA O DOWNSAMPLING
    % t2_downsample=0:time_signal_s/length(y2_downsample):time_signal_s-(time_signal_s/length(y2_downsample));%cria vetor de tempo do downsampling
    %
    y3 = filtfilt(SOS_PB,G_PB,y); %FILTRAR COM O PB SEM PASSAR PELA DECOMPOSICAO
    y3_downsample = downsample(y3,4);
    t3_downsample=0:time_signal_s/length(y3_downsample):time_signal_s-(time_signal_s/length(y3_downsample));%cria vetor de tempo do downsampling
    
    % %% FILTRO PB MUITO SELETIVO
    %
    % load('filtro_pf.mat');
    %
    % y_rejeita_faixa = filtfilt(SOS_2,G_2,y3_downsample);
    %
    % figure()
    % plot (t3_downsample,y3_downsample)
    % hold on
    % plot (t3_downsample,y_rejeita_faixa)
    %
    %
    %
    % L=length(y_rejeita_faixa);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(y_rejeita_faixa);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    % figure()
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do FILTRO')
    % % xlim([0 100])
    % % ylim([0 0.08])
    
    %% DENOISE------------------------------
    % %%%apply Wavelet Transform
    % [C,L]=wavedec(y2,1,'sym5');
    % %[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
    % % %%%Denoise
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % cleanecg=wdencmp('gbl',C,L,'sym5',1,thr,sorh,keepapp);
    
    %% Comparação entre denoiseds
    % wavelet='db1';
    % ord=3;
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % clean1=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='db2';
    %
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % clean2=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='db3';
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % clean3=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='db4';
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % clean4=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='db5';
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % clean5=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='db6';
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % clean6=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % figure()
    % plot (t,y,'k')
    % hold on
    % plot (t,clean1)
    % plot (t,clean2)
    % plot (t,clean3)
    % plot (t,clean4)
    % plot (t,clean5)
    % plot (t,clean6)
    % plotbrowser('on');
    % grid on;
    % xlim([0 10])
    
    %%
    % wavelet='db4';
    % ord=3;
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % compara1=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='sym4';
    % ord=3;
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % compara2=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % wavelet='haar';
    % ord=3;
    % [C,L]=wavedec(y2,ord,wavelet);
    % [thr,sorh,keepapp]=ddencmp('den','wv',y2);
    % compara3=wdencmp('gbl',C,L,wavelet,ord,thr,sorh,keepapp);
    %
    % figure()
    % plot (t,y,'k')
    % hold on
    % plot (t,compara1)
    % plot (t,compara2)
    % plot (t,compara3)
    % plotbrowser('on');
    % grid on;
    % xlim([0 10])
    %
    
    
    
    %% comparação entre original e subamostrado
    % figure()
    % plot(t,y2)
    % hold on;
    % plot (t3_downsample,y3_downsample)
    
    %% plots
    % figure()
    % plot (t,y,'k')
    % hold on
    % plot (t,reconst,'b')
    % plotbrowser('on');
    % grid on;
    % plot(t,y2)
    % plot(t2_downsample,y2_downsample);
    % plot(t,y3);
    % plot(t,cleanecg)
    % legend('original','a3','a3 filtrado altas frequencias','downsampling a3 filtrado','sem wavelet','Denoised')
    % xlim([0 10])
    % ylim([min(y) max(y)])
    
    %% ----------Plots das decomposicoes------------------------------
    %
    % t1=0:time_signal_s/length(a1):time_signal_s-(time_signal_s/length(a1));
    % t2=0:time_signal_s/length(a2):time_signal_s-(time_signal_s/length(a2));
    % t3=0:time_signal_s/length(a3):time_signal_s-(time_signal_s/length(a3));
    %
    % janela=30; %define o intervalo dos prints
    % figure()
    % subplot(3,2,1);
    % plot(t1,a1);
    % title('a1')
    % xlim([0 janela])
    % subplot(3,2,2);
    % plot(t1,d1);
    % title('d1')
    % xlim([0 janela])
    % xlabel('Seconds')
    % ylabel('Amplitude')
    %
    % subplot(3,2,3);
    % plot(t2,a2);
    % title('a2')
    % xlim([0 janela])
    % subplot(3,2,4);
    % plot(t2,d2);
    % title('d2')
    % xlim([0 janela])
    %
    % subplot(3,2,5);
    % plot(t3,a3);
    % title('a3')
    % xlim([0 janela])
    % subplot(3,2,6);
    % plot(t3,d3);
    % title('d3')
    % xlim([0 janela])
    
    %% TRANSFORMADA DE FOURIER
    % L=length(sinal_original);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(sinal_original);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    % figure()
    % subplot(3,1,1)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do sinal original')
    % xlim([0 100])
    % ylim([0 0.08])
    % %---------------------------------------------------------
    %
    % L=length(y3_downsample);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(y3_downsample);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    % subplot(3,1,2)
    % plot(f,P1)
    % grid on
    % title('Espectro de frequencia do sinal pre-processado PA+PB')
    % xlim([0 100])
    % ylim([0 0.08])
    %
    % %----------------------------
    %
    % L=length(cleanecg);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(cleanecg);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    % subplot(3,1,3)
    % plot(f,P1)
    % grid on
    % title('Espectro de frequencia do sinal denoised wavelet')
    % xlim([0 100])
    % ylim([0 0.08])
    
    
    %% FFT das decomposições
    
    % L=length(a1);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(a1);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    % figure()
    % subplot(3,2,1)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do a1')
    % xlim([0 100])
    % ylim([0 0.08])
    %
    % %---------------------
    %
    % L=length(d1);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(d1);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    %
    % subplot(3,2,2)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do d1')
    % xlim([0 100])
    % ylim([0 0.08])
    %
    % %---------------------
    %
    % L=length(a2);
    % Fs=L/time_signal_s;
    % T=1/Fs;
    % Y=fft(a2);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    %
    % subplot(3,2,3)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do a2')
    % xlim([0 100])
    % ylim([0 0.08])
    %
    % %---------------------
    %
    %  L=length(d2);
    %  Fs=L/time_signal_s;
    %  T=1/Fs;
    % Y=fft(d2);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    %
    % subplot(3,2,4)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do d2')
    % xlim([0 100])
    % ylim([0 0.08])
    %
    % %---------------------
    %
    %  L=length(a3);
    %  Fs=L/time_signal_s;
    %  T=1/Fs;
    % Y=fft(a3);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    %
    % subplot(3,2,5)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do a3')
    % xlim([0 100])
    % ylim([0 0.08])
    %
    % %---------------------
    %
    %  L=length(d3);
    %  Fs=L/time_signal_s;
    %  T=1/Fs;
    % Y=fft(d3);
    %
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    %
    % f = Fs*(0:(L/2))/L;
    %
    % subplot(3,2,6)
    % plot(f,P1)
    % grid on
    %
    % title('Espectro de frequencia do d3')
    % xlim([0 100])
    % ylim([0 0.08])
    
    
    
    
    %% Encontra picos R sinal y3_downsample
    
    m=max(y3_downsample);
    [peak_y, peak_x] = findpeaks(y3_downsample,t3_downsample,'minpeakheight',m*0.3,'MinPeakDistance',0.3);
    
    % analise dos intervalos RR
    peak_aux2 = peak_x(2:end);
    intervalo = peak_aux2 - peak_x(1:end-1);
    media_RR_s(cont)=mean(intervalo);
    bpm= (60./intervalo);
    Media_bpm(cont)= sum(bpm)/length(bpm);
    variancia_RR(cont)=var(intervalo);
    % figure()
    % plot(peak_x(2:end),intervalo);
    % title('intervalos de tempo RR')
    % xlabel('Segundos')
    % ylabel('RR (seg)')
    % grid on;
    
    
    
    
    %-------------------- CRIA ALERTA PARA ONDA R FORA DO RITMO
    intervalo_aux= intervalo(2:end);
    dif_intervalo=intervalo_aux - intervalo(1:end-1);
    t_dif=peak_x(3:end);
    
    pontos_de_alerta=find(dif_intervalo>0.1);
    
    alerta=t_dif(pontos_de_alerta-1);
    
    n_alertas(cont)=length(alerta);
    
    
    
    
    %% Encontra onda Q e S do sinal PREPROCESSADO
    
    t_matriz=repmat(t3_downsample,length(peak_x),1); %cria uma matriz de tamanho: length(peak_x) x length(t3) apenas repetindo o t3
    peak_x_matrix= repmat(peak_x',1,length(t3_downsample)); %cria uma matriz de tamanho: length(peak_x3') x length(t) apenas repetindo o peak_x3'
    
    S=sum(t_matriz>peak_x_matrix & t_matriz<(peak_x_matrix+0.15)); %localiza todas as amostras que estão no intervalo de 0,15s após cada onda R
    Q=sum(t_matriz<peak_x_matrix & t_matriz>(peak_x_matrix-0.15));%localiza todas as amostras que estão no intervalo de 0,15s antes de cada onda R
    
    
    pos_S=find(S); %pega a posição das amostras
    pos_Q=find(Q); %pega a posição das amostras
    
    
    [Sy,Sx]  = findpeaks(- y3_downsample(pos_S), t3_downsample(pos_S),'MinPeakDistance',0.5); %encontra os picos S
    [Qy,Qx]  = findpeaks(- y3_downsample(pos_Q), t3_downsample(pos_Q),'MinPeakDistance',0.5); %encontra os picos Q
    
    
%    Q_S=sum(Sx([1:length(Qx)])-Qx)/length(Qx); %calcula tempo do intervalo QS
    Q_S=sum(Sx - Qx([1:length(Sx)]))/length(Qx); %calcula tempo do intervalo QS
    %string = ['intervalo QS= ',num2str(Q_S*1000),' ms'];disp(string);
    intervalo_QS_ms(cont)=Q_S*1000;
    %% Encontra ondas P e T
    t_matriz_qx=repmat(t3_downsample,length(Qx),1);
    t_matriz_sx=repmat(t3_downsample,length(Sx),1);
    Qx_matriz= repmat(Qx',1,length(t3_downsample));
    Sx_matriz= repmat(Sx',1,length(t3_downsample));
    
    media_Intervalo_QQ=median(intervalo);
    
    intervalo_P=sum(t_matriz_qx<Qx_matriz & t_matriz_qx>(Qx_matriz - (media_Intervalo_QQ/3)));
    intervalo_T=sum(t_matriz_sx>(Sx_matriz + 0.10)& t_matriz_sx<(Sx_matriz + (media_Intervalo_QQ/2)));
    
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
    y_alerta=repmat(m/2,1,length(alerta));
    text(alerta,y_alerta,'Alerta')
    
    
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
end
%% 
snr_atual
media_RR_s
Media_bpm
variancia_RR
intervalo_QS_ms
n_alertas

figure()
plot(snr_atual,intervalo_QS_ms)
title('intervalo QS')
xlabel('snr')
ylabel('ms')
grid on;