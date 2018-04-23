%Calculating Beat Rate of heart Using Discrete Wavelet Transform 
%Databases are available in:
%http://www.physionet.org
%you can  convert ECG signal by using PhysioNet_ECG_Exporter_2.m
%If you have any question mail me at:
%ehsan.mirrahimi@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;
load 100m.mat
%%%Eliminate Baseline Drift
%s1=ECG_1;s2=smooth(s1,150);ecgsmooth=s1-s2;
ecgsmooth=val;
%%%apply Wavelet Transform
[C,L]=wavedec(ecgsmooth,8,'db4');
[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
% %%%Denoise 
 [thr,sorh,keepapp]=ddencmp('den','wv',ecgsmooth);
 cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);


 %%%thresholding1
% max_value=max(cleanecg);
% mean_value=mean(cleanecg);
% threshold=(max_value-mean_value)/2;
% %%%R detection algorithm
% a5=appcoef(C,L,'db4',5);
% C1=[a5;d5;d4;d3];


% L1=[length(a5);length(d5);length(d4);length(d3);length(cleanecg)];
% R_detect_signal=waverec(C1,L1,'db4');
% R_detect_squared=R_detect_signal.^2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Beat_Rate_Extraction_Algorithm
% for a=1:length(R_detect_squared)
%     if R_detect_squared(a)>threshold
%         R_detect_new(a)=R_detect_squared(a);
%     else
%         R_detect_new(a)=0;
%     end
% end
% mean_R_detect=5*mean(R_detect_new);
% for q=1:length( R_detect_new)-1
%     if  R_detect_new(q)< mean_R_detect
%         R_detect_new(q)=0;
%     end
% 
% end
% %%%%%%%%%%%%%%%%%%
% d=0;
% for b=1:length( R_detect_new)-1
%         if ( R_detect_new(b)==0) & ( R_detect_new(b+1)~=0)
%         d=d+1;
%         indext(d)= b+1;
%         end
% end
% fs_R_deetect=length(R_detect_new)/20;
% time=indext.*1/fs_R_deetect;
% ind=0;
% for z=1:length(time)-1
%     ind=ind+1;
%     time_diff(ind)=time(z+1)-time(z);
% end
% av_time=mean(time_diff);
% Square_Number=av_time/.2;
% beat_Rate=300/Square_Number;
% high=max(R_detect_new);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Plot the Orginal Signal and Eliminating Baseline Drift signal
% subplot(411);plot(s1);title('Orginal Signal');
% subplot(412);plot(s1-s2);title('Baseline drift Elimination');
% subplot(413);plot(cleanecg);title('Main Signal');
% subplot(414);plot(R_detect_new);title('R detected Signal');
% text(length(R_detect_new)/2,high,['Beat Rate = ',num2str(fix(beat_Rate))],'EdgeColor','red');