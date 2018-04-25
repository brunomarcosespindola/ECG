clc; clear all;
q = load('16265m');

EKG1 = q.val(1,:);
EKG2 = q.val(1,:);
Fs   = 128;
t    = [0:size(EKG1,2)-1]/Fs;
[R1,TR1]  = findpeaks( EKG1, t, 'MinPeakHeight',200);
[Q1,TQ1]  = findpeaks(-EKG1, t, 'MinPeakHeight',100);              % NOTE: No ‘S’ Waves In EKG1
[R2,TR2]   = findpeaks( EKG2, t, 'MinPeakHeight', 50);
[QS2,TQS2] = findpeaks(-EKG2, t, 'MinPeakHeight', 75);
figure(1)
subplot(2,1,1)
plot(t, EKG1)
hold on
plot(TR1, R1, '^r')
plot(TQ1, -Q1, 'vg')
hold off
grid
axis([0  2    ylim])
legend('EKG', 'R', 'Q')
subplot(2,1,2)
plot(t, EKG2)
hold on
plot(TR2, R2, '^r')
plot(TQS2(1:2:end), -QS2(1:2:end), 'vg')
plot(TQS2(2:2:end), -QS2(2:2:end), 'vb')
hold off
grid
axis([0  10    ylim])
legend('EKG', 'R', 'Q', 'S')