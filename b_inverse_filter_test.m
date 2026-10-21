%% b_inverse_filter
% the purpose of this file is to build and validate an inverse filter to 
% account for the earphone delivery path (DAC -> ED1 -> EC1 (potentially with horn 
% -> ear canal), so sound reaching eardrum matches the target HRTFs 
% (from the allimps_percep_ILDcorrected var).

% we use the left ear only as an example.

%% setup
data = 'D:\Data\Catherine\hrtf\2018-08-15_1000samp_ILDITDmod_correct_noitd_ex_ILD.mat';
m  = matfile(data);
fs = 500e3;

% get earphone response
figure;
subplot(1, 2, 1);
plot(m.f_l, m.fr_l, 'b'); % left ear
subplot(1, 2, 2);
plot(m.f_r, m.fr_r, 'r'); % right ear

% convert dB to linear
mag_l = 10.^(m.fr_l/20);

%% inverse filter
N  = 2048; % 2048
gL = make_inverse_filter(m.f_l, mag_l, fs, N, [5e3 80e3], 25);

%% apply & validation
% ito provides allimps_percep and allimps_present (_percep but put thru their
% earphone inverse filter). to validate this code, we need to transform _percep
% match _present.

P = m.allimps_percep;  Q = m.allimps_present;
e = 1; a = 10; ear = 1;        % m.leftid
hP = squeeze(P(ear,e,a,:));
hQ = squeeze(Q(ear,e,a,:));
hMine = conv(hP, gL, 'same');

nfft = 4096; f = (0:nfft/2)*fs/nfft;
dB = @(x) 20*log10(abs(fft(x, nfft)) + eps);
sP = dB(hP); sQ = dB(hQ); sM = dB(hMine);
k = 1:nfft/2+1;
% figure;
% plot(f/1e3, sP(k) - mean(sP(k)), f/1e3, sQ(k) - mean(sQ(k)), f/1e3, sM(k) - mean(sM(k)), '--');
% xlim([0 100]); grid on; xlabel('kHz'); ylabel('dB (mean removed)');

inb = f >= 5e3 & f <= 80e3;
nP = sP(k) - mean(sP(inb));
nQ = sQ(k) - mean(sQ(inb));
nM = sM(k) - mean(sM(inb));

figure;
plot(f/1e3, nP, f/1e3, nQ, f/1e3, nM, '--');
xlim([5 80]); grid on; xlabel('kHz'); ylabel('dB (band mean removed)');
legend('percep','present (theirs)','percep * my inverse');

%% qualification
inb = f >= 5e3 & f <= 80e3;
d = (sQ(k) - mean(sQ(k))) - (sM(k) - mean(sM(k)));
figure; plot(f(inb)/1e3, d(inb)); grid on;
xlabel('kHz'); ylabel('dB'); title('present (theirs) − mine, 5–80 kHz');
fprintf('RMS difference in band: %.1f dB\n', rms(d(inb)));
fprintf('RMS difference (band mean removed), 5–80 kHz: %.2f dB\n', rms(nQ(inb) - nM(inb)));



