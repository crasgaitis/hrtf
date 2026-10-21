%% a_preview
% the purpose of this file is to double check the file structure of the
% provided hrtf data from the ito, 2022 paper and analyze ild, itd, and spectral
% cues.

%% check file structure

data = 'D:\Data\Catherine\hrtf\2018-08-15_1000samp_ILDITDmod_correct_noitd_ex_ILD.mat';

d = dir(data);
fprintf('FILE: %s (%.0f MB)\n\n', d.name, d.bytes/1e6);

m = matfile(data);
whos(m);

%% set up

% check allimps_percep_ILD_corrected
% 2 10 200 2001
% ear, elevation, azimuth, hrir (in time samps)

H = m.allimps_percep_ILDcorrected;
fs = 500e3;
[nEar, nEl, nAz, nT] = size(H);

nfft = 4096;
f    = (0:nfft/2) * fs / nfft;
band = f >= 5e3 & f <= 80e3;
magdB = @(x) 20*log10(abs(fft(squeeze(x), nfft)) + eps);
half  = @(s) s(1:nfft/2+1);

%% one dir example
e = 1;  a = 50;                          % elevation index, azimuth index
t = ((1:nT) - ceil(nT/2)) / fs * 1e3;    % ms, centered
figure('Name','One direction');
subplot(2,1,1);
plot(t, squeeze(H(1,e,a,:)), t, squeeze(H(2,e,a,:)));
xlim([-0.5 1]); xlabel('Time (ms)'); ylabel('Amplitude');
legend('ear 1','ear 2'); title(sprintf('HRIR  (elev idx %d, az idx %d)', e, a));
subplot(2,1,2);
plot(f/1e3, half(magdB(H(1,e,a,:))), f/1e3, half(magdB(H(2,e,a,:))));
xlim([0 100]); grid on; xlabel('Frequency (kHz)'); ylabel('dB');
legend('ear 1','ear 2'); title('HRTF magnitude');

%% ild v azimuth
ild = zeros(1, nAz);
for a = 1:nAz
    s1 = half(magdB(H(1,1,a,:)));  s2 = half(magdB(H(2,1,a,:)));
    ild(a) = mean(s1(band)) - mean(s2(band));
end
figure('Name','ILD vs azimuth');
plot(1:nAz, ild, 'o-'); grid on; yline(0,'k--');
xlabel('Azimuth index'); ylabel('ILD (dB), ear1 − ear2');
title('ILD vs azimuth, elevation 0  (zero crossings ≈ front/back)');

%% itd v azimuth
itd_xc = zeros(1, nAz);
for a = 1:nAz
    h1 = squeeze(H(1,1,a,:));  h2 = squeeze(H(2,1,a,:));
    [c, lags] = xcorr(h1, h2, 50);          % ±50 samples = ±100 µs
    [~, i] = max(c);
    itd_xc(a) = lags(i) / fs * 1e6;         % µs
end

pt = m.peaktimes;                            % 2 x 10 x 200
itd_pt = squeeze(pt(1,1,:) - pt(2,1,:))';   

az = (0:nAz-1) * 1.8;                        % assumed azimuth in degrees
figure;
yyaxis left;  plot(az, itd_xc, 'o-'); ylabel('ITD, xcorr (µs)');
yyaxis right; plot(az, itd_pt, '.-');  ylabel('peaktimes diff (file units)');
xlabel('Azimuth (deg, assumed)'); grid on; xlim([0 360]);
title('ITD vs azimuth, elevation 0');

%% spectral map
S = zeros(nnz(band), nAz);
for a = 1:nAz
    s = half(magdB(H(1,1,a,:)));  S(:,a) = s(band);
end
figure('Name','Spectral map');
imagesc(1:nAz, f(band)/1e3, S); axis xy; colorbar;
xlabel('Azimuth index'); ylabel('Frequency (kHz)');
title('Ear 1 HRTF (dB) vs azimuth, elevation 0  (notches = spectral cues)');
