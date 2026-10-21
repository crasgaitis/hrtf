function g = make_inverse_filter(f_meas, mag_meas, fs, N, band, maxBoost_dB)
% Zero-phase (linear-phase after centering) inverse FIR filter.
%   f_meas, mag_meas : measured path, frequency (Hz) and LINEAR magnitude
%   fs               : sampling rate of the filter (Hz)
%   N                : filter length (even), e.g. 2048
%   band             : [fLow fHigh] range to correct, e.g. [5e3 80e3]
%   maxBoost_dB      : cap on boost, e.g. 20

    fgrid = (0:N/2) * fs / N;
    mag = interp1(f_meas(:)', mag_meas(:)', fgrid, 'linear', 'extrap');
    mag = max(mag, eps);

    inb = fgrid >= band(1) & fgrid <= band(2);
    mag = mag / median(mag(inb));                     % normalize: in-band ≈ 1

    G = zeros(size(fgrid));
    G(inb) = min(1 ./ mag(inb), 10^(maxBoost_dB/20)); % invert, cap boost

    % smooth half-octave tapers outside the band (no hard edges)
    lo = fgrid >= band(1)/sqrt(2) & fgrid < band(1);
    G(lo) = G(find(inb,1)) * 0.5*(1 - cos(pi*log2(fgrid(lo)/(band(1)/sqrt(2)))/0.5));
    hi = fgrid > band(2) & fgrid <= min(band(2)*sqrt(2), fs/2);
    G(hi) = G(find(inb,1,'last')) * 0.5*(1 + cos(pi*log2(fgrid(hi)/band(2))/0.5));

    Gfull = [G, G(end-1:-1:2)];                       % symmetric spectrum, zero phase
    g = real(ifft(Gfull));
    g = fftshift(g) .* hann(N)';                      % center (delay N/2) + window
    g = g(:);
end