%% GENERATING MESSAGE SIGNAL

% We define message amplitude and frequency: Am and fm
Am = 1;
fm = 50;

fc = 1000; % Carrier frequency
fs = 10*fc; % Sampling frequency
cycles = 4;

t = 0:1/fs:(cycles/fm - 1/fs);
 
% Defining message signal
msg_signal = Am*(cos(2*pi*fm*t) + sin(2*pi*4*fm*t));

% Getting Fourier transform of the message signal
msg_signal_freq = fftshift(fft(msg_signal))/numel(msg_signal);

Nfft = numel(msg_signal_freq);
f = -fs/2: fs/Nfft :(fs/2 - fs/Nfft);

% Plotting the message signal in time and frequency domain
figure(1)
subplot(2,1,1)
plot(t,msg_signal,'b');
title('Message Signal');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2)
plot(f, abs(msg_signal_freq));
title('Message Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%% CONVENTIONAL AM OF MESSAGE SIGNAL

A = 1; % a constant value
a_mod = 0.5; % modulation index
Ac = A*abs(min(msg_signal))/a_mod; % carrier wave amplitude

% A general carrier for the signals
carrier_wave = cos(2*pi*fc*t);       

% Plotting the carreir wave in time and frequency domain
figure(2)
subplot(2,1,1)
plot(t,carrier_wave,'b');
title('Carrier wave in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2)
plot(f, abs(fftshift(fft(carrier_wave))/numel(carrier_wave)));
title('Carrier wave in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

% Generating passband signal for conventional AM
am_passband_signal = A*(msg_signal.*carrier_wave) + Ac*carrier_wave;

% Getting the Fourier transform of the passband signal
am_passband_freq = fftshift(fft(am_passband_signal))/numel(am_passband_signal);
%f = -fs/2:(fs/2 - 1);

figure(3)
subplot(2,1,1);
hold on
plot(t,msg_signal,'r');
plot(t,am_passband_signal, 'b');
hold off
title('Conventional AM Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');
legend('Message signal','Conventional AM passband signal')

subplot(2,1,2);
plot(f, abs(am_passband_freq));
title('Conventional AM Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');

%% DEMODULATION FOR CONVENTIONAL AM

temp_signal = am_passband_signal.*(2*carrier_wave);
temp_signal = temp_signal - Ac;

temp_signal_freq = fftshift(fft(temp_signal))/numel(temp_signal);

[b, a] = butter(10, fc/(fs/2));

demod_signal = filter(b, a, temp_signal);
demod_signal_freq = fftshift(fft(demod_signal))/numel(demod_signal);

figure(4)
subplot(2,1,1)
plot(demod_signal);
title('Conventional AM Demodulated Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2)
plot(f, abs(demod_signal_freq));
title('Conventional AM Demodulated Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');

%% DSB-SC MODULATION OF MESSAGE SIGNAL

% Generating the passband signal
dsbsc_passband_signal = A*msg_signal.*carrier_wave;  

% Getting the Fourier Transform
dsbsc_passband_freq = fftshift(fft(dsbsc_passband_signal))/numel(dsbsc_passband_signal);

% Plotting the modulated signal along with message signal
figure(5)
subplot(2,1,1);
hold on
plot(t,msg_signal,'r');
plot(t,dsbsc_passband_signal, 'b');
hold off
title('DSB-SC Modulated Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');
legend('Message signal','DSB-SC Passband Signal')

subplot(2,1,2);
plot(f, abs(dsbsc_passband_freq));
title('DSB-SC Modulated Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');

%% DEMODULATION FOR DSB-SC

% Multiplying the received passband signal with 2 times the carrier wave
dsbsc_temp_signal = dsbsc_passband_signal.*(2*carrier_wave);

% Applying the low pass Butterworth filter generated above
dsbsc_demod_signal = filter(b, a, dsbsc_temp_signal);

% Getting Fourier transform of demodulated signal
dsbsc_demod_signal_freq = fftshift(fft(dsbsc_demod_signal))/numel(dsbsc_demod_signal);

figure(6)
subplot(2,1,1)
plot(dsbsc_demod_signal);
title('DSB-SC Demodulated Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2)
plot(f, abs(dsbsc_demod_signal_freq));
title('DSB-SC Demodulated Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');

%% SSB-SC MODULATION OF MESSAGE SIGNAL

% We use the upper sideband (USB) signal and remove the lower sideband
% (LSB) signal using the Hilbert transform. We can do vice-versa as well.
% USB signal = I - Q
% LSB signal = I + Q
% I => inphase component and Q => quadrature component

% Generating the inphase and quadrature components of the carrier wave
carrier_wave_inphase = cos(2*pi*fc*t);
carrier_wave_quadrature = sin(2*pi*fc*t);

% Generating the Hilbert transform of message signal
hilbert_msg_signal = imag(hilbert(msg_signal));

% Modulating the message signal using SSB-SC technique
ssbsc_passband_signal = msg_signal.*carrier_wave_inphase - hilbert_msg_signal.*carrier_wave_quadrature;

% Getting the Fourier transform of modulated signal
ssbsc_passband_freq = fftshift(fft(ssbsc_passband_signal))/numel(ssbsc_passband_signal);

% Plotting the modulated message signal in time and frequency domain
figure(7)
subplot(2,1,1);
hold on
plot(t,msg_signal,'r');
plot(t,ssbsc_passband_signal, 'b');
hold off
title('SSB-SC Modulated Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');
legend('Message signal','SSB-SC Passband Signal')

subplot(2,1,2);
plot(f, abs(ssbsc_passband_freq));
title('SSB-SC Modulated Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');

%% DEMODULATION FOR SSB-SC

% Multiplying the passband signal by 2cos(2*pi*fc*t)
ssbsc_temp_signal = ssbsc_passband_signal.*(2*carrier_wave_inphase);

% Applying the low pass Butterworth filter generated above
ssbsc_demod_signal = filter(b, a, ssbsc_temp_signal);

% Getting Fourier transform of demodulated signal
ssbsc_demod_signal_freq = fftshift(fft(ssbsc_demod_signal))/numel(ssbsc_demod_signal);

% Plotting the demodulated signal in time and frequency domain
figure(8)
subplot(2,1,1)
plot(ssbsc_demod_signal);
title('SSB-SC Demodulated Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2)
plot(f, abs(ssbsc_demod_signal_freq));
title('SSB-SC Demodulated Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');
