%% GENERATING AN ANALOG MESSAGE SIGNAL
% This is the experiment number 2 of Communications Lab course.

% Define the amplitude and frequency for message signal
Am = 4;
Fm = 100;

% Define a sampling frequency
Fs = 100000;

% Define number of cycles and time vector
ncycles = 4;
t = 0: 1/Fs : (ncycles/Fm - 1/Fs);

% Generating a message signal (uncomment the message to be used)
msg_signal = Am*(cos(2*pi*Fm*t) + sin(2*pi*2*Fm*t)); % sinusoid with 2 frequecncies
% msg_signal = Am*cos(2*pi*Fm*t) % sinusoid with 1 frequency

% Getting the length of message signal
msg_length = numel(msg_signal);

% Getting the Fourier transform of the message signal
msg_signal_freq = fftshift(fft(msg_signal))/msg_length;

% Define frequency vector
f = -Fs/2 : Fs/msg_length : (Fs/2 - Fs/msg_length);

figure(1);
subplot(211);
plot(t, msg_signal);
title('Message Signal in Time Domain');
xlabel('Time');
ylabel('Magnitude');

subplot(212);
plot(f, abs(msg_signal_freq));
title('Message Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Magnitude');
xlim([-5000 5000]);

%% FREQUENCY MODULATION

% Generating the modulated phase shift for the passband signal
theta = zeros(1, msg_length);

% Define a constant
kf = 100;

% Define the magnitude and frequency of carrier wave for passband signal
Ac = 5;
Fc = 2000;

for i = 1:msg_length
    % Numerically integrating the msg_signal, that is, summation
    theta(i) = 2*pi*kf*sum(msg_signal(1:i))/(Fs);
end

% Generating the passband signal
passband_signal = Ac*cos(2*pi*Fc*t + theta);

% Getting the Fourier Transform of the passband signal
passband_signal_freq = fftshift(fft(passband_signal))/numel(passband_signal);

figure(2);
subplot(211);
plot(t, passband_signal);
title('Frequency Modulated Passband Signal in Time Domain');
xlabel('Time');
ylabel('Magnitude');

subplot(212);
plot(f, abs(passband_signal_freq));
title('Frequency Modulated Passband Signal in Frequency Domain');
xlabel('Time');
ylabel('Magnitude');
xlim([-5000 5000]);

%% FREQUENCY DEMODULATION

% Differentiating the passband signal
derivative = diff(passband_signal)*Fs;

% Appending a zero to remove length mismatch
derivative = [derivative 0];

figure(3);
plot(t, derivative);
title('Derivative of Passband Signal');
xlabel('Time');
ylabel('Magnitude');

% Detecting the upper envelope using Hilbert transform
envelope = abs(hilbert(derivative));

% Removing scaling factors and DC value
demod_signal = (envelope/(2*pi) - Ac*Fc)/(kf*Ac);

% Getting the Fourier transform of demodulated signal
demod_signal_freq = fftshift(fft(demod_signal))/numel(demod_signal);

figure(4);
subplot(211);
plot(t, demod_signal);
title('Demodulated Signal in Time Domain');
xlabel('Time');
ylabel('Magnitude');

subplot(212);
plot(f, abs(demod_signal_freq));
title('Demodulated Signal in Frequency Domain');
xlabel('Frequency');
ylabel('Magnitude');
xlim([-5000 5000]);