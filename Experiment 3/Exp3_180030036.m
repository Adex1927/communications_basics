%% TRANSMISSION OF A RANDOM SIGNAL USING BPSK
% This is the experiment number 3 of Communications Lab course.

% In this first section an example of BPSK modulation and demodulation is
% demonstrated with randomly generated bits (using an SNR of 5 dBW).

% We choose a length of our random signal to be transmitted
signal_length = 10000;

% We choose 20 samples for generating pulse
samples = 20;

% ############### BPSK DIGITAL MODULATION ###############

%Generating a random signal with binary bits (0 and 1)
signal_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(signal_bits(1:20));

% Mapping the signal bits to the BPSK symbols
% 0 -> -1
% 1 -> 1
% Generating BPSK symbols from the bits
tx_signal_symbols = bits_to_symbols(signal_bits);

[tx_symbols_pulse, time] = pulse_shaping(tx_signal_symbols, samples);

figure(1)
subplot(311)
plot(time(1:2000),tx_symbols_pulse(1:2000))
title('Transmitted Symbols Pulse (First 200)');
xlabel('Time/T');
ylabel('BPSK Symbols');

% ############### ADDING WHITE GAUSSIAN NOISE ###############

% Definition a SNR
snr = 5;

% To get the received noisy signal
noisy_signal = add_noise(tx_symbols_pulse, snr);

subplot(312)
plot(time(1:2000),noisy_signal(1:2000))
title('Received Noisy Symbols Pulse (First 200)');
xlabel('Time/T');
ylabel('Noisy BPSK Symbols');

% ############### DEMODULATION ###############

rx_signal_symbols = rxsignal_to_symbols(noisy_signal, samples);

%Generating a pulse from the demodulated symbols
[rx_symbols_pulse, time] = pulse_shaping(rx_signal_symbols, samples);

subplot(313)
plot(time(1:2000),rx_symbols_pulse(1:2000))
title('Demodulated Symbols Pulse (First 200)');
xlabel('Time/T');
ylabel('Noisy BPSK Symbols');

%% BIT ERROR RATE CALCULATION
% In this second section we calculate the bit error rate of our system.
% We have an error resolution of 10^-8
% Please reduce the number of iterations to get lower run time.

% Defining the signal length
signal_length = 100000;

% Defining number of samples for making pulse
samples = 1;

% Defining the number of iterations the signal is tested for each SNR
iterations = 1000;

% We define SNR array in dBW
snr_array = 0:0.2:15;

%Vector to store bit error rate corresponding to each SNR value
bit_error = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        tx_signal_bits = (randi([0,1],signal_length,1))';

        % Generating BPSK symbols from the bits
        tx_signal_symbols = bits_to_symbols(tx_signal_bits);

        % Generating a pulse from the symbols
        [tx_symbols_pulse, ~] = pulse_shaping(tx_signal_symbols, samples);
        
        % Adding noise and demodulating for each SNR
        noisy_signal = add_noise(tx_symbols_pulse, snr_array(i));
        rx_signal_symbols = rxsignal_to_symbols(noisy_signal, samples);

        % To get the binary bits corresponding to demodulated symbols
        rx_signal_bits = symbols_to_bits(rx_signal_symbols);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = tx_signal_bits - rx_signal_bits;

        bit_error(i) = bit_error(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error = bit_error/(signal_length*iterations);

% Plotting the BER vs. SNR graph
figure(2);
semilogy(snr_array, bit_error, '-bo')
title('BER vs. SNR Plot');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% THEORETICAL GRAPH USING Q-FUNCTION
% In this third section we first compute a theoretical curve defining the
% error probability variation with SNR for a BPSK communication system.
% We then plot this theoretical curve with our practical curve.

% According to the textbook, the error probability is given as:
% P_e = Q(sqrt(2Eb/No))

% We define signal to noise ratio as:
% SNR = Eb/No

% In our case of BPSK communication system, we have
% Eb = 1
% Also, the textbook defines
% sigma = sqrt(No/2)        (sigma => standard deviation of noise)

% So, we get the relation between sigma and SNR as:
% sigma = 1/sqrt(2*SNR)

ber_theoretical = zeros(1, numel(snr_array));

% Calculation of theoretical error probability using Q-function
for i = 1:numel(snr_array)
    ber_theoretical(i) = qfunc(sqrt(2*(10^(snr_array(i)/10))));
end

figure(3);
semilogy(snr_array, bit_error, '-bo')
hold on
semilogy(snr_array, ber_theoretical, '-rx')
hold off
title('BER vs. SNR Plot');
xlabel('SNR (dBW)');
ylabel('BER');
legend('Experimental curve', 'Theoretical curve using Q-function');
legend
grid on;

%% FUNCTIONS

function signal_symbols = bits_to_symbols(signal_bits)
    signal_symbols = signal_bits*2 - 1;
end

function signal_bits = symbols_to_bits(signal_symbols)
    signal_bits = (signal_symbols + 1)/2;
end

%Function to make a pulse, given symbols and samples for each symbol
function [tx_symbols_pulse, time] = pulse_shaping(tx_symbols, samples)
    tx_symbols_pulse = [];
    time = [];
    time_period = 1/samples;
    for i = 1:numel(tx_symbols)
        for j = 1:samples
            if tx_symbols(i) == 1
                tx_symbols_pulse((i-1)*samples+j) = 1;
            else
                tx_symbols_pulse((i-1)*samples+j) = -1;
            end
            time((i-1)*samples+j) = ((i-1)*samples+j)*time_period;
        end
    end
end

% Function to add noise to the given signal
function noisy_signal = add_noise(tx_signal, snr)
    % Generating a random number with normal distribution
    noise = sqrt(1/(2*(10^(snr/10))))*randn(1,numel(tx_signal));
    
    noisy_signal = tx_signal + noise;
end

%Function to get BPSK symbols from the received noisy signal
function rx_symbols = rxsignal_to_symbols(rx_signal_noisy, samples)
    rx_symbols = [];
    for i = 1:(numel(rx_signal_noisy)/samples)
        count1 = 0;
        for j = 1:samples
            if rx_signal_noisy((i-1)*samples+j) > 0
                count1 = count1+1;
            end
        end
        
        if count1 > samples/2
            rx_symbols(i) = 1;
        else
            rx_symbols(i) = -1;
        end
    end
end

%Function to calculate the number of erroneous bits
function bit_error = calculate_ber(bitwise_difference)
    bit_error = 0;
    for i=1:numel(bitwise_difference)
        if bitwise_difference(i) ~= 0
            bit_error = bit_error + 1;
        end
    end
end