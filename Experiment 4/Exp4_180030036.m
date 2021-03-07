%% QPSK MODULATION AND DEMODULATION
% This is the experiment number 4 of Communications Lab course.

% In this first section an example of QPSK modulation and demodulation is
% demonstrated with randomly generated bits (using an SNR of 5 dBW).

% ############### GENERATING MESSAGE ###############

% We choose a length of our random signal to be transmitted
signal_length = 10000; % Chosen to be a multiple of 2

% A random signal of binary digits (0 and 1) is generated
signal_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(signal_bits(1:20));

% ############### QPSK DIGITAL MODULATION ###############

% We take 2 bits from the message signal at a time and match it to a QPSK
% symbol. The mapping can be done in two ways:

% 1. With Gray Labeling
% The mapping of bits is done as follows:
% 00 -> (1 + i)/sqrt(2)
% 01 -> (-1 + i)/sqrt(2)
% 11 -> (-1 - i)/sqrt(2)
% 10 -> (1 - i)/sqrt(2)

% Second parameter is 1 => Modulation with Gray labeling
tx_qpsk_symbols_gray = bits_to_qpsk(signal_bits.', 1);

% 2. Without Gray Labeling
% The mapping of bits is done as follows:
% 00 -> (1 + i)/sqrt(2)
% 01 -> (-1 + i)/sqrt(2)
% 10 -> (-1 - i)/sqrt(2)
% 11 -> (1 - i)/sqrt(2)

% Second parameter is 0 => Modulation without Gray labeling
tx_qpsk_symbols_nogray = bits_to_qpsk(signal_bits.', 0);

disp('First 10 Transmitted Symbols: With Gray Labeling');
disp(tx_qpsk_symbols_gray(1:10).');

disp('First 10 Transmitted Symbols: Without Gray Labeling');
disp(tx_qpsk_symbols_nogray(1:10).');

% ############### ADDING WHITE GAUSSIAN NOISE ###############

% Defining a SNR
snr = 5;

% To get the received noisy signal for gray labeling
noisy_signal_gray = add_noise(tx_qpsk_symbols_gray, snr);

% To get the received noisy signal for without gray labeling
noisy_signal_nogray = add_noise(tx_qpsk_symbols_nogray, snr);

figure(1)
scatter(real(noisy_signal_gray(1:100)), imag(noisy_signal_gray(1:100)))
title('Received Noisy Symbols (First 100)');
xlabel('Re');
ylabel('Im');
grid on;

% ############### DEMODULATION ###############

% Calling the demodulation function for both cases
demod_signal_gray = demodulation(noisy_signal_gray);
demod_signal_nogray = demodulation(noisy_signal_nogray);

% Second parameter is 1 => with Gray labeling
rx_signal_bits_gray = qpsk_to_bits(demod_signal_gray, 1);

% Second parameter is 0 => without Gray labeling
rx_signal_bits_nogray = qpsk_to_bits(demod_signal_nogray, 0);

disp('First 10 Demodulated Symbols: With Gray Labeling');
disp(demod_signal_gray(1:10).');

disp('First 10 Demodulated Symbols: Without Gray Labeling');
disp(demod_signal_nogray(1:10).');

% ############### BIT ERROR RATE ###############

% BER in case of gray labeling
bitwise_difference = rx_signal_bits_gray - signal_bits.';
ber_gray = calculate_ber(bitwise_difference)/signal_length;

% BER in case of no gray labeling
bitwise_difference = rx_signal_bits_nogray - signal_bits.';
ber_nogray = calculate_ber(bitwise_difference)/signal_length;

disp('BER in case of Gray Labeling:');
disp(ber_gray);

disp('BER in case of no Gray Labeling:');
disp(ber_nogray);

%% BER CALCULATION WITH GRAY LABELING

% In this second section we calculate the bit error rate with gray labeling
% We have an error resolution of 10^-7
% Please reduce the number of iterations to get lower run time.

% Defining the signal length
signal_length = 100000;

% Defining the number of iterations the signal is tested for each SNR
iterations = 100;

% We define SNR array in dBW
snr_array = 0:0.2:15;

%Vector to store bit error rate corresponding to each SNR value
bit_error_gray = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        signal_bits = randi([0,1],signal_length,1);

        % Generating QPSK symbols with gray labeling
        tx_qpsk_symbols_gray = bits_to_qpsk(signal_bits.', 1);
        
        % Adding noise and demodulating for each SNR
        noisy_signal_gray = add_noise(tx_qpsk_symbols_gray, snr_array(i));
        
        % Demodulation
        demod_signal_gray = demodulation(noisy_signal_gray);

        % To get the binary bits corresponding to demodulated symbols
        rx_signal_bits_gray = qpsk_to_bits(demod_signal_gray, 1);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = rx_signal_bits_gray - signal_bits.';

        bit_error_gray(i) = bit_error_gray(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error_gray = bit_error_gray/(signal_length*iterations);

% Plotting the BER vs. SNR graph
figure(2);
semilogy(snr_array, bit_error_gray, '-bo')
title('BER vs. SNR Plot (With Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% BER CALCULATION WITHOUT GRAY LABELING

% In this second section we calculate the BER without gray labeling
% We have an error resolution of 10^-7
% Please reduce the number of iterations to get lower run time.

% Defining the signal length
signal_length = 100000;

% Defining the number of iterations the signal is tested for each SNR
iterations = 100;

% We define SNR array in dBW
snr_array = 0:0.2:15;

%Vector to store bit error rate corresponding to each SNR value
bit_error_nogray = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        signal_bits = randi([0,1],signal_length,1);

        % Generating QPSK symbols with gray labeling
        tx_qpsk_symbols_nogray = bits_to_qpsk(signal_bits.', 0);
        
        % Adding noise and demodulating for each SNR
        noisy_signal_nogray = add_noise(tx_qpsk_symbols_nogray, snr_array(i));
        
        % Demodulation
        demod_signal_nogray = demodulation(noisy_signal_nogray);

        % To get the binary bits corresponding to demodulated symbols
        rx_signal_bits_nogray = qpsk_to_bits(demod_signal_nogray, 0);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = rx_signal_bits_nogray - signal_bits.';

        bit_error_nogray(i) = bit_error_nogray(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error_nogray = bit_error_nogray/(signal_length*iterations);

% Plotting the BER vs. SNR graph
figure(3);
semilogy(snr_array, bit_error_nogray, '-bo')
title('BER vs. SNR Plot (Without Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% THEORETICAL GRAPH USING Q-FUNCTION

% In this third section we first compute a theoretical curves defining the
% error probability variation with SNR for a QPSK communication system.
% We then plot the theoretical curves with our practical curves.

% We define signal to noise ratio as:
% SNR = Eb/No

% In our case of QPSK communication system, we have
% Eb = 1/2
% Also, the textbook defines
% sigma = sqrt(No/2)        (sigma => standard deviation of noise)

% So, we get the relation between sigma and SNR as:
% sigma = 1/sqrt(4*SNR)

% We use the result:
% P_e = Q(d_ij/2*sigma)
% d_ij = sqrt(2) or 2.....(distance between ith and jth symbols)

% ########## WITH GRAY LABELING ##########

% On analysing the probability of error with gray labeling, we get:
% P_e = Q(sqrt(2*SNR))

ber_theo_gray = zeros(1, numel(snr_array));

% Calculation of theoretical error probability using Q-function
for i = 1:numel(snr_array)
    ber_theo_gray(i) = qfunc(sqrt(2*(10^(snr_array(i)/10))));
end

% ########## WITHOUT GRAY LABELING ##########

% On analysing the probability of error without gray labeling, we get:
% P_e = (3/2)*Q(sqrt(2*SNR)) - [Q(sqrt(2*SNR))]^2

ber_theo_nogray = zeros(1, numel(snr_array));

% Calculation of theoretical error probability using Q-function
for i = 1:numel(snr_array)
    q_value = qfunc(sqrt(2*(10^(snr_array(i)/10))));
    ber_theo_nogray(i) = 1.5*q_value - q_value*q_value;
end

figure(4);
semilogy(snr_array, ber_theo_gray, '-r')
hold on
semilogy(snr_array, ber_theo_nogray, '-b')
hold off
title('Theoretical BER vs. SNR Plot (With/Without Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
legend('With Gray Labeling', 'Without Gray Labeling');
legend
grid on;

figure(5);
semilogy(snr_array, ber_theo_gray, '-rx')
hold on
semilogy(snr_array, bit_error_gray, '-bo')
hold off
title('BER vs. SNR Plot (With Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
legend('Theoretical (Gray Labeling)', 'Practical (Gray Labeling)');
legend
grid on;

figure(6);
semilogy(snr_array, ber_theo_nogray, '-rx')
hold on
semilogy(snr_array, bit_error_nogray, '-bo')
hold off
title('BER vs. SNR Plot (Without Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
legend('Theoretical (No Gray Labeling)', 'Practical (No Gray Labeling)');
legend
grid on;

figure(7);
semilogy(snr_array, ber_theo_gray, '-rx')
hold on
semilogy(snr_array, bit_error_gray, '-co')

semilogy(snr_array, ber_theo_nogray, '-bx')
semilogy(snr_array, bit_error_nogray, '-go')
hold off
title('BER vs. SNR Plot (With/Without Gray Labeling, Theoretical and Practical)');
xlabel('SNR (dBW)');
ylabel('BER');
xlim([0 12]);
legend('Theoretical (Gray Labeling)', 'Practical (Gray Labeling)', 'Theoretical (No Gray Labeling)', 'Practical (No Gray Labeling)');
legend
grid on;

%% FUNCTIONS

% Function for QPSK modulation with/without gray labeling
function tx_qpsk_symbols = bits_to_qpsk(signal_bits, gray_label)
    
    tx_qpsk_symbols = zeros(1, numel(signal_bits)/2);
    symbols = [1+1i, -1+1i, -1-1i, 1-1i]/sqrt(2);

    if gray_label == 1 % With Gray Labeling
        for i = 1:2:numel(signal_bits)-1
            if isequal(signal_bits(i:i+1), [0 0])
                tx_qpsk_symbols(ceil(i/2)) = symbols(1);

            elseif isequal(signal_bits(i:i+1), [0 1])
                tx_qpsk_symbols(ceil(i/2)) = symbols(2);

            elseif isequal(signal_bits(i:i+1), [1 1])
                tx_qpsk_symbols(ceil(i/2)) = symbols(3);

            elseif isequal(signal_bits(i:i+1), [1 0])
                tx_qpsk_symbols(ceil(i/2)) = symbols(4);
            end
        end
    elseif gray_label == 0 % Without Gray Labeling
        for i = 1:2:numel(signal_bits)-1
            if isequal(signal_bits(i:i+1), [0 0])
                tx_qpsk_symbols(ceil(i/2)) = symbols(1);

            elseif isequal(signal_bits(i:i+1), [0 1])
                tx_qpsk_symbols(ceil(i/2)) = symbols(2);

            elseif isequal(signal_bits(i:i+1), [1 1])
                tx_qpsk_symbols(ceil(i/2)) = symbols(4);

            elseif isequal(signal_bits(i:i+1), [1 0])
                tx_qpsk_symbols(ceil(i/2)) = symbols(3);
            end
        end
    else
        disp('Invalid Input for Gray Labeling');
    end
end

% Function to add complex white gaussian noise to a signal
function noisy_signal = add_noise(tx_signal, snr)
    % Generating a random complex number with normal distribution
    Re_noise = sqrt(1/(4*(10^(snr/10))))*randn(1,numel(tx_signal));
    Im_noise = sqrt(1/(4*(10^(snr/10))))*randn(1,numel(tx_signal));
    noise = Re_noise + 1i*Im_noise;
    
    noisy_signal = tx_signal + noise;
end

function demod_signal = demodulation(noisy_signal)
    
    demod_signal = zeros(1, numel(noisy_signal));
    symbols = [1+1i, -1+1i, -1-1i, 1-1i]/sqrt(2);
    
    for i = 1:numel(noisy_signal)
        % To get distance of signal point from each of the symbols
        distance = abs(symbols - noisy_signal(i));
        
        % To get the index of the minimum distance
        [~, index] = min(distance);
        
        % The closest (distance wise) symbol is chosen
        demod_signal(i) = symbols(index);
    end
end

function signal_bits = qpsk_to_bits(signal, gray_label)
    
    signal_bits = zeros(1, numel(signal)*2);
    symbols = [1+1i, -1+1i, -1-1i, 1-1i]/sqrt(2);

    if gray_label == 1 % With gray labeling
        for i = 1:2:2*numel(signal)-1
            if isequal(signal(ceil(i/2)), symbols(1))
                signal_bits(i:i+1) = [0 0]; 
                
            elseif isequal(signal(ceil(i/2)), symbols(2))
                signal_bits(i:i+1) = [0 1];
                
            elseif isequal(signal(ceil(i/2)), symbols(3))
                signal_bits(i:i+1) = [1 1];
                
            elseif isequal(signal(ceil(i/2)), symbols(4))
                signal_bits(i:i+1) = [1 0];
            end
        end
    elseif gray_label == 0 % Without gray labeling
        for i = 1:2:2*numel(signal)-1
            if isequal(signal(ceil(i/2)), symbols(1))
                signal_bits(i:i+1) = [0 0]; 
                
            elseif isequal(signal(ceil(i/2)), symbols(2))
                signal_bits(i:i+1) = [0 1];
                
            elseif isequal(signal(ceil(i/2)), symbols(3))
                signal_bits(i:i+1) = [1 0];
                
            elseif isequal(signal(ceil(i/2)), symbols(4))
                signal_bits(i:i+1) = [1 1];
            end
        end
    else
        disp('Invalid Input for Gray Labeling');
    end
end

function bit_error = calculate_ber(bitwise_difference)
    bit_error = 0;
    for i=1:numel(bitwise_difference)
        if bitwise_difference(i) ~= 0
            bit_error = bit_error + 1;
        end
    end
end