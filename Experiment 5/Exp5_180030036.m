%% 16QAM MODULATION AND DEMODULATION
% This is the experiment number 5 of Communications Lab course.

% In this first section an example of 16QAM modulation and demodulation is
% demonstrated with randomly generated bits (using an SNR of 5 dBW).

% ############### GENERATING MESSAGE ###############

% We choose a length of our random signal to be transmitted
signal_length = 10000; % Chosen to be a multiple of 4

% A random signal of binary digits (0 and 1) is generated
signal_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(signal_bits(1:20));

% ############### 16QAM DIGITAL MODULATION ###############

% We take 4 bits from the message signal at a time and match it to a QPSK
% symbol. The mapping can be done in two ways:

% 1. With Gray Labeling
% The mapping of bits is done as follows:

%                         16QAM With Gray Labelling
%                                    |
%                   *          *    3|     *          *
%                  0010       0011   |    0001       0000
%                                    |
%                                    |
%                   *          *    1|     *          *
%                  0110       0111   |    0101       0100   
%           _________________________|___________________________
%                  -3         -1     |     1          3 
%                                    |
%                   *          *   -1|     *          *
%                  1110       1111   |    1101       1100
%                                    | 
%                                    |
%                   *          *   -3|     *          *
%                  1010       1011   |    1001       1000
%                                    |

% Second parameter is 1 => Modulation with Gray labeling
tx_16qam_symbols_gray = bits_to_16qam(signal_bits.', 1);

% 2. Without Gray Labeling
% The mapping of bits is done as follows:

%                         16QAM Without Gray Labelling
%                                    |
%                   *          *    3|     *          *
%                  0011       0010   |    0001       0000
%                                    |
%                                    |
%                   *          *    1|     *          *
%                  0111       0110   |    0101       0100   
%           _________________________|___________________________
%                  -3         -1     |     1          3 
%                                    |
%                   *          *   -1|     *          *
%                  1011       1010   |    1001       1000
%                                    | 
%                                    |
%                   *          *   -3|     *          *
%                  1111       1110   |    1101       1100
%                                    |

% Second parameter is 0 => Modulation without Gray labeling
tx_16qam_symbols_nogray = bits_to_16qam(signal_bits.', 0);

disp('First 10 Transmitted Symbols: With Gray Labeling');
disp(tx_16qam_symbols_gray(1:10).');

disp('First 10 Transmitted Symbols: Without Gray Labeling');
disp(tx_16qam_symbols_nogray(1:10).');

% ############### ADDING WHITE GAUSSIAN NOISE ###############

% Defining a SNR
snr = 10;

% To get the received noisy signal for gray labeling
noisy_signal_gray = add_noise(tx_16qam_symbols_gray, snr);

% To get the received noisy signal for without gray labeling
noisy_signal_nogray = add_noise(tx_16qam_symbols_nogray, snr);

figure(1)
scatter(real(noisy_signal_gray(1:500)), imag(noisy_signal_gray(1:500)))
title('Received Noisy Symbols (First 500) (With Gray Labelling)');
xlabel('Re');
ylabel('Im');
grid on;

figure(2)
scatter(real(noisy_signal_nogray(1:500)), imag(noisy_signal_nogray(1:500)))
title('Received Noisy Symbols (First 500) (Without Gray Labelling)');
xlabel('Re');
ylabel('Im');
grid on;

% ############### DEMODULATION ###############

% Calling the demodulation function for both cases
demod_signal_gray = demodulation(noisy_signal_gray);
demod_signal_nogray = demodulation(noisy_signal_nogray);

% Second parameter is 1 => with Gray labeling
rx_signal_bits_gray = qam16_to_bits(demod_signal_gray, 1);

% Second parameter is 0 => without Gray labeling
rx_signal_bits_nogray = qam16_to_bits(demod_signal_nogray, 0);

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
symbol_error_gray = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        signal_bits = randi([0,1],signal_length,1);

        % Generating QPSK symbols with gray labeling
        tx_16qam_symbols_gray = bits_to_16qam(signal_bits.', 1);
        
        % Adding noise and demodulating for each SNR
        noisy_signal_gray = add_noise(tx_16qam_symbols_gray, snr_array(i));
        
        % Demodulation
        demod_signal_gray = demodulation(noisy_signal_gray);

        % To get the binary bits corresponding to demodulated symbols
        rx_signal_bits_gray = qam16_to_bits(demod_signal_gray, 1);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = rx_signal_bits_gray - signal_bits.';
        symbolwise_diff = demod_signal_gray - tx_16qam_symbols_gray;
        
        symbol_error_gray(i) = symbol_error_gray(i) + calculate_ber(symbolwise_diff);

        bit_error_gray(i) = bit_error_gray(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error_gray = bit_error_gray/(signal_length*iterations);
symbol_error_gray = symbol_error_gray/((signal_length/4)*iterations);

% Plotting the BER vs. SNR graph
figure(3);
semilogy(snr_array, bit_error_gray, '-bo')
title('BER vs. SNR Plot (With Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% BER CALCULATION WITHOUT GRAY LABELING

% In this third section we calculate the BER without gray labeling
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
symbol_error_nogray = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        signal_bits = randi([0,1],signal_length,1);

        % Generating QPSK symbols with gray labeling
        tx_16qam_symbols_nogray = bits_to_16qam(signal_bits.', 0);
        
        % Adding noise and demodulating for each SNR
        noisy_signal_nogray = add_noise(tx_16qam_symbols_nogray, snr_array(i));
        
        % Demodulation
        demod_signal_nogray = demodulation(noisy_signal_nogray);

        % To get the binary bits corresponding to demodulated symbols
        rx_signal_bits_nogray = qam16_to_bits(demod_signal_nogray, 0);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = rx_signal_bits_nogray - signal_bits.';
        symbolwise_diff = demod_signal_nogray - tx_16qam_symbols_nogray;
        
        symbol_error_nogray(i) = symbol_error_nogray(i) + calculate_ber(symbolwise_diff);

        bit_error_nogray(i) = bit_error_nogray(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error_nogray = bit_error_nogray/(signal_length*iterations);
symbol_error_nogray = symbol_error_nogray/((signal_length/4)*iterations);

% Plotting the BER vs. SNR graph
figure(4);
semilogy(snr_array, bit_error_nogray, '-bo')
title('BER vs. SNR Plot (Without Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% THEORETICAL GRAPH USING Q-FUNCTION

% In this fourth section we first compute a theoretical curves defining the
% error probability variation with SNR for a 16QAM communication system.
% We then plot the theoretical curves with our practical curves.

% We define signal to noise ratio as:
% SNR = Eb/No

% In our case of 16QAM communication system, we have
% Eb = 5/2
% Also, the textbook defines
% sigma = sqrt(No/2)        (sigma => standard deviation of noise)

% So, we get the relation between sigma and SNR as:
% sigma = sqrt(5/(4*SNR))

% As the calculation of bit error rate is difficult for 16QAM
% We calculate the symbol error rate as follows:
% P_s = 3*Q(sqrt(4*SNR/5)) - (9/4)*[Q(sqrt(4*SNR/5))]^2

% Symbol Error Rate (SER) is the same for both of the above cases
symbol_error_theo = zeros(1, numel(snr_array));
for i = 1:numel(snr_array)
    qvalue = qfunc(sqrt(4*(10^(snr_array(i)/10))/5));
    symbol_error_theo(i) = 3*qvalue - (9/4)*qvalue*qvalue;
end

% ########## WITH GRAY LABELING ##########

% We get an approximation for the BER in this case:
% P_e = (3/4)*Q(sqrt(4*SNR/5))

ber_theo_gray = zeros(1, numel(snr_array));

% Calculation of theoretical error probability using Q-function
for i = 1:numel(snr_array)
    qvalue = qfunc(sqrt(4*(10^(snr_array(i)/10))/5));
    ber_theo_gray(i) = (3/4)*qvalue;
end

% ########## WITHOUT GRAY LABELING ##########

% We were unable to derive an expression for BER in this case.

figure(5);
semilogy(snr_array, ber_theo_gray, '-rx')
hold on
semilogy(snr_array, bit_error_gray, '-bo')
hold off
title('BER vs. SNR Plot (With Gray Labeling)');
xlabel('SNR (dBW)');
ylabel('BER');
legend('Theoretical BER (Gray Labeling)', 'Practical BER (Gray Labeling)');
legend
grid on;

figure(6);
semilogy(snr_array, symbol_error_theo, '-rx')
hold on
semilogy(snr_array, symbol_error_gray, '-bo')
semilogy(snr_array, symbol_error_nogray, '-go')
hold off
title('Symbol Error Rate vs. SNR Plot');
xlabel('SNR (dBW)');
ylabel('Symbol Error Rate');
legend('Theoretical Symbol Error Rate', 'Practical SER (Gray Labeling)', 'Practical SER (No Gray Labeling)');
legend
grid on;

figure(7);
semilogy(snr_array, ber_theo_gray, '-rx')
hold on
semilogy(snr_array, bit_error_gray, '-bo')
semilogy(snr_array, bit_error_nogray, '-go')
hold off
title('BER vs. SNR Plot (With/Without Gray Labeling, Theoretical and Practical)');
xlabel('SNR (dBW)');
ylabel('BER');
legend('Theoretical (Gray Labeling)', 'Practical BER (Gray Labeling)', 'Practical BER (No Gray Labeling)');
legend
grid on;

%% FUNCTIONS

% Function for 16QAM modulation with/without gray labeling
function tx_16qam_symbols = bits_to_16qam(signal_bits, gray_label)
    
    tx_16qam_symbols = zeros(1, numel(signal_bits)/4);
    symbols = [3+3i, 1+3i, -1+3i, -3+3i;
               3+1i, 1+1i, -1+1i, -3+1i;
               3-1i, 1-1i, -1-1i, -3-1i;
               3-3i, 1-3i, -1-3i, -3-3i];
           
    bit_combos_gray = ["0 0 0 0", "0 0 0 1", "0 0 1 1", "0 0 1 0";
                       "0 1 0 0", "0 1 0 1", "0 1 1 1", "0 1 1 0";
                       "1 1 0 0", "1 1 0 1", "1 1 1 1", "1 1 1 0";
                       "1 0 0 0", "1 0 0 1", "1 0 1 1", "1 0 1 0"];
                   
    bit_combos_nogray = ["0 0 0 0", "0 0 0 1", "0 0 1 0", "0 0 1 1";
                         "0 1 0 0", "0 1 0 1", "0 1 1 0", "0 1 1 1";
                         "1 0 0 0", "1 0 0 1", "1 0 1 0", "1 0 1 1";
                         "1 1 0 0", "1 1 0 1", "1 1 1 0", "1 1 1 1"];

    if gray_label == 1 % With Gray Labeling
        for i = 1:4:numel(signal_bits)-3
            % Finding the symbols corresponding to 4 bits at a time
            bits = strjoin(string(signal_bits(i:i+3)));
            index_mat = contains(bit_combos_gray, bits);
            tx_16qam_symbols(ceil(i/4)) = symbols(index_mat);
        end
    elseif gray_label == 0 % Without Gray Labeling
        for i = 1:4:numel(signal_bits)-3
            % Finding the symbols corresponding to 4 bits at a time
            bits = strjoin(string(signal_bits(i:i+3)));
            index_mat = contains(bit_combos_nogray, bits);
            tx_16qam_symbols(ceil(i/4)) = symbols(index_mat);
        end
    else
        disp('Invalid Input for Gray Labeling');
    end
end

% Function to add complex white gaussian noise to a signal
function noisy_signal = add_noise(tx_signal, snr)
    % Generating a random complex number with normal distribution
    Re_noise = sqrt(5/(4*(10^(snr/10))))*randn(1,numel(tx_signal));
    Im_noise = sqrt(5/(4*(10^(snr/10))))*randn(1,numel(tx_signal));
    noise = Re_noise + 1i*Im_noise;
    
    noisy_signal = tx_signal + noise;
end

% Function to perform demodulation using decision boundaries
function demod_signal = demodulation(noisy_signal)
    
    demod_signal = zeros(1, numel(noisy_signal));
    symbols = [3+3i, 1+3i, -1+3i, -3+3i;
               3+1i, 1+1i, -1+1i, -3+1i;
               3-1i, 1-1i, -1-1i, -3-1i;
               3-3i, 1-3i, -1-3i, -3-3i];
    
    for i = 1:numel(noisy_signal)
        % To get distance of signal point from each of the symbols
        distance = abs(symbols - noisy_signal(i));
        
        % To get the index of the minimum distance
        [~, index] = min(distance, [], 'all', 'linear');
        
        % The closest (distance wise) symbol is chosen
        demod_signal(i) = symbols(index);
    end
end

function signal_bits = qam16_to_bits(signal, gray_label)
    
    signal_bits = zeros(1, numel(signal)*4);
    
    symbols = [3+3i, 1+3i, -1+3i, -3+3i;
               3+1i, 1+1i, -1+1i, -3+1i;
               3-1i, 1-1i, -1-1i, -3-1i;
               3-3i, 1-3i, -1-3i, -3-3i];
           
    bit_combos_gray = ["0 0 0 0", "0 0 0 1", "0 0 1 1", "0 0 1 0";
                       "0 1 0 0", "0 1 0 1", "0 1 1 1", "0 1 1 0";
                       "1 1 0 0", "1 1 0 1", "1 1 1 1", "1 1 1 0";
                       "1 0 0 0", "1 0 0 1", "1 0 1 1", "1 0 1 0"];
                   
    bit_combos_nogray = ["0 0 0 0", "0 0 0 1", "0 0 1 0", "0 0 1 1";
                         "0 1 0 0", "0 1 0 1", "0 1 1 0", "0 1 1 1";
                         "1 0 0 0", "1 0 0 1", "1 0 1 0", "1 0 1 1";
                         "1 1 0 0", "1 1 0 1", "1 1 1 0", "1 1 1 1"];

    if gray_label == 1 % With gray labeling
        for i = 1:4:4*numel(signal)-3
            % Finding corresponding bit combinations
            index_mat = (symbols == signal(ceil(i/4)));
            signal_bits(i:i+3) = str2num(bit_combos_gray(index_mat));
        end
    elseif gray_label == 0 % Without gray labeling
        for i = 1:4:4*numel(signal)-3
            % Finding corresponding bit combinations
            index_mat = (symbols == signal(ceil(i/4)));
            signal_bits(i:i+3) = str2num(bit_combos_nogray(index_mat));
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
