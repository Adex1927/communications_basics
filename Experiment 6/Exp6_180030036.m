%% DBPSK & DQPSK MODULATION AND DEMODULATION
% This is the experiment number 6 of Communications Lab course.

% In this first section an example of DBPSK & DQPSK modulation and 
% demodulation is demonstrated with randomly generated bits.

% ############### GENERATING MESSAGE ###############

% We choose a length of our random signal to be transmitted
signal_length = 10000; % Chosen to be a multiple of 2

% A random signal of binary digits (0 and 1) is generated
signal_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(signal_bits(1:20));

% ############### BPSK & QPSK DIGITAL MODULATION ###############

% We take 1 bit from the message signal at a time for BPSK,
% and 2 bits at a time for QPSK (we use mapping with Gray labelling).

% FOR BPSK mapping is done as follows:
% 0 -> -1
% 1 -> 1
% Generating BPSK symbols from the bits
bpsk_symbols = bpsk_modulation(signal_bits);

% FOR QPSK mapping is done as follows:
% 00 -> (1 + i)/sqrt(2)
% 01 -> (-1 + i)/sqrt(2)
% 11 -> (-1 - i)/sqrt(2)
% 10 -> (1 - i)/sqrt(2)
% Generating QPSK symbols from the bits
% Second parameter is 1 => Modulation with Gray labeling
qpsk_symbols = qpsk_modulation(signal_bits.', 1);

% ############### DIFFERENTIAL ENCODING ###############

% We map every possible value of phase difference between consecutive
% symbols to some value and transmit these new sequence of symbols.

% FOR BPSK
% 0 degree -> 1
% 180/-180 degrees -> -1
% We assume the first encoded symbol to be -1
dbpsk_symbols = dbpsk_encoding(bpsk_symbols);

% FOR QPSK
% 0 degree -> (1 + i)/sqrt(2)
% 90/-270 degrees -> (-1 + i)/sqrt(2)
% -90/270 degrees -> (1 - i)/sqrt(2)
% 180/-180 degrees -> (-1 - i)/sqrt(2)
% We assume the first encoded symbol to be (1 + i)/sqrt(2)
dqpsk_symbols = dqpsk_encoding(qpsk_symbols);

% ############### ADDING WHITE GAUSSIAN NOISE ###############

% Defining a SNR
snr = 5;

% To get the received noisy signal for DBPSK
dbpsk_noisy = add_noise_bpsk(dbpsk_symbols.', snr);

% To get the received noisy signal for DQPSK
dqpsk_noisy = add_noise_qpsk(dqpsk_symbols.', snr);

% ############### DEMODULATION ###############

% Using demodulation for BPSK
dbpsk_rx_symbols = bpsk_demodulation(dbpsk_noisy, 1);

% Using demodulation for QPSK
dqpsk_rx_symbols = qpsk_demodulation(dqpsk_noisy);

% ############### DIFFERENTIAL DECODING ###############

% Decoding for DBPSK
bpsk_rx_symbols = dbpsk_decoding(dbpsk_rx_symbols);

% Decoding for DQPSK
qpsk_rx_symbols = dqpsk_decoding(dqpsk_rx_symbols);
qpsk_rx_symbols = qpsk_demodulation(qpsk_rx_symbols);

% Converting symbols to signal bits for both cases
bpsk_signal_bits = bpsk_symbols_to_bits(bpsk_rx_symbols);
qpsk_signal_bits = qpsk_symbols_to_bits(qpsk_rx_symbols, 1);

% ############### BIT ERROR RATE ###############

% BER in case of DBPSK
bitwise_difference = bpsk_signal_bits - signal_bits;
ber_dbpsk = calculate_ber(bitwise_difference)/signal_length;

% BER in case of DQPSK
bitwise_difference = qpsk_signal_bits - signal_bits.';
ber_dqpsk = calculate_ber(bitwise_difference)/signal_length;

disp('BER in case of DBPSK:');
disp(ber_dbpsk);

disp('BER in case of DQPSK:');
disp(ber_dqpsk);

%% BER CALCULATION FOR DBPSK

% In this second section we calculate the bit error rate for DBPSK
% We have an error resolution of 10^-7
% Please reduce the number of iterations to get lower run time.

% Defining the signal length
signal_length = 100000;

% Defining the number of iterations the signal is tested for each SNR
iterations = 100;

% We define SNR array in dBW
snr_array = 0:0.2:15;

%Vector to store bit error rate corresponding to each SNR value
bit_error_dbpsk = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        signal_bits = randi([0,1],signal_length,1);

        % Generating QPSK symbols with gray labeling
        bpsk_symbols = bpsk_modulation(signal_bits);
        
        % Encoding
        dbpsk_symbols = dbpsk_encoding(bpsk_symbols);
        
        % Adding noise and demodulating for each SNR
        dbpsk_noisy = add_noise_bpsk(dbpsk_symbols.', snr_array(i));
        
        % Demodulation
        dbpsk_rx_symbols = bpsk_demodulation(dbpsk_noisy, 1);
        
        % Decoding
        bpsk_rx_symbols = dbpsk_decoding(dbpsk_rx_symbols);

        % To get the binary bits corresponding to demodulated symbols
        bpsk_signal_bits = bpsk_symbols_to_bits(bpsk_rx_symbols);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = bpsk_signal_bits - signal_bits;

        bit_error_dbpsk(i) = bit_error_dbpsk(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error_dbpsk = bit_error_dbpsk/(signal_length*iterations);

% Plotting the BER vs. SNR graph
figure(1);
semilogy(snr_array, bit_error_dbpsk, '-bo')
title('BER vs. SNR Plot for DBPSK');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% BER CALCULATION FOR DQPSK

% In this second section we calculate the bit error rate for DQPSK
% We have an error resolution of 10^-7
% Please reduce the number of iterations to get lower run time.

% Defining the signal length
signal_length = 100000;

% Defining the number of iterations the signal is tested for each SNR
iterations = 100;

% We define SNR array in dBW
snr_array = 0:0.2:15;

%Vector to store bit error rate corresponding to each SNR value
bit_error_dqpsk = zeros(numel(snr_array),1);

for i = 1:numel(snr_array)
    
    for j = 1:iterations
        % Generating a random signal with binary bits (0 and 1)
        signal_bits = randi([0,1],signal_length,1);

        % Generating QPSK symbols with gray labeling
        qpsk_symbols = qpsk_modulation(signal_bits.', 1);
        
        % Encoding
        dqpsk_symbols = dqpsk_encoding(qpsk_symbols);
        
        % Adding noise and demodulating for each SNR
        dqpsk_noisy = add_noise_qpsk(dqpsk_symbols.', snr_array(i));
        
        % Demodulation
        dqpsk_rx_symbols = qpsk_demodulation(dqpsk_noisy);
        
        % Decoding
        qpsk_rx_symbols = dqpsk_decoding(dqpsk_rx_symbols);
        qpsk_rx_symbols = qpsk_demodulation(qpsk_rx_symbols);

        % To get the binary bits corresponding to demodulated symbols
        qpsk_signal_bits = qpsk_symbols_to_bits(qpsk_rx_symbols, 1);

        % To get a bitwise difference between the tx and rx signal
        bitwise_difference = qpsk_signal_bits - signal_bits.';

        bit_error_dqpsk(i) = bit_error_dqpsk(i) + calculate_ber(bitwise_difference);
    end
end

% To get the fractional error
bit_error_dqpsk = bit_error_dqpsk/(signal_length*iterations);

% Plotting the BER vs. SNR graph
figure(2);
semilogy(snr_array, bit_error_dqpsk, '-bo')
title('BER vs. SNR Plot for DQPSK');
xlabel('SNR (dBW)');
ylabel('BER');
grid on;

%% THEORETICAL GRAPH USING Q-FUNCTION

% We define signal to noise ratio as:
% SNR = Eb/No

% We get the BER for the case of DBPSK as follows:
% P_e = 2*Q(sqrt(2Eb/No))*(1 - Q(sqrt(2Eb/No)))

ber_theoretical_dbpsk = zeros(1, numel(snr_array));

% Calculation of theoretical error probability for DBPSK
for i = 1:numel(snr_array)
    qvalue = qfunc(sqrt(2*(10^(snr_array(i)/10))));
    ber_theoretical_dbpsk(i) = 2*qvalue*(1-qvalue);
end

% We are unable to get a theoretical expression for DQPSK

figure(3);
semilogy(snr_array, bit_error_dbpsk, '-bo')
hold on
semilogy(snr_array, ber_theoretical_dbpsk, '-rx')
hold off
title('BER vs. SNR Plot for DBPSK');
xlabel('SNR (dBW)');
ylabel('BER');
xlim([0 10]);
legend('Experimental curve', 'Theoretical curve using Q-function');
legend
grid on;

%% FUNCTIONS

function signal_symbols = bpsk_modulation(signal_bits)
    signal_symbols = signal_bits*2 - 1;
end

function tx_qpsk_symbols = qpsk_modulation(signal_bits, gray_label)
    
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

function dbpsk_symbols = dbpsk_encoding(bpsk_symbols)
    dbpsk_symbols = zeros(numel(bpsk_symbols)+1, 1);
    dbpsk_symbols(1) = -1; % First symbol is -1
    
    bpsk_phase = angle(bpsk_symbols)*(180/pi);
    
    for i = 2:numel(bpsk_symbols)+1
        dbpsk_phase = angle(dbpsk_symbols(i-1))*(180/pi);
        
        if abs(bpsk_phase(i-1) - dbpsk_phase) == 0
            dbpsk_symbols(i) = 1;
            
        elseif abs(bpsk_phase(i-1) - dbpsk_phase) == 180
            dbpsk_symbols(i) = -1;
            
        end
    end
end

function dqpsk_symbols = dqpsk_encoding(qpsk_symbols)
    dqpsk_symbols = zeros(numel(qpsk_symbols)+1, 1);
    dqpsk_symbols(1) = (1 + 1i)/sqrt(2); % First symbol is (1 + i)/sqrt(2)
    
    qpsk_phase = angle(qpsk_symbols)*(180/pi);
    
    for i = 2:numel(qpsk_symbols)+1
        dqpsk_phase = angle(dqpsk_symbols(i-1))*(180/pi);
        
        if abs(qpsk_phase(i-1) - dqpsk_phase) == 0
            dqpsk_symbols(i) = (1 + 1i)/sqrt(2);
            
        elseif qpsk_phase(i-1) - dqpsk_phase == 90
            dqpsk_symbols(i) = (-1 + 1i)/sqrt(2);
            
        elseif qpsk_phase(i-1) - dqpsk_phase == -270
            dqpsk_symbols(i) = (-1 + 1i)/sqrt(2);
            
        elseif qpsk_phase(i-1) - dqpsk_phase == -90
            dqpsk_symbols(i) = (1 - 1i)/sqrt(2);
            
        elseif qpsk_phase(i-1) - dqpsk_phase == 270
            dqpsk_symbols(i) = (1 - 1i)/sqrt(2);
            
        elseif abs(qpsk_phase(i-1) - dqpsk_phase) == 180
            dqpsk_symbols(i) = (-1 - 1i)/sqrt(2);
            
        end
    end
end

% Function to add noise to DQPSK symbols
function noisy_signal = add_noise_qpsk(tx_signal, snr)
    % Generating a random complex number with normal distribution
    Re_noise = sqrt(1/(4*(10^(snr/10))))*randn(1,numel(tx_signal));
    Im_noise = sqrt(1/(4*(10^(snr/10))))*randn(1,numel(tx_signal));
    noise = Re_noise + 1i*Im_noise;
    
    noisy_signal = tx_signal + noise;
end

% Function to add noise to DBPSK symbols
function noisy_signal = add_noise_bpsk(tx_signal, snr)
    % Generating a random number with normal distribution
    noise = sqrt(1/(2*(10^(snr/10))))*randn(1,numel(tx_signal));
    
    noisy_signal = tx_signal + noise;
end

function rx_symbols = bpsk_demodulation(rx_signal_noisy, samples)
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

function demod_signal = qpsk_demodulation(noisy_signal)
    
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

function qpsk_rx_symbols = dqpsk_decoding(dqpsk_rx_symbols)
    qpsk_rx_symbols = zeros(numel(dqpsk_rx_symbols)-1, 1);
    
    phases = [0 90 -90 180];
    symbols = [(1 + 1i) (-1 + 1i) (1 - 1i) (-1 - 1i)]/sqrt(2);
    
    for i = 2:numel(dqpsk_rx_symbols)-1
        previous_phase = angle(dqpsk_rx_symbols(i-1));
        
        if dqpsk_rx_symbols(i) == symbols(1)
            current_phase = phases(1)*(pi/180);
            
        elseif dqpsk_rx_symbols(i) == symbols(2)
            current_phase = phases(2)*(pi/180);
            
        elseif dqpsk_rx_symbols(i) == symbols(3)
            current_phase = phases(3)*(pi/180);
            
        elseif dqpsk_rx_symbols(i) == symbols(4)
            current_phase = phases(4)*(pi/180);            
        end
        
        qpsk_rx_symbols(i-1) = exp(1i*(previous_phase+current_phase));
    end
end

function bpsk_rx_symbols = dbpsk_decoding(dbpsk_rx_symbols)
    bpsk_rx_symbols = zeros(numel(dbpsk_rx_symbols)-1, 1);
    
    for i = 2:numel(dbpsk_rx_symbols)-1
        previous_phase = angle(dbpsk_rx_symbols(i-1));
        
        if dbpsk_rx_symbols(i) == 1
            current_phase = 0*(pi/180);
            
        elseif dbpsk_rx_symbols(i) == -1
            current_phase = 180*(pi/180);          
        end
        
        bpsk_rx_symbols(i-1) = real(exp(1i*(previous_phase+current_phase)));
    end
end

function signal_bits = bpsk_symbols_to_bits(signal_symbols)
    signal_bits = (signal_symbols + 1)/2;
end

function signal_bits = qpsk_symbols_to_bits(signal, gray_label)
    
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