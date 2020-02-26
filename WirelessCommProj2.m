% Wireless Comm Project 2
% Simulation of IEEE 802.11b 2Mbps Transmitter & Reciever for 1 ms
clear;close all;clc
nSym = 1000; % 2 Mbps is 1 Msymbols/s, in a period of 1 ms 1000 symbols are sent
numIter = 2000; % The simulation is run 1000 times to generate a better BER curve
SNR_Vec = 0:2:16; % Vector containing values of SNR to iterate through
lenSNR = length(SNR_Vec);
M = 4; % We are using QPSK for 2 Mbps
ber = zeros(1, lenSNR); % Vector to hold mean bit error rates
berVec = zeros(numIter, lenSNR); % Vector for holding ber at each SNR for all iterations
barkcode = transpose([1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1]); % The Barker Code
% Run the simulation numIter amount of times
for i = 1:numIter
bits = randi([0 1], 1, nSym*log2(M)); % Generate random bits, twice the number of symbols
bitLen = length(bits);
% Convert pairs of binary to decimal before using QPSK
msg = zeros(1, bitLen/log2(M));
msgindex = 1;
for bitindex = 1:log2(M):(bitLen)
msg(msgindex) = bi2de(bits(bitindex:bitindex+log2(M)-1), 'left-msb');
msgindex = msgindex+1;
end
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
    tx = qammod(msg,M); % QPSK modulate the signal
    txCoded = tx.*barkcode; % Multiply by Barker code
    txNoisy = awgn(txCoded,SNR_Vec(j)+ 10*log10(log2(M)/11)); % AWGN Channel, scaled to fit chip rate of Barker Code
    rxDecoded = mean(txNoisy.*barkcode,1); % Multiply again by Barker Code and find the mean
    rxDemod = qamdemod(rxDecoded,M); % Demodulate
    % Converting the symbols back to bits
    rxMSG = zeros(1, bitLen);
    rxMSGindex = 1;
    for rxDemindex = 1:1:length(rxDemod)
    rxMSG(rxMSGindex:rxMSGindex+log2(M)-1) = de2bi(rxDemod(rxDemindex), log2(M), 2, 'left-msb');
    rxMSGindex = rxMSGindex + log2(M);
    end
    % Compute and store the BER for this iteration
    [~, berVec(i,j)] = biterr(bits(1:end), rxMSG(1:end)); % We're interested in the BER, which is the 2nd output of BITERR
    end % End SNR iteration
    
end % End numIter iteration
% Compute and plot the mean BER
ber(1,:) = mean(berVec,1);
figure;
semilogy(SNR_Vec, ber(1,:));
title('BER of 802.11b through AWGN Channel');
xlabel('SNR');
ylabel('BER');
% Compute the theoretical BER for an AWGN channel
M = 4;
berTheory = berawgn(SNR_Vec + 10*log10(log2(M)),'psk',M,'nondiff');
figure;
semilogy(SNR_Vec,berTheory,'r')
title('Theoretical BER for AWGN Channel');
xlabel('SNR');
ylabel('BER');