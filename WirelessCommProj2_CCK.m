% CCK transmitter only implementation according to 802.11b methods
clear;close all;clc
numIter = 1; % The number of iterations of the simulation
nSym = 11000/2; % Let's say this is a simulation of 11 Mbps running for 0.01 s, the number of symbols is half of number of bits
bit2phase = [0, pi, pi/2, -pi/2]; % Phase array for CCK
chan = 1; % Only applying AWGN channel
M = 4;
bits = randi([0 1], 1, nSym*log2(M)); % Generate random bits
bitLen = length(bits);
chips = zeros(1, bitLen); % cck chips
dibits = zeros(1,4); % storage array for current iteration of 8 bits, converted to decimal
phi = zeros(1,4); % storage array for current iteration of phases
% Code the signal using CCK
for cckindex = 1:8:(bitLen-7) % Convert 8 bits to four phases. Each two bits can represent one phase
    for dibitindex = 1:4 % Find the decimal number associated with each two bits, convert to phase
        dibits(dibitindex) = bi2de(bits(cckindex+2*(dibitindex-1):cckindex+1+2*(dibitindex-1)), 'left-msb');
        phi(dibitindex) = bit2phase(dibits(dibitindex)+1); % Decimal is converted to phase
    end
    chips(cckindex:cckindex+7)=[exp(1j*sum(phi)),exp(1j*(phi(1)+phi(3)+phi(4))),...
        exp(1j*(phi(1)+phi(2)+phi(4))),-exp(1j*(phi(1)+phi(4))),exp(1j*(phi(1)+phi(2)+phi(3))),...
        exp(1j*(phi(1)+phi(3))),-exp(1j*(phi(1)+phi(2))),exp(1j*phi(1))]; % Use the cck formula to turn phases into 8 chips
end