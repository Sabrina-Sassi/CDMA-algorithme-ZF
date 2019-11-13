% Gold sequence generator
% Parameter descriptions
% r = number of binary shift register stage = 5
% N = period of m-sequence = 31
% g(D) = 1 + D^2 + D^5 = m-sequence generator function
% S0 = initial state of shift register 
% G = r x r m-sequence generator matrix
% b = output of Galois shift-register generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gold_waveform]=gold(r)
r = 5; % Number of stages
N = 2^r-1; % m-sequence period
% Memory allocation
S0 = zeros(1,r);
G_1 = zeros(r,r); 
G_2 = zeros(r,r);
G_temp_1 = zeros(r,r); % Temporary generator function for matrix multiplication
G_temp_2 = zeros(r,r); % Temporary generator function for matrix multiplication
b_1 = zeros(1,N);
b_2 = zeros(1,N);
% Initial state of shift register
S0 = [ 1 0 0 0 0 ]';
% Generator matrix of g(D)
G_1 = [ 0 1 0 0 1 ; 1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ;  ]';
G_2 = [ 0 1 1 1 1 ; 1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ;  ]';
G_temp_1 = G_1;
G_temp_2 = G_2;
% Find the N-1 outputs of the generator
b_1(1) = 1; % The first element of the output which is [ 1 0 0 0 0 ] x [ 1 0 0 0 0 ]' = 1
b_2(1) = 1; % The first element of the output which is [ 1 0 0 0 0 ] x [ 1 0 0 0 0 ]' = 1
for i = 2 : N,
    b_1(i) = mod([ 1 0 0 0 0 ] * G_temp_1 * S0 , 2);
    b_2(i) = mod([ 1 0 0 0 0 ] * G_temp_2 * S0 , 2);
    G_temp_1 = mod(G_temp_1 * G_1,2);
    G_temp_2 = mod(G_temp_2 * G_2,2);
end
% Generating the family of gold codes
gold = zeros(N , N+2); % There are N+2 number of gold codes with length N
gold_waveform = zeros(N , N+2); % Waveform consists of 1's and -1's
gold(:,1) = b_1';
gold(:,2) = b_2';
boy = length(gold(:,2)); % Find the length of the gold(:,2) function for circular shifting
tempGold2 = gold(:,2)';
% For the generation of the other gold codes, shift the second m-sequence
% and add to the first m-sequence
for i = 3 : N+2,
    gold2Shifted = [tempGold2(boy) tempGold2(1:(boy-1))]';
    gold(:,i) = mod(gold(:,1) + gold2Shifted, 2);
    tempGold2 = gold2Shifted';   
end
% Convert the spreading code to spreading waveform: 0 -> -1, 1 -> 1
gold_waveform = (gold.*2)-1;

% % Find the X-correlation values of the generated waveforms for proof
% correlationMat = zeros(33,33);
% for i = 1 : 33,
%     for j = 1 : 33,
%         correlationMat(i,j)=(gold_waveform(:,i)'*gold_waveform(:,j))./31;
%     end
% end
% 
%     


    
    