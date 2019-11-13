



clc;
clear vars;
format long;
ZF_MMSE=0;                     % Choix de l'algorithme (0 pour ZF/ 1 pour MMSE)
N = 100000;                    % Nombre de bits transmis
K = 30;                        % Nombre d'utilisateurs
g=5;                           % ordre du registre à décalage (fonction Gold)
C=2^g-1;                       % Facteur de gain
a=gold(g);                     % Matrice des codes (séquences de Gold)
S = a(:,(2:(K + 1)))./sqrt(C); % C x K = Matrice des codes normalisée 
Eb = 1; 
SNR_in_dB =0:2:20;             % Rapport signal à bruit
R = corrcoef(S);               % Coeffs de corralation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des coeffs de correlation
RR = zeros( K , K ); 
for i = 1 : K
for j = 1 : K
RR( i , j ) = S(:,i)' * S(:,j);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAI_in_dB = 0;                 % Puissance relative des interférants (en dB)
MAI = 10 ^ (MAI_in_dB/10);
A = sqrt(MAI) .* eye(K,K);     % K x K = matrice des amplitudes des utilisateurs
d = randsrc(K,N,[-1 1]);       % K x N = matrice des données des utilisateurs
d_estimate = zeros(K,N);       % K x N = matrice des données estimées des utilisateurs
y = zeros(K,N);                % K x N = matrice des données estimées aprés filter banks
t = zeros(C,N);                % C x N = matrice des données transmises
r = zeros(C,N);                % C x N = matrice des données reçues

for j=1:length(SNR_in_dB) 
SNR = 10 ^ (SNR_in_dB(j)/10);
sgma = sqrt(1/(2*SNR));        % Puissance du bruit (2*sgma^2= 1/SNR)
n = sgma * randn(C,N);         % Génération du bruit AWGN
L_mmse = inv(RR^1 + ZF_MMSE *(sgma ^ 2) .* (A ^ -2)); 
for i = 1 : N
t(:,i) = S * A * d(:,i);       % signal transmis
r(:,i) = t(:,i) + n(:,i);      % signal reçu
end
%%%%%%%%%%%%%%%%%%%%
% Filter banks
for i = 1 : N
for ii = 1 : K
y (ii,i) = S(:,ii)' * r(:,i); 
end
end
%%%%%%%%%%%%%%%%%%%%%
%estimation
for i=1 : N
d_estimate(:,i) = sign(L_mmse * y(:,i));
end
%%%%%%%%%%%%%%%%%%%%%
%Calcul du BER
error_number = length(find(d(2,:)-d_estimate(2,:))); % calcul du nombre de bits erronés
BER(j)=error_number/N;
end
semilogy (SNR_in_dB,BER,'sr-'); % Plot du BER Vs SNR
hold on;
title('Détection  CDMA');
xlabel('Rapport signal sur bruit (SNR)');
ylabel('BER');
legend('mf','zf','Mmse');
grid on;


CDMA_ZF_MMSE.m
Affichage de CDMA_ZF_MMSE.m en cours...