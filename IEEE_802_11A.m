%%%%% IEEE 802.11a %%%%%
%Designed by Haochen Yin
%It may take several hours to run the whole program in 10k loops
%Part 3 and Part 4 may be crushed, better to Comment one when run other

clear all
close all
clc

%%%%% Parameters %%%%%

N = 64;
Ncp = 16;
Ns = N+Ncp;
Symbols = 12;
M = 4; 
fs = 20000000;
Ts = 1/fs;
delta_f = fs/N;
looptimes = 10000;
% SNR = 0:2:20; % Use in Part 1-3
SNR = [0, 10, 20]; %Only use in Part 4
h = [1, 0.9, 0.5];
H_err = zeros(1,64);
% alphas = 1; % Use in Part 1-3
alphas = 0:0.1:1; %Only use in Part 4
R_1 = zeros(11,3);
R_2 = zeros(11,3);

%%%%% DATA %%%%%

msg = randsrc(N, Symbols, (0:M-1));
x_data = qammod(msg, M)/(sqrt(2));
x_data_t = ifft(x_data, 64, 1);
x_data_t_cp = [x_data_t(N - Ncp + 1:end,:);x_data_t];
x = reshape(x_data_t_cp, 1, Ns*Symbols); 

%$%%% STF %%%%%

j= sqrt(-1);
S_26_26 = sqrt(1/2)* [0,0,1+j, 0,0,0,-1-j,0,0,0,1+j,0,0,0,-1-j,0,0,0,-1-j,0,0,0,1+j,0,0,0,0,0,0,0,-1-j,0,0,0,-1-j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0];      
x_stf = zeros(1,160);
for n = 0:159
    temp_sum_1 = 0;
    for k = -26:26
        temp_sum_1 = temp_sum_1 + S_26_26(k + 27) * exp(j * 2 * pi * k * delta_f * n * Ts);
    end
    x_stf(n+1) = sqrt(1/12) * temp_sum_1;
end

%$%%% LTF %%%%%

L_26_26 = [0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,0,0,0,0,0,0];
x_ltf = ifft(L_26_26, 64, 2);
x_ltf = [x_ltf, x_ltf, x_ltf(1:32)];

%%%% Simulate %%%%%
x_n = [x_stf, x_ltf, x];
z_stf = x_stf(1:16);

n = 0:15;
figure;
plot(n,abs(z_stf));
xlabel('n');
ylabel('|z_S_T_F(n)|');
title('|z_S_T_F(n)| versus n');

figure;
plot(n,angle(z_stf));
xlabel('n');
ylabel('\angle z_S_T_F(n)');
title('\angle z_S_T_F(n) versus n');
   
H = fft(h, 64, 2);
y = conv(x_n, h);
count = 0;
for i = SNR 
    det = 0;
    count = count + 1;
    H_error_1 = 0;
    H_error_2 = 0;
    sym_error_1 = 0;
    sym_error_2 = 0;
    for alpha = alphas
        for loop = 1 : looptimes
            clear peak;
%             peak = 0;
            y_n1 = awgn(y(1:160), i, 'measured');
            y_n2 = awgn(y(161:320), i, 'measured');
            y_n3 = awgn(y(321:1282), i, 'measured');
            y_n = [y_n1, y_n2, y_n3];
            
            %%%%% Part 3: Channel estimate %%%%%
            
            v = y_n(225:288);
            L = fft(v, 64, 2);
            for k = [10, 22]
                H_est(k) = (L(k) / L_26_26(k));
                H_err(k) = (abs(H_est(k) - H(k))^2)/(abs(H(k))^2);         
            end
  
            H_error_1 = H_error_1 + H_err(10);
            H_error_2 = H_error_2 + H_err(22);
            
            %%%%% Part 2: Packet detection %%%%%
                    
            r = zeros(1,160);
            for n = 1:160
                for m = 0:15
                    r(n) = r(n) + (conj(y_n(n + m)) * z_stf(m + 1));
                end
            end
            r = abs(r);
            peak = findpeaks(r,'minpeakheight',mean(r) + 2 * std(r),'minpeakdistance',15);
            if length(peak) > 8
                det = det + 1;
            end
            
            %%%%% Part 3: BER %%%%%
        
            y_data_t_cp = reshape(y_n(321:1280),Ns,Symbols);
            y_data_t = y_data_t_cp(Ncp+1:end,:);
            y_data = fft(y_data_t, 64, 1);
             
            for symbol = 1:Symbols
                for sb = [10,22]
                    data_r(sb,symbol) = y_data(sb,symbol) / (H_est(sb));
                    msg_r(sb,symbol) = qamdemod(data_r(sb,symbol), M);
                end
                if (msg_r(10,symbol) ~= msg(10,symbol))
                    sym_error_1 = sym_error_1 + 1;
                end
                if (msg_r(22,symbol) ~= msg(22,symbol))
                    sym_error_2 = sym_error_2 + 1;
                end
            end
            
            %%%% Part 4: power allocation %%%%%

                    
            eta_1 = log2(1 + (abs(H_est(10)))^2 * alpha * i);
            eta_2 = log2(1 + (abs(H_est(22)))^2 * (1 - alpha) * i);
            M_1_1 = floor(2.^eta_1);
            M_1_2 = floor(2.^eta_2);

            if ((M_1_1 < 2) || (M_1_2 < 2))
                M_1_1 = 1;
                M_1_2 = 1;
            end
            if ((M_1_1 > 1) && (M_1_2 > 1))
                M_1_1 = 2;
                M_1_2 = 2;
            end
            if ((M_1_1 > 3) && (M_1_2 > 3))
                M_1_1 = 4;
                M_1_2 = 4;
            end
            if ((M_1_1 > 15) && (M_1_2 > 15))
                M_1_1 = 16;
                M_1_2 = 16;
            end
            if ((M_1_1 > 63) && (M_1_2 > 63))
                M_1_1 = 64;
                M_1_2 = 64;
            end

            R_1(alpha * 10 + 1, count) = R_1(alpha * 10 + 1, count) + 3/4 * 2 * log2(M_1_1) * delta_f / (10^6);
            
            M_2_1 = floor(2.^eta_1);
            M_2_2 = floor(2.^eta_2);

            if ((M_2_1 > 1) && (M_2_1 < 4))
                M_2_1 = 2;
            end
            if ((M_2_1 > 3) && (M_2_1 < 16))
                M_2_1 = 4;
            end
            if ((M_2_1 > 15) && (M_2_1 < 64))
                M_2_1 = 16;
            end
            if (M_2_1 > 63) 
                M_2_1 = 64;
            end
            if ((M_2_2 > 1) && (M_2_2 < 4))
                M_2_2 = 2;
            end
            if ((M_2_2 > 3) && (M_2_2 < 16))
                M_2_2 = 4;
            end
            if ((M_2_2 > 15) && (M_2_2 < 64))
                M_2_2 = 16;
            end
            if (M_2_2 > 63) 
                M_2_2 = 64;
            end
            
            R_2(alpha * 10 + 1, count) = R_2(alpha * 10 + 1, count) + 3/4 * 2 * (log2(M_2_1) + log2(M_2_2)) * delta_f / (10^6);
        end
    
    %%%%% Part 2 reslut %%%%%
    dets(count) = det / length(alphas) / looptimes;
    %%%%% Part 3.1 reslut %%%%%
    H_error_set_1(count) = H_error_1 / length(alphas) / looptimes;
    H_error_set_2(count) = H_error_2 / length(alphas) / looptimes;
    %%%%% Part 3.3 reslut %%%%%
    sym_error_set_1(count) = sym_error_1 / length(alphas) / (Symbols * looptimes);
    sym_error_set_2(count) = sym_error_2 / length(alphas) / (Symbols * looptimes);
    end
end
    
figure;
plot(SNR, dets);
xlabel('SNR/dB');
ylabel('probability of packet detectio');
title('probability of packet detectio versus SNR');  

figure;
plot(SNR, H_error_set_1, SNR, H_error_set_2);
xlabel('SNR/dB');
ylabel('mean of normalized squared');
title('estimation error versus SNR'); 

figure;
semilogy(SNR, sym_error_set_1, SNR, sym_error_set_2);
grid on;
xlabel('SNR/dB');
ylabel('BER');
title('BER versus SNR');  

figure
plot(H)
xlabel('k');
ylabel('H');
title('Frequency response of the Channel');  


for i = SNR
    figure
    plot(alphas, R_1(:,i/10+1) / looptimes);
    xlabel('\alpha');
    ylabel('R/Mbps');
    title(['mean sum rate for SNR = ', num2str(i),'dB']);
end



for i = SNR
    figure
    plot(alphas, R_2(:,i/10+1) / looptimes);
    xlabel('\alpha');
    ylabel('R/Mbps');
    title(['mean sum rate for SNR = ', num2str(i),'dB']);
end