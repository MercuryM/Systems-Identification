%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%
N = 1000;
batch_num = 100;
a = 1/2;
lamda_sqr = 9;
mu_sqr = 1;
y_initial = 0;
y_initial_2 = 0;
w_initial = 0;
S_initial = 0;
alpha = 1;
alpha_1 = 0.001;
alpha_2 = 100;
para_AR1_rec = zeros(batch_num,1); 
para_AR2_rec = zeros(2,1,batch_num); 
para_AR1_R_rec = zeros(batch_num,1); 
para_AR2_R_rec = zeros(2,1,batch_num); 
para_AR1_V_rec = zeros(batch_num,1); 
para_AR1_V_1_rec = zeros(batch_num,1);
para_AR1_V_2_rec = zeros(batch_num,1);
para_AR2_V_rec = zeros(2,1,batch_num);
para_AR1_B = zeros(batch_num,1); 
para_AR2_B = zeros(2,1,batch_num);
for k = 1: 1: batch_num
    y = zeros(N,1);
    delay_y = zeros(N,1);
    delay_y_AR2 = zeros(2,1,N); 
    temp_mat_data_output = zeros(2,1,N); 
    temp_mat_data_data = zeros(2,2,N);
    %1st
    S_AR1 = zeros(N,1);
    error_AR1 = zeros(N,1);
    K_AR1 = zeros(N,1);
    para_AR1 = zeros(N,1);
    %2nd
    R_AR1 = zeros(N,1);
    error_AR1_R = zeros(N,1);
    K_AR1_R = zeros(N,1);
    para_AR1_R = zeros(N,1);
    %3rd
    %alpha
    V_AR1 = zeros(N,1);
    V_initial = alpha * 1;
    error_AR1_V = zeros(N,1);
    K_AR1_V = zeros(N,1);
    para_AR1_V = zeros(N,1);
    para_AR1_V_initial = 0;
    %alpha_1
    V_AR1_1 = zeros(N,1);
    V_initial_1 = alpha_1 * 1;
    error_AR1_V_1 = zeros(N,1);
    K_AR1_V_1 = zeros(N,1);
    para_AR1_V_1 = zeros(N,1);
    para_AR1_V_initial_1 = 0;
    %alpha_2
    V_AR1_2 = zeros(N,1);
    V_initial_2 = alpha_2 * 1;
    error_AR1_V_2 = zeros(N,1);
    K_AR1_V_2 = zeros(N,1);
    para_AR1_V_2 = zeros(N,1);
    para_AR1_V_initial_2 = 0;
    %1st
    S_AR2 = zeros(2,2,N);
    error_AR2 = zeros(N,1);
    K_AR2 = zeros(2,1,N);
    para_AR2 = zeros(2,1,N);
    %2nd
    R_AR2 = zeros(2,2,N);
    error_AR2_R = zeros(N,1);
    K_AR2_R = zeros(2,1,N);
    para_AR2_R = zeros(2,1,N);
    %3rd
    V_AR2 = zeros(2,2,N);
    V_initial_AR2 = alpha * eye(2,2); 
    error_AR2_V = zeros(N,1);
    K_AR2_V = zeros(2,1,N);
    para_AR2_V = zeros(2,1,N); 
    para_AR2_V_intial = [0;0];
    % for i = 1: 1: N
        % ifi == 1
            % e(i) = randn(1);
        % else
            % e(i) = -1/2 * e(i-1) + randn(1); 
        % end
    % end
    e = sqrt(lamda_sqr) * randn(N,1);
    w = sqrt(mu_sqr) * randn(N,1);
%%%%%%%%%%%%%%%%%%% Generate output and delay output %%%%%%%%%%%%%%%%%%%
    for i = 1: 1: N
        if i == 1
            y(i) = a * y_initial + e(i) + w(i) - a * w_initial; 
            delay_y(i) = y_initial;
            delay_y_AR2(:,:,i) = [y_initial; y_initial_2];
            %3rd_1
            temp = 1 + delay_y(i)' * V_initial * delay_y(i); 
            V_AR1(i) = V_initial - (1 / temp) * V_initial * delay_y(i) * delay_y(i)' * V_initial; 
            error_AR1_V(i) = y(i) - delay_y(i)' * para_AR1_V_initial;
            K_AR1_V(i) = V_AR1(i) * delay_y(i);
            para_AR1_V(i) = para_AR1_V_initial + K_AR1_V(i) * error_AR1_V(i);
            %3rd_2
            temp_1 = 1 + delay_y(i)' * V_initial_1 * delay_y(i); V_AR1_1(i) = V_initial_1 - (1 / temp_1) * V_initial_1 * delay_y(i) * delay_y(i)' * V_initial_1; 
            error_AR1_V_1(i) = y(i) - delay_y(i)' * para_AR1_V_initial_1;
            K_AR1_V_1(i) = V_AR1_1(i) * delay_y(i);
            para_AR1_V_1(i) = para_AR1_V_initial_1 + K_AR1_V_1(i) * error_AR1_V_1(i);
            %3rd_3
            temp_2 = 1 + delay_y(i)' * V_initial_2 * delay_y(i); V_AR1_2(i) = V_initial_2 - (1 / temp_2) * V_initial_2 * delay_y(i) * delay_y(i)' * V_initial_2; 
            error_AR1_V_2(i) = y(i) - delay_y(i)' * para_AR1_V_initial_2;
            K_AR1_V_2(i) = V_AR1_2(i) * delay_y(i);
            para_AR1_V_2(i) = para_AR1_V_initial_2 + K_AR1_V_2(i) * error_AR1_V_2(i);
            %3rd
            temp_AR2 = 1 + delay_y_AR2(:, :, i)' * V_initial_AR2 * delay_y_AR2(:, :, i);
            V_AR2(:, :, i) = V_initial_AR2 - (1 / temp_AR2) * V_initial_AR2 * delay_y_AR2(:, :, i) * delay_y_AR2(:, :, i)' * V_initial_AR2;
            error_AR2_V(i) = y(i) - delay_y_AR2(:, :, i)' * para_AR2_V_intial; K_AR2_V(:, :, i) = V_AR2(:, :, i) * delay_y_AR2(:, :, i);
            para_AR2_V(:, :, i) = para_AR2_V_intial + K_AR2_V(:, :, i) * error_AR2_V(i);
        elseif i == 2
            y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); 
            delay_y(i) = y(i-1);
            delay_y_AR2(:,:,i) = [y(i-1); y_initial];
            %1st
            S_AR1(i) = S_AR1(i-1) + delay_y(i) * delay_y(i)'; 
            error_AR1(i) = y(i) - delay_y(i)' * para_AR1(i-1);
            K_AR1(i) = S_AR1(i) \ delay_y(i);
            para_AR1(i) = para_AR1(i-1) + K_AR1(i) * error_AR1(i);
            %2nd
            R_AR1(i) = R_AR1(i-1) + 1/i *(delay_y(i) * delay_y(i)' - R_AR1(i-1)); 
            error_AR1_R(i) = y(i) - delay_y(i)' * para_AR1_R(i-1);
            K_AR1_R(i) = 1/i * (R_AR1(i) \ delay_y(i));
            para_AR1_R(i) = para_AR1_R(i-1) + K_AR1_R(i) * error_AR1_R(i);
            %3rd_1
            temp = 1 + delay_y(i)' * V_AR1(i-1) * delay_y(i); 
            V_AR1(i) = V_AR1(i-1) - (1 / temp) * V_AR1(i-1) * delay_y(i) * delay_y(i)' * V_AR1(i-1); 
            error_AR1_V(i) = y(i) - delay_y(i)' * para_AR1_V(i-1);
            K_AR1_V(i) = V_AR1(i) * delay_y(i);
            para_AR1_V(i) = para_AR1_V(i-1) + K_AR1_V(i) * error_AR1_V(i);
            %3rd_2
            temp_1 = 1 + delay_y(i)' * V_AR1_1(i-1) * delay_y(i); 
            V_AR1_1(i) = V_AR1_1(i-1) - (1 / temp_1) * V_AR1_1(i-1) * delay_y(i) * delay_y(i)' * V_AR1_1(i-1);
            error_AR1_V_1(i) = y(i) - delay_y(i)' * para_AR1_V_1(i-1);
            K_AR1_V_1(i) = V_AR1_1(i) * delay_y(i);
            para_AR1_V_1(i) = para_AR1_V_1(i-1) + K_AR1_V_1(i) * error_AR1_V_1(i);
            %3rd_3
            temp_2 = 1 + delay_y(i)' * V_AR1_2(i-1) * delay_y(i); 
            V_AR1_2(i) = V_AR1_2(i-1) - (1 / temp_2) * V_AR1_2(i-1) * delay_y(i) * delay_y(i)' * V_AR1_2(i-1); 
            error_AR1_V_2(i) = y(i) - delay_y(i)' * para_AR1_V_2(i-1);
            K_AR1_V_2(i) = V_AR1_2(i) * delay_y(i);
            para_AR1_V_2(i) = para_AR1_V_2(i-1) + K_AR1_V_2(i) * error_AR1_V_2(i);
            %3rd
            temp_AR2 = 1 + delay_y_AR2(:, :, i)' * V_AR2(:, :, i-1) * delay_y_AR2(:, :, i);
            V_AR2(:, :, i) = V_AR2(:, :, i-1) - (1 / temp_AR2) * V_AR2(:, :, i-1) * delay_y_AR2(:, :, i) * delay_y_AR2(:, :, i)' * V_AR2(:, :, i-1);
            error_AR2_V(i) = y(i) - delay_y_AR2(:, :, i)' * para_AR2_V(:, :, i-1); 
            K_AR2_V(:, :, i) = V_AR2(:, :, i) * delay_y_AR2(:, :, i);
            para_AR2_V(:, :, i) = para_AR2_V(:, :, i-1) + K_AR2_V(:, :, i) * error_AR2_V(i);
        else
            y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); 
            delay_y(i) = y(i-1);
            delay_y_AR2(:,:,i) = [y(i-1); y(i-2)];
            %1st
            S_AR1(i) = S_AR1(i-1) + delay_y(i) * delay_y(i)';
            error_AR1(i) = y(i) - delay_y(i)' * para_AR1(i-1);
            K_AR1(i) = delay_y(i) / S_AR1(i);
            para_AR1(i) = para_AR1(i-1) + K_AR1(i) * error_AR1(i);
            %2nd
            R_AR1(i) = R_AR1(i-1) + 1/i *(delay_y(i) * delay_y(i)' - R_AR1(i-1)); 
            error_AR1_R(i) = y(i) - delay_y(i)' * para_AR1_R(i-1);
            K_AR1_R(i) = 1/i * (R_AR1(i) \ delay_y(i));
            para_AR1_R(i) = para_AR1_R(i-1) + K_AR1_R(i) * error_AR1_R(i);
            %3rd_1
            temp = 1 + delay_y(i)' * V_AR1(i-1) * delay_y(i); 
            V_AR1(i) = V_AR1(i-1) - (1 / temp) * V_AR1(i-1) * delay_y(i) * delay_y(i)' * V_AR1(i-1); error_AR1_V(i) = y(i) - delay_y(i)' * para_AR1_V(i-1);
            K_AR1_V(i) = V_AR1(i) * delay_y(i);
            para_AR1_V(i) = para_AR1_V(i-1) + K_AR1_V(i) * error_AR1_V(i);
            %3rd_2
            temp_1 = 1 + delay_y(i)' * V_AR1_1(i-1) * delay_y(i); 
            V_AR1_1(i) = V_AR1_1(i-1) - (1 / temp_1) * V_AR1_1(i-1) * delay_y(i) * delay_y(i)' * V_AR1_1(i-1); 
            error_AR1_V_1(i) = y(i) - delay_y(i)' * para_AR1_V_1(i-1);
            K_AR1_V_1(i) = V_AR1_1(i) * delay_y(i);
            para_AR1_V_1(i) = para_AR1_V_1(i-1) + K_AR1_V_1(i) * error_AR1_V_1(i);
            %3rd_3
            temp_2 = 1 + delay_y(i)' * V_AR1_2(i-1) * delay_y(i); 
            V_AR1_2(i) = V_AR1_2(i-1) - (1 / temp_2) * V_AR1_2(i-1) * delay_y(i) * delay_y(i)' * V_AR1_2(i-1); 
            error_AR1_V_2(i) = y(i) - delay_y(i)' * para_AR1_V_2(i-1);
            K_AR1_V_2(i) = V_AR1_2(i) * delay_y(i);
            para_AR1_V_2(i) = para_AR1_V_2(i-1) + K_AR1_V_2(i) * error_AR1_V_2(i);
            %1st
            S_AR2(:,:,i) = S_AR2(:,:,i-1) + delay_y_AR2(:,:,i) * delay_y_AR2(:,:,i)';
            error_AR2(i) = y(i) - delay_y_AR2(:,:,i)' * para_AR2(:,:,i-1);
            K_AR2(:,:,i) = S_AR2(:,:,i) \ delay_y_AR2(:,:,i);
            para_AR2(:,:,i) = para_AR2(:,:,i-1) + K_AR2(:,:,i) .* error_AR2(i);
            %2nd
            R_AR2(:,:,i) = R_AR2(:,:,i-1) + 1/i * (delay_y_AR2(:,:,i) * delay_y_AR2(:,:,i)' - R_AR2(:,:,i-1));
            error_AR2_R(i) = y(i) - delay_y_AR2(:,:,i)' * para_AR2_R(:,:,i-1);
            K_AR2_R(:,:,i) = 1/i * (R_AR2(:,:,i) \ delay_y_AR2(:,:,i)); 
            para_AR2_R(:,:,i) = para_AR2_R(:,:,i-1) + K_AR2_R(:,:,i) .* error_AR2_R(i);
            %3rd
            temp_AR2 = 1 + delay_y_AR2(:, :, i)' * V_AR2(:, :, i-1) * delay_y_AR2(:, :, i);
            V_AR2(:, :, i) = V_AR2(:, :, i-1) - (1 / temp_AR2) * V_AR2(:, :, i-1) * delay_y_AR2(:, :, i) * delay_y_AR2(:, :, i)' * V_AR2(:, :, i-1);
            error_AR2_V(i) = y(i) - delay_y_AR2(:, :, i)' * para_AR2_V(:, :, i-1);
            K_AR2_V(:, :, i) = V_AR2(:, :, i) * delay_y_AR2(:, :, i);
            para_AR2_V(:, :, i) = para_AR2_V(:, :, i-1) + K_AR2_V(:, :, i) * error_AR2_V(i);
        end
        temp_mat_data_output(:,:,i) = delay_y_AR2(:,:,i) * y(i);
        temp_mat_data_data(:,:,i) = delay_y_AR2(:,:,i) * delay_y_AR2(:,:,i)';
    end
    % AR(1)
    para_AR1_rec(k) = para_AR1(N);
    para_AR1_R_rec(k) = para_AR1_R(N);
    para_AR1_V_rec(k) = para_AR1_V(N);
    para_AR1_V_1_rec(k) = para_AR1_V_1(N);
    para_AR1_V_2_rec(k) = para_AR1_V_2(N);
    para_AR1_B(k) = (sum(delay_y .* delay_y)) \ sum(delay_y .* y);
    % AR(2)
    para_AR2_rec(:, :, k) = para_AR2(:,:,N);
    para_AR2_R_rec(:, :, k) = para_AR2_R(:,:,N);
    para_AR2_V_rec(:, :, k) = para_AR2_V(:,:,N);
    para_AR2_B(:, :, k) = (sum(temp_mat_data_data, 3)) \ sum(temp_mat_data_output, 3);
end
%%%%%%%%%%%%%%%%%%% Obtain the parameters %%%%%%%%%%%%%%%%%%%
para_AR1_rec_average = mean(para_AR1_rec); 
para_AR1_R_rec_average = mean(para_AR1_R_rec); 
para_AR1_V_rec_average = mean(para_AR1_V_rec); 
para_AR1_V_1_rec_average = mean(para_AR1_V_1_rec); 
para_AR1_V_2_rec_average = mean(para_AR1_V_2_rec); 
para_AR1_B_average = mean(para_AR1_B);
R_mean = mean(delay_y .* delay_y);
y_var = var(y);

para_AR2_rec_average = mean(para_AR2_rec,3); 
para_AR2_R_rec_average = mean(para_AR2_R_rec,3); 
para_AR2_V_rec_average = mean(para_AR2_V_rec,3); 
para_AR2_B_average = mean(para_AR2_B,3);

%%%%%%%%%%%%%%%%%%% the variance of the parameters %%%%%%%%%%%%%%%%%%%
para_AR1_rec_var = var(para_AR1_rec); 
para_AR1_R_rec_var = var(para_AR1_R_rec);
para_AR1_V_rec_var = var(para_AR1_V_rec);
para_AR1_B_var = var(para_AR1_B);
para_AR1_reference_1 = zeros(N,1); 
para_AR1_reference_1(:) = 0.4615; 
para_AR2_temp_1 = reshape(para_AR2(1,1,:),[N,1]); 
para_AR2_reference_1 = zeros(N,1); 
para_AR2_reference_1(:) = 0.4511; 
para_AR2_temp_2 = reshape(para_AR2(2,1,:),[N,1]); 
para_AR2_reference_2 = zeros(N,1); 
para_AR2_reference_2(:) = 0.0226;
R_mreference = zeros (N,1);
R_mreference(:) = R_mean;
y_var_reference = zeros(N,1);
y_var_reference(:)=y_var;
figure(1)
x = 1:1:N;
plot(x,para_AR1,'o','Linewidth',0.25)
hold on
plot(x,para_AR1_R,'-','Linewidth',1)
hold on
plot(x,para_AR1_V,':','Linewidth',1)
hold on
plot(x,para_AR1_reference_1,'Linewidth',2);
legend('RLS 1st','RLS 2nd','RLS 3rd','Analyzed')
xlabel('\itN')
ylabel('\theta_1^o')
title('Simulated \theta_1^o versus \itN')
grid on
figure(2)
plot(x, para_AR2_temp_1,'Linewidth',1)
hold on
plot(x, para_AR2_reference_1,'Linewidth',2)
hold on
plot(x, para_AR2_temp_2,'Linewidth',1)
hold on
plot(x, para_AR2_reference_2,'Linewidth',2);
legend('Simulated \theta_2^o(1)','Analyzed \theta_2^o(1)','Simulated \theta_2^o(2)','Analyzed \theta_2^o(2)')
xlabel('\itN')
ylabel('\theta_2^o')
title('Simulated \theta_2^o versus \itN')
grid on
figure(3)
plot(x,para_AR1_V,x,para_AR1_V_1,x,para_AR1_V_2,'Linewidth',1)
hold on
plot(x,para_AR1_reference_1,'Linewidth',2); 
legend('alpha = 1','alpha = 0.001','alpha = 100');
xlabel('\itN')
ylabel('\theta_1^o')
title('Simulated \theta_1^o versus \itN \rmusing RLS 3rd')
grid on
figure(4)
plot(x,S_AR1,'Linewidth',2);
xlabel('\itN')
ylabel('\itS')
title('\itS \rmversus \itN \rm(RLS 1st)')
grid on
figure(5)
plot(x,R_AR1,x,y_var_reference,'Linewidth',2);
xlabel('\itN')
ylabel('\itR')
title('\itR \rmversus \itN \rm(RLS 2nd)')
grid on
figure(6)
plot(x,V_AR1,x,V_AR1_1,x,V_AR1_2,'Linewidth',2); 
legend('alpha = 1','alpha = 0.001','alpha = 100');
xlabel('\itN')
ylabel('\itV')
xlim([0,200])
ylim([0,0.01])
title('\itV \rmversus \itN \rm(RLS 3rd)')
grid on
% figure(7)
% plot(x,para_AR1_R,'Linewidth',1)
% hold on
% 
% 
% n = 1000;
% para1_1 = zeros(100,1);
% para2 = zeros(2,1,100);
% m = 1;
% for N = 1:1:n
% batch_num = 100;
% a = 1/2;
% lamda_sqr = 9;
% mu_sqr = 1;
% y_initial = 0;
% y_initial_2 = 0;
% w_initial = 0;
% y = zeros(N,1);
% delay_y = zeros(N,1);
% delay_y_AR2 = zeros(2,1,N); 
% temp_mat_data_output = zeros(2,1,N); 
% temp_mat_data_data = zeros(2,2,N);
% A = zeros(N,1);
% A(:) = 1;
% B = [];
% para_AR1_1 = zeros(batch_num,1); 
% para_AR2 = zeros(2,1,batch_num);
% for k = 1: 1: batch_num
% e = sqrt(lamda_sqr) * randn(N,1);
% w = sqrt(mu_sqr) * randn(N,1);
% %%%%%%%%%%%%%%%%%%% Generate output and delay output %%%%%%%%%%%%%%%%%%%
% for i = 1: 1: N
% if i == 1
% y(i) = a * y_initial + e(i) + w(i) - a * w_initial; 
% delay_y(i) = y_initial;
% % delay_y_AR2(:,:,i) = [y_initial; y_initial_2];
% elseif i == 2
% y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); 
% delay_y(i) = y(i - 1);
% % delay_y_AR2(:,:,i) = [y(i - 1); y_initial];
% else
% y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); 
% delay_y(i) = y(i - 1);
% % delay_y_AR2(:,:,i) = [y(i - 1); y(i - 2)];
% end
% % temp_mat_data_output(:,:,i) = delay_y_AR2(:,:,i) * y(i);
% % temp_mat_data_data(:,:,i) = delay_y_AR2(:,:,i) * delay_y_AR2(:,:,i)';
% end
% B = [B y];
% %%%%%%%%%%%%%%%%%%% Least-square Algorithm %%%%%%%%%%%%%%%%%%%
% % AR(1)
% para_AR1_1(k) = sum(delay_y .* y)./ (sum(delay_y .* delay_y));
% % AR(2)
% % para_AR2(:,:,k) = (sum(temp_mat_data_data, 3)) \ sum(temp_mat_data_output, 3);
% end
% %%%%%%%%%%%%%%%%%%% Obtain the parameters %%%%%%%%%%%%%%%%%%%
% para_AR1_1_average = mean(para_AR1_1);
% % para_AR2_average = mean(para_AR2,3);
% %%%%%%%%%%%%%%%%%%% the variance of the parameters %%%%%%%%%%%%%%%%%%%
% % para_AR1_var = var(para_AR1);
% % temp1 = mean((para_AR2(1,:,:) - para_AR2_average(1,:)).^ 2); 
% % temp2 = mean((para_AR2(1,:,:) - para_AR2_average(1,:)) .*(para_AR2(2,:,:) - para_AR2_average(2,:)));
% % temp3 = mean((para_AR2(2,:,:) - para_AR2_average(2,:)).^ 2); 
% % para_AR2_var = [temp1, temp2; temp2, temp3];
% %%%%%%%%%%%%%%%%%%% the variance of prediction errors %%%%%%%%%%%%%%%%%%%
% % prediction_error_AR1 = mean(y - para_AR1(batch_num).* delay_y); 
% % prediction_error_AR1_var = var((y - para_AR1(batch_num).* delay_y));
% % delay_y_t_2 = [0; delay_y(1:N-1)];
% % prediction_error_AR2 = mean(y - para_AR2(1,:,batch_num) .* delay_y - para_AR2(2,:,batch_num) .* delay_y_t_2);
% % prediction_error_AR2_var = var((y - para_AR2(1,:,batch_num) .* delay_y - para_AR2(2,:,batch_num) .* delay_y_t_2));
% % epsilon = y - para_AR1(batch_num).* delay_y; cov_fun = [];
% % for i = 0: 1: N-1
% % temp = (1 / N) * (sum(epsilon(1: N-i) .* epsilon(1+i: N)));
% %     cov_fun = [cov_fun; temp];
% % end
% % adtest(cov_fun)
% para1_1(m) = para_AR1_1_average;
% % para2(:,:,m) = para_AR2_average;
% m = m + 1
% end
% para1_ref = zeros(1000,1);
% para2_1_ref = zeros(1000,1);
% para2_2_ref = zeros(1000,1);
% para1_ref(:) = 0.4615;
% para2_1_ref(:) = 0.4511;
% para2_2_ref(:) = 0.0226;
% x = 1:1:n;
% plot(x,para1_1,'LineWidth',1)
% hold on 
% plot(x,para1_ref,'LineWidth',2)
% grid on;
% legend('Simulated \theta_1^o (RLS 2nd)','Simulated \theta_1^o (batch num = 100)', 'Analyzed \theta_1^o')
% xlabel('\itN')
% ylabel('\theta_1^o')
% title('Simulated \theta_1^o versus \itN')