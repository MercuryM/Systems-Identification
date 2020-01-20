%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%
N = 2000;
batch_num = 100;
a = 1/2;
lamda_sqr = 9;
mu_sqr = 1;
y_initial = 0;
y_initial_2 = 0;
w_initial = 0;
y = zeros(N,1);
delay_y = zeros(N,1);
delay_y_AR2 = zeros(2,1,N); 
temp_mat_data_output = zeros(2,1,N); 
temp_mat_data_data = zeros(2,2,N);
A = zeros(N,1);
A(:) = 1;
B = [];
para_AR1 = zeros(batch_num,1); 
para_AR2 = zeros(2,1,batch_num);
for k = 1: 1: batch_num
e = sqrt(lamda_sqr) * randn(N,1);
w = sqrt(mu_sqr) * randn(N,1);
%%%%%%%%%%%%%%%%%%% Generate output and delay output %%%%%%%%%%%%%%%%%%%
for i = 1: 1: N
if i == 1
y(i) = a * y_initial + e(i) + w(i) - a * w_initial; delay_y(i) = y_initial;
delay_y_AR2(:,:,i) = [y_initial; y_initial_2];
elseif i == 2
y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); delay_y(i) = y(i - 1);
delay_y_AR2(:,:,i) = [y(i - 1); y_initial];
else
y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); delay_y(i) = y(i - 1);
delay_y_AR2(:,:,i) = [y(i - 1); y(i - 2)];
end
temp_mat_data_output(:,:,i) = delay_y_AR2(:,:,i) * y(i);
temp_mat_data_data(:,:,i) = delay_y_AR2(:,:,i) * delay_y_AR2(:,:,i)';
end
B = [B y];
%%%%%%%%%%%%%%%%%%% Least-square Algorithm %%%%%%%%%%%%%%%%%%%
% AR(1)
para_AR1(k) = sum(delay_y .* y)./ (sum(delay_y .* delay_y));
% AR(2)
para_AR2(:,:,k) = (sum(temp_mat_data_data, 3)) \ sum(temp_mat_data_output, 3);
end
%%%%%%%%%%%%%%%%%%% Obtain the parameters %%%%%%%%%%%%%%%%%%%
para_AR1_average = mean(para_AR1);
para_AR2_average = mean(para_AR2,3);
%%%%%%%%%%%%%%%%%%% the variance of the parameters %%%%%%%%%%%%%%%%%%%
para_AR1_var = var(para_AR1);
temp1 = mean((para_AR2(1,:,:) - para_AR2_average(1,:)).^ 2); 
temp2 = mean((para_AR2(1,:,:) - para_AR2_average(1,:)) .*(para_AR2(2,:,:) - para_AR2_average(2,:)));
temp3 = mean((para_AR2(2,:,:) - para_AR2_average(2,:)).^ 2); 
para_AR2_var = [temp1, temp2; temp2, temp3];
%%%%%%%%%%%%%%%%%%% the variance of prediction errors %%%%%%%%%%%%%%%%%%%
prediction_error_AR1 = mean(y - para_AR1(batch_num).* delay_y); 
prediction_error_AR1_var = var((y - para_AR1(batch_num).* delay_y));
delay_y_t_2 = [0; delay_y(1:N-1)];
prediction_error_AR2 = mean(y - para_AR2(1,:,batch_num) .* delay_y - para_AR2(2,:,batch_num) .* delay_y_t_2);
prediction_error_AR2_var = var((y - para_AR2(1,:,batch_num) .* delay_y - para_AR2(2,:,batch_num) .* delay_y_t_2));
epsilon = y - para_AR1(batch_num).* delay_y; cov_fun = [];
for i = 0: 1: N-1
temp = (1 / N) * (sum(epsilon(1: N-i) .* epsilon(1+i: N)));
    cov_fun = [cov_fun; temp];
end
adtest(cov_fun)




%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%
% n = 100000;
% para1 = zeros(100,1);
% para2 = zeros(2,1,100);
% m = 1;
% for N = 1000: 1000: n
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
% para_AR1 = zeros(batch_num,1); 
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
% para_AR1(k) = sum(delay_y .* y)./ (sum(delay_y .* delay_y));
% % AR(2)
% % para_AR2(:,:,k) = (sum(temp_mat_data_data, 3)) \ sum(temp_mat_data_output, 3);
% end
% %%%%%%%%%%%%%%%%%%% Obtain the parameters %%%%%%%%%%%%%%%%%%%
% para_AR1_average = mean(para_AR1);
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
% para1(m) = para_AR1_average;
% para2(:,:,m) = para_AR2_average;
% m = m + 1
% end
% para1_ref = zeros(100,1);
% para2_1_ref = zeros(100,1);
% para2_2_ref = zeros(100,1);
% para1_ref(:) = 0.4615;
% para2_1_ref(:) = 0.4511;
% para2_2_ref(:) = 0.0226;
% x = 1000:1000:100000;
% figure (1)
% plot(x,para1,'LineWidth',1)
% hold on
% plot(x,para1_ref,'LineWidth',2)
% grid on;
% legend('Simulated \theta_1^o (batch num = 100)', 'Analyzed \theta_1^o')
% xlabel('\itN')
% ylabel('\theta_1^o')
% title('Simulated \theta_1^o versus \itN')
% figure(2)
% plot(x,para2(1,:))
% hold on
% plot(x,para2_1_ref,'LineWidth',2)
% hold on
% plot(x,para2(2,:))
% hold on
% plot(x,para2_2_ref,'LineWidth',2)
% grid on;
% legend('Simulated \theta_2^o(1)','Analyzed \theta_2^o(1)','Simulated \theta_2^o(2)','Analyzed \theta_2^o(2)')
% xlabel('\itN')
% ylabel('\theta_2^o')
% title('Simulated \theta_2^o versus \itN')

% n = 100000;
% para1_1 = zeros(100,1);
% para2 = zeros(2,1,100);
% m = 1;
% for N = 1000: 1000: n
% batch_num = 1;
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
% para2(:,:,m) = para_AR2_average;
% m = m + 1
% end
% para1_ref = zeros(100,1);
% para2_1_ref = zeros(100,1);
% para2_2_ref = zeros(100,1);
% para1_ref(:) = 0.4615;
% para2_1_ref(:) = 0.4511;
% para2_2_ref(:) = 0.0226;
% x = 1000:1000:100000;
% plot(x,para1_1,'LineWidth',1)
% hold on 
% plot(x,para1_ref,'LineWidth',2)
% % plot(x,para1_ref,'LineWidth',2)
% grid on;
% legend('Simulated \theta_1^o (batch num = 100)','Simulated \theta_1^o (batch num = 1)', 'Analyzed \theta_1^o')
% xlabel('\itN')
% ylabel('\theta_1^o')
% title('Simulated \theta_1^o versus \itN')

