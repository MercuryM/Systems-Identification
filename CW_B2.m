N = 2000;
batch_num = 100;
a = 1/2;
lamda_sqr = 9;
mu_sqr = 1;
y_initial = 0;
y_initial_2 = 0;
w_initial = 0;
e = zeros(N,1);
y = zeros(N,1);
delay_y = zeros(N,1);
delay_y_AR2 = zeros(2,1,N); 
temp_mat_data_output = zeros(2,1,N); 
temp_mat_data_data = zeros(2,2,N);
B = [];
para_AR1 = zeros(batch_num,1); 
para_AR2 = zeros(2,1,batch_num);
for k = 1: 1: batch_num
    for i = 1: 1: N
        if i == 1
            e(i) = randn(1);
        else
            e(i) = -1/2 * e(i-1) + randn(1);
        end
    end
    w = sqrt(mu_sqr) * randn(N,1);
% %%%%%%%%%%%%%%%%%%% Generate output and delay output %%%%%%%%%%%%%%%%%%%
for i = 1: 1: N
if i == 1
y(i) = a * y_initial + e(i) + w(i) - a * w_initial; 
delay_y(i) = y_initial;
delay_y_AR2(:,:,i) = [y_initial; y_initial_2];
elseif i == 2
y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); 
delay_y(i) = y(i - 1);
delay_y_AR2(:,:,i) = [y(i - 1); y_initial];
else
y(i) = a * y(i - 1) + e(i) + w(i) - a * w(i - 1); 
delay_y(i) = y(i - 1);
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
temp2 = mean((para_AR2(1,:,:) - para_AR2_average(1,:)) .* (para_AR2(2,:,:) - para_AR2_average(2,:)));
temp3 = mean((para_AR2(2,:,:) - para_AR2_average(2,:)).^ 2); 
para_AR2_var = [temp1, temp2; temp2, temp3];
%%%%%%%%%%%%%%%%%%% the variance of prediction errors %%%%%%%%%%%%%%%%%%%
prediction_error_AR1 = mean(y - para_AR1(batch_num).* delay_y); prediction_error_AR1_var = var((y - para_AR1(batch_num).* delay_y));
delay_y_t_2 = [0;delay_y(1:N-1)];
prediction_error_AR2 = mean(y - para_AR2(1,:,batch_num) .* delay_y - para_AR2(2,:,batch_num) .* delay_y_t_2); prediction_error_AR2_var = var((y - para_AR2(1,:,batch_num) .* delay_y -para_AR2(2,:,batch_num) .* delay_y_t_2));
epsilon = y - para_AR1(batch_num).* delay_y; 
cov_fun = [];
for i = 0: 1: N-1
temp = (1 / N) * (sum(epsilon(1: N-i) .* epsilon(1+i: N)));
    cov_fun = [cov_fun; temp];
end
adtest(cov_fun, 'Alpha', 0.01)