% Load the .mat file attached to this assignment in Matlab. Then write your own code that fits a line
% to angledata_Q1 and rangedata_Q1 using the weighted least squares solution. Submit your code,
% and plots, and clearly label the alpha and r values found. The plots should include the given sensor
% data, the weights assigned as error bars, and the fitted line.
clc;
clear;
close all;
load('HW2data.mat');
p = rangedata_Q1;
Q = angledata_Q1 * (pi / 180);
weights = 1 ./ rangedata_Q1.^2;
N = size(rangedata_Q1, 2);
%% Numerator Part
tot_num1 = 0;
for i = 1:N
    num1 = weights(i) * p(i)*p(i) * sin(2 * Q(i));
    tot_num1 = tot_num1 + num1;
end

total = 0;  
for i = 1:N
    for k = 1:N
        num2 = weights(i) * weights(k) * p(i) * p(k) * cos(Q(i)) * sin(Q(k));
        total = total + num2;  
    end
end
tot_weights = sum(weights);
final_result = tot_num1 - (2*total /(tot_weights));

%% Denominator Part
detot_num2=0;
for i=1:N
    num3 = weights(i) * p(i)* p(i) * cos(2 * Q(i));
    detot_num2 = detot_num2 + num3;
end

detot_num3 = 0;  
for i = 1:N
    for k = 1:N
        num3 = weights(i) * weights(k) * p(i) * p(k) * cos(Q(i)+Q(k));
        detot_num3 = detot_num3 + num3;  
    end
end
tot_weights = sum(weights);
De_final_result = detot_num2 - (detot_num3 /(1* tot_weights));

%% alfa
inside=final_result/De_final_result;
alfa= 0.5*atan2(final_result,De_final_result)+pi/2;

%% R Calculation
R = sum(weights .* p .* cos(Q - alfa)) / sum(weights);
%% Plottings
range_errors = sqrt(1 ./ weights);  

angle_errors = (std(Q) / mean(Q)) * Q;

x = p .* cos(Q);
y = p .* sin(Q);

x_errors = sqrt((range_errors .* cos(Q)).^2 + (p .* angle_errors .* -sin(Q)).^2); 
y_errors = sqrt((range_errors .* sin(Q)).^2 + (p .* angle_errors .* cos(Q)).^2); 

error_scale = 0.1;  
x_errors_scaled = x_errors * error_scale;
y_errors_scaled = y_errors * error_scale;

figure;
errorbar(x, y, y_errors_scaled, x_errors_scaled, 'o');  
title('Sensor Data With Fitted Line Using Weighted Least Squares');
xlabel('x');
ylabel('y');
hold on;

start_x = R * cos(alfa);
start_y = R * sin(alfa);

coefficients = polyfit(x, y, 1);  
x_fit = linspace(min(x), max(x), 100);

y_fit = coefficients(1) * x_fit + coefficients(2);

plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

legend('Sensor Data', 'Linear Fit');

text(mean(x), max(y), sprintf('Alpha: %.4f radians', alfa), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(mean(x), max(y), sprintf('R: %.4f m', R), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

hold off;



