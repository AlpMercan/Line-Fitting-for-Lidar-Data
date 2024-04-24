clc;
clear;
close all;

load('HW2data.mat'); 

x_pos = rangedata_Q2 .* cos(angledata_Q2 * pi / 180);
y_pos = rangedata_Q2 .* sin(angledata_Q2 * pi / 180);

N = length(rangedata_Q2);
points = [x_pos', y_pos'];% to get the points on x-y plane
%% clustering 
% Calculate differences in radial distances and angles, to group the large
% data points
radial_diff = [0; diff(rangedata_Q2(:))]; 
angle_diff = [0; diff(angledata_Q2(:))]; 


radial_threshold = prctile(abs(radial_diff), 93); % 93 based on trial and error
angle_threshold = prctile(abs(angle_diff), 93); 


segments = ones(N, 1);
for i = 2:N
    if abs(radial_diff(i)) > radial_threshold || abs(angle_diff(i)) > angle_threshold
        segments(i) = segments(i-1) + 1; 
    else
        segments(i) = segments(i-1); %
    end
end
segment_points = cell(max(segments), 1);


for k = 1:max(segments)

    segment_indices = find(segments == k);
    
    % Extract the [x, y] points for the current segment
    segment_points{k} = points(segment_indices, :);
end


figure;
hold on;
colors = hsv(max(segments)); 
for k = 1:max(segments)
    segment_points_k = segment_points{k};
    plot(segment_points_k(:, 1), segment_points_k(:, 2), '.', 'MarkerEdgeColor', colors(k, :));
end
%% Line with dp
max_distance_threshold = 0.12; % dp


segments_to_check = [1; N]; 


finished_segments = []; % the finished s


while ~isempty(segments_to_check)

    segment_indices = segments_to_check(:, 1);
    segments_to_check(:, 1) = [];
    start_idx = segment_indices(1);
    end_idx = segment_indices(2);
    start_point = points(start_idx, :);
    end_point = points(end_idx, :);
    
    segment_point_indices = start_idx:end_idx;
    segment_points_subset = points(segment_point_indices, :);

    distances = arrayfun(@(idx) pointLineDistance(segment_points_subset(idx - start_idx + 1, :), start_point, end_point), segment_point_indices);

    [max_distance, max_idx_rel] = max(distances);
    max_idx = segment_point_indices(max_idx_rel); 
    
    if max_distance > max_distance_threshold
        segments_to_check = [segments_to_check, [start_idx; max_idx], [max_idx; end_idx]];
    else
        finished_segments = [finished_segments, [start_idx; end_idx]];
    end
end

figure;
hold on;
colors = hsv(max(segments));
for k = 1:max(segments)
    scatter(points(segments == k, 1), points(segments == k, 2), [], colors(k, :), 'filled');
end

for i = 1:size(finished_segments, 2)
    start_idx = finished_segments(1, i);
    end_idx = finished_segments(2, i);
    plot(points([start_idx, end_idx], 1), points([start_idx, end_idx], 2), 'k-', 'LineWidth', 2);
end

title('Segmented Fitting Line to Clustered Data');
xlabel('X Position');
ylabel('Y Position');
hold off;

num_segments = max(segments);

segment_radials = cell(num_segments, 1);
segment_angles = cell(num_segments, 1);

for k = 1:num_segments

    cart_points = segment_points{k}; 
    

    radials = sqrt(sum(cart_points.^2, 2));  
    

    angles = atan2d(cart_points(:, 2), cart_points(:, 1)); 
    segment_radials{k} = radials;
    segment_angles{k} = angles;
end

for k = 1:num_segments
    fprintf('Segment %d:\n', k);
    fprintf('Radial Distances:\n');
    disp(segment_radials{k});
    [alpha, R] = computeLineParameters(segment_angles{k}, segment_radials{k});
    fprintf('alpha:\n');
    disp(alpha);
    fprintf('R:\n');
    disp(R);
    fprintf('Angles (degrees):\n');
    disp(segment_angles{k});
end

figure;
hold on;
colors = hsv(max(segments));
for k = 1:max(segments)
    segment_points_k = segment_points{k};
    plot(segment_points_k(:, 1), segment_points_k(:, 2), '.', 'MarkerEdgeColor', colors(k, :));


    [alpha, R] = computeLineParameters(atan2d(segment_points_k(:,2), segment_points_k(:,1)), sqrt(sum(segment_points_k.^2, 2)));
    str = sprintf('Segment %d: \\alpha=%.2f, R=%.2f', k, alpha, R);
    text(mean(segment_points_k(:, 1)), mean(segment_points_k(:, 2)), str, 'FontSize', 8, 'BackgroundColor', 'white');
end

title('alfa and R values');
xlabel('X Position');
ylabel('Y Position');
hold off;

function d = pointLineDistance(point, lineStart, lineEnd)
    A = lineEnd(2) - lineStart(2);
    B = lineStart(1) - lineEnd(1);
    C = -A * lineStart(1) - B * lineStart(2);
    d = abs(A * point(1) + B * point(2) + C) / sqrt(A^2 + B^2);
end
function [alfa, R] = computeLineParameters(angledata, rangedata)

    Q = angledata ;

    weights = 1 ./ rangedata.^2;

    N = numel(rangedata);

    tot_num1 = sum(weights .* rangedata.^2 .* sin(2 * Q));
    total = 0;
    for i = 1:N
        for k = 1:N
            num2 = weights(i) * weights(k) * rangedata(i) * rangedata(k) * cos(Q(i)) * sin(Q(k));
            total = total + num2;
        end
    end
    tot_weights = sum(weights);
    final_result = tot_num1 - (2 * total / tot_weights);

    detot_num2 = sum(weights .* rangedata.^2 .* cos(2 * Q));
    detot_num3 = 0;
    for i = 1:N
        for k = 1:N
            num3 = weights(i) * weights(k) * rangedata(i) * rangedata(k) * cos(Q(i) + Q(k));
            detot_num3 = detot_num3 + num3;
        end
    end
    De_final_result = detot_num2 - (detot_num3 / tot_weights);

    alfa = 0.5 * atan2(final_result, De_final_result) + pi / 2;

    R = sum(weights .* rangedata .* cos(Q - alfa)) / sum(weights);

end
