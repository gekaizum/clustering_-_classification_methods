%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q 1.0 Creating dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
m = 3; % length of vector (dimension)
c = 4; % number of clusters
colors = ['r', 'g', 'b','m','k']; % colors for clusters

% Init matrix and vectors
all_vectors = [];
labels = [];
cluster_centers = [];

for i = 1:c
    n = randi([10,100]);% number of points in cluster
    mean_val=randi([-20,20]);% mean value for cluster
    std_devs=2*rand+0.1;% standart deviation for cluster
    % Creating vector of cluster
    cluster_vectors = mean_val + std_devs * randn(n, m);
    % Adding cluster to space matrix
    all_vectors = [all_vectors; cluster_vectors];
    % Adding labels for clusters
    labels = [labels; repmat([i,n,mean_val,std_devs], n, 1)]; %cluster ID/number of points in cluster/mean value/standart deviation
    % Calculating center of cluster
    cluster_centers = [cluster_centers; mean(cluster_vectors, 1)];
end
% Sort labels by the third column (center coordinates)
%all_vectors = sortrows(labels, 3);

% Vizualization
figure;
hold on;
for i = 1:c
    cluster_points = all_vectors(labels(:,1) == i, :);
    scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(i), 'filled');
    % Grid line from center to points in cluster
    for j = 1:size(cluster_points, 1)
        plot3([cluster_points(j, 1), cluster_centers(i, 1)], ...
              [cluster_points(j, 2), cluster_centers(i, 2)], ...
              [cluster_points(j, 3), cluster_centers(i, 3)], colors(i));
    end
end
hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Clusters Visualization with Lines to Cluster Centers');
grid on; 
if m==3
    view(3); % 3D vizualization
end
xlim auto;
ylim auto;
zlim auto;

%filename=['results' var '.mat'];
save('gaussian_cluster_3D.mat');
