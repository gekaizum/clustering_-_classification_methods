%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q 1.0 Creating dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%m = 3; % length of vector (dimension)
c = 6; % number of clusters
colors = ['r', 'g', 'b','c','m','o','k']; % colors for clusters

for m = 2:6
    for clusnum = 1:4
        % Init matrix and vectors
        all_vectors = [];
        labels = [];
        cluster_centers = [];
        for i = 1:c
            n = randi([10,100]); % number of points in cluster
            mean_val = randi([-20,20], 1, m); % mean values for cluster
            cov_matrix = wishrnd(eye(m), m); % random covariance matrix using Wishart distribution

            % Creating vector of cluster
            cluster_vectors = mvnrnd(mean_val, cov_matrix, n);

            % Adding cluster to space matrix
            all_vectors = [all_vectors; cluster_vectors];
            % Adding labels for clusters
            labels = [labels; repmat([i, n, mean_val, cov_matrix(:)'], n, 1)]; % cluster ID / number of points in cluster / mean values / variances

            % Calculating center of cluster
            cluster_centers = [cluster_centers; mean(cluster_vectors, 1)];
        end
        
        % Vizualization
        if m == 3
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:, 1) == i, :);
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
            view(3); % 3D visualization
            xlim auto;
            ylim auto;
            zlim auto;
        end
        if m == 2
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:, 1) == i, :);
                scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors(i), 'filled');
                % Grid line from center to points in cluster
                for j = 1:size(cluster_points, 1)
                    plot([cluster_points(j, 1), cluster_centers(i, 1)], ...
                          [cluster_points(j, 2), cluster_centers(i, 2)], ...
                          colors(i));
                end
            end
            hold off;
            xlabel('X');
            ylabel('Y');
            title('Clusters Visualization with Lines to Cluster Centers');
            grid on; 
            xlim auto;
            ylim auto;
        end

        filename = ['gaussian_cluster_', num2str(clusnum), '_', num2str(m), 'D', '.mat'];
        save(filename);
    end
end
