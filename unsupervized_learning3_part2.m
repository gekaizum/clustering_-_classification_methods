function accuracy = evaluate_clustering(labels, estimated_labels, c)
    % Evaluate the accuracy of clustering
    % labels - original cluster labels
    % estimated_labels - estimated cluster labels
    % c - number of clusters
    
    correct = 0;
    for i = 1:c
        % Find the most frequent label in the estimated cluster
        true_label = mode(labels(estimated_labels == i));
        % Count the number of correct predictions
        correct = correct + sum(labels(estimated_labels == i) == true_label);
    end
    % Accuracy as a percentage
    accuracy = (correct / length(labels)) * 100;
end
function clusters = stepwise_optimal_clustering(data, num_clusters)
    % Stepwise Optimal Hierarchical Clustering
    % data - input data points
    % num_clusters - desired number of clusters
    
    [n, ~] = size(data); % Get the number of data points
    clusters = cell(n, 1); % Initialize each point as its own cluster
    for i = 1:n
        clusters{i} = data(i, :);
    end
    
    % Initialize the number of clusters
    current_num_clusters = n;
    
    while current_num_clusters > num_clusters
        min_criterion = inf;
        merge_pair = [1, 2];
        
        % Search for the pair of clusters that results in the smallest change in criterion function
        for i = 1:length(clusters)
            for j = i+1:length(clusters)
                ni = size(clusters{i}, 1);
                nj = size(clusters{j}, 1);
                mi = mean(clusters{i}, 1);
                mj = mean(clusters{j}, 1);
                criterion = sqrt((ni * nj) / (ni + nj)) * norm(mi - mj);
                
                if criterion < min_criterion
                    min_criterion = criterion;
                    merge_pair = [i, j];
                end
            end
        end
        
        % Merge the pair of clusters
        clusters{merge_pair(1)} = [clusters{merge_pair(1)}; clusters{merge_pair(2)}];
        clusters(merge_pair(2)) = []; % Remove the merged cluster
        current_num_clusters = current_num_clusters - 1;
    end
end

% Integration with unsupervized_learning0.m
method_list={'min', 'max', 'avg', 'mean'};
for dimensionNum=2:6 % Choose topology dimension
    for topologyNum=1:4 % Choose topology number
        filename=['gaussian_cluster_', num2str(topologyNum),'_',num2str(dimensionNum),'D', '.mat'];
        load(filename);
        
        % Visualization of Original Clusters
        if dimensionNum == 3
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:,1) == i, :);
                scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(i), 'filled');
                % Original cluster center lines
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
            title({'Original Clusters Visualization', ...
                sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)});
            grid on; 
            view(3); % 3D visualization
            xlim auto;
            ylim auto;
            zlim auto;
        elseif dimensionNum == 2
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:,1) == i, :);
                scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors(i), 'filled');
                % Original cluster center lines
                for j = 1:size(cluster_points, 1)
                    plot([cluster_points(j, 1), cluster_centers(i, 1)], ...
                          [cluster_points(j, 2), cluster_centers(i, 2)], ...
                          colors(i));
                end
            end
            hold off;
            xlabel('X');
            ylabel('Y');
            title({'Original Clusters Visualization', ...
                sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)});
            grid on; 
            xlim auto;
            ylim auto;
        end
        
        for mth=1:length(method_list)
            % Run the stepwise optimal clustering algorithm
            optimal_clusters = stepwise_optimal_clustering(all_vectors, c);
            % Determine estimated cluster labels
            estimated_labels = zeros(size(all_vectors, 1), 1);
            for i = 1:length(optimal_clusters)
                for j = 1:size(optimal_clusters{i}, 1)
                    index = find(ismember(all_vectors, optimal_clusters{i}(j,:), 'rows'));
                    estimated_labels(index) = i;
                end
            end
            
            % Evaluate accuracy
            accuracy = evaluate_clustering(labels, estimated_labels, c);
            fprintf('Accuracy for method %s with topology %d, %dD: %.2f%%\n', method_list{mth}, topologyNum, dimensionNum, accuracy);
            % Visualization of Stepwise Optimal Clusters
            if dimensionNum == 3
                figure;
                hold on;
                for i = 1:length(optimal_clusters)
                    scatter3(optimal_clusters{i}(:, 1), optimal_clusters{i}(:, 2), optimal_clusters{i}(:, 3), 36, colors(i), 'filled');
                    cluster_center = mean(optimal_clusters{i}, 1);
                    for j = 1:size(optimal_clusters{i}, 1)
                        plot3([optimal_clusters{i}(j, 1), cluster_center(1)], ...
                              [optimal_clusters{i}(j, 2), cluster_center(2)], ...
                              [optimal_clusters{i}(j, 3), cluster_center(3)], colors(i));
                    end
                end
                hold off;
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                title({['Stepwise Optimal Clusters Visualization ', method_list{mth}], ...
                       sprintf('Topology %d, %dD', topologyNum, dimensionNum)});
                grid on; 
                view(3); % 3D visualization
                xlim auto;
                ylim auto;
                zlim auto;
            elseif dimensionNum == 2
                figure;
                hold on;
                for i = 1:length(optimal_clusters)
                    scatter(optimal_clusters{i}(:, 1), optimal_clusters{i}(:, 2), 36, colors(i), 'filled');
                    cluster_center = mean(optimal_clusters{i}, 1);
                    for j = 1:size(optimal_clusters{i}, 1)
                        plot([optimal_clusters{i}(j, 1), cluster_center(1)], ...
                              [optimal_clusters{i}(j, 2), cluster_center(2)], ...
                              colors(i));
                    end
                end
                hold off;
                xlabel('X');
                ylabel('Y');
                title({['Stepwise Optimal Clusters Visualization ', method_list{mth}], ...
                       sprintf('Topology %d, %dD', topologyNum, dimensionNum)});
                grid on; 
                xlim auto;
                ylim auto;
            end
        end
    end
end
