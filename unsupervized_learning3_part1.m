function dist = calc_distance(D1, D2, method)
    % Calculate the distance between two clusters D1 and D2 based on the specified method
    switch method
        case 'min'
            % Minimum distance between any two points from the clusters
            dist = min(pdist2(D1, D2, 'euclidean'), [], 'all');
        case 'max'
            % Maximum distance between any two points from the clusters
            dist = max(pdist2(D1, D2, 'euclidean'), [], 'all');
        case 'avg'
            % Average distance between all points from the clusters
            dist = mean(pdist2(D1, D2, 'euclidean'), 'all');
        case 'mean'
            % Euclidean distance between the centroids of the clusters
            m1 = mean(D1, 1);
            m2 = mean(D2, 1);
            dist = norm(m1 - m2);
        otherwise
            error('Unknown method.');
    end
end

function clusters = agglomerative_clustering(data, num_clusters, method)
    % Perform agglomerative clustering on the data
    % data - input data
    % num_clusters - desired number of clusters
    % method - method to calculate distances ('min', 'max', 'avg', 'mean')

    % Initialization
    [n, ~] = size(data);
    clusters = cell(n, 1);
    for i = 1:n
        clusters{i} = data(i, :);  % Each data point starts as its own cluster
    end

    % Clustering process
    while length(clusters) > num_clusters
        min_dist = inf;  % Initialize minimum distance to a very large value
        merge_pair = [1, 2];  % Initialize the pair of clusters to be merged

        % Find the closest pair of clusters
        for i = 1:length(clusters)
            for j = i+1:length(clusters)
                dist = calc_distance(clusters{i}, clusters{j}, method);
                if dist < min_dist
                    min_dist = dist;
                    merge_pair = [i, j];
                end
            end
        end

        % Merge the closest pair of clusters
        clusters{merge_pair(1)} = [clusters{merge_pair(1)}; clusters{merge_pair(2)}];
        clusters(merge_pair(2)) = [];  % Remove the merged cluster
    end
end

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

method_list = {'min', 'max', 'avg', 'mean'};  % List of methods
colors = 'rgbcmyk';  % Colors for clusters

for dimensionNum = 2:6  % Loop over different dimensions
    for topologyNum = 1:4  % Loop over different topologies
        filename = ['gaussian_cluster_', num2str(topologyNum), '_', num2str(dimensionNum), 'D', '.mat'];
        load(filename);

        % Visualization of original clusters
        if dimensionNum == 3
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:,1) == i, :);
                scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(i), 'filled');
                % Draw lines from the cluster center to the points
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
            title({'Original Clusters', sprintf('Topology %d, %dD', topologyNum, dimensionNum)});
            grid on;
            view(3);
        elseif dimensionNum == 2
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:,1) == i, :);
                scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors(i), 'filled');
                % Draw lines from the cluster center to the points
                for j = 1:size(cluster_points, 1)
                    plot([cluster_points(j, 1), cluster_centers(i, 1)], ...
                          [cluster_points(j, 2), cluster_centers(i, 2)], colors(i));
                end
            end
            hold off;
            xlabel('X');
            ylabel('Y');
            title({'Original Clusters', sprintf('Topology %d, %dD', topologyNum, dimensionNum)});
            grid on;
        end

        % Visualization of agglomerative clustering results
        for mth = 1:length(method_list)
            % Perform clustering
            clusters = agglomerative_clustering(all_vectors, c, method_list{mth});
            
            % Determine estimated cluster labels
            estimated_labels = zeros(size(all_vectors, 1), 1);
            for i = 1:length(clusters)
                for j = 1:size(clusters{i}, 1)
                    index = find(ismember(all_vectors, clusters{i}(j,:), 'rows'));
                    estimated_labels(index) = i;
                end
            end
            
            % Evaluate accuracy
            accuracy = evaluate_clustering(labels, estimated_labels, c);
            fprintf('Accuracy for method %s with topology %d, %dD: %.2f%%\n', method_list{mth}, topologyNum, dimensionNum, accuracy);

            % Visualization of clustering results
            if dimensionNum == 3
                figure;
                hold on;
                for i = 1:length(clusters)
                    cluster_points = clusters{i};
                    scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(i), 'filled');
                    % Draw lines from the cluster center to the points
                    cluster_center = mean(cluster_points, 1);
                    for j = 1:size(cluster_points, 1)
                        plot3([cluster_points(j, 1), cluster_center(1)], ...
                              [cluster_points(j, 2), cluster_center(2)], ...
                              [cluster_points(j, 3), cluster_center(3)], colors(i));
                    end
                end
                hold off;
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                title({['Agglomerative Clustering with Method ', method_list{mth}], ...
                       sprintf('Topology %d, %dD', topologyNum, dimensionNum)});
                grid on;
                view(3);
            elseif dimensionNum == 2
                figure;
                hold on;
                for i = 1:length(clusters)
                    cluster_points = clusters{i};
                    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors(i), 'filled');
                    % Draw lines from the cluster center to the points
                    cluster_center = mean(cluster_points, 1);
                    for j = 1:size(cluster_points, 1)
                        plot([cluster_points(j, 1), cluster_center(1)], ...
                              [cluster_points(j, 2), cluster_center(2)], colors(i));
                    end
                end
                hold off;
                xlabel('X');
                ylabel('Y');
                title({['Agglomerative Clustering with Method ', method_list{mth}], ...
                       sprintf('Topology %d, %dD', topologyNum, dimensionNum)});
                grid on;
            end
        end
    end
end
