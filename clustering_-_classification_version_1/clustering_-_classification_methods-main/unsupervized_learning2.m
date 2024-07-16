%%Evaluation scores
function evaluation_scores = evaluate_clusters(data, labels, centers)
    num_clusters = size(centers, 1);
    evaluation_scores = struct();

    % Davies-Bouldin Index
    db_index = davies_bouldin(data, labels, centers);
    evaluation_scores.db_index = db_index;

    % Silhouette Index
    silhouette_avg = silhouette(data, labels);
    evaluation_scores.silhouette_avg = mean(silhouette_avg);

    % Sum of Squared Errors (SSE)
    sse = sum_of_squared_errors(data, labels, centers);
    evaluation_scores.sse = sse;

    % Within-Cluster Scatter (WCS)
    wcs = within_cluster_scatter(data, labels, centers);
    evaluation_scores.wcs = wcs;

    % Between-Cluster Scatter (BCS)
    bcs = between_cluster_scatter(data, labels, centers);
    evaluation_scores.bcs = bcs;

    % Calinski-Harabasz Index
    ch_index = calinski_harabasz(data, labels, centers);
    evaluation_scores.ch_index = ch_index;
end

function db_index = davies_bouldin(data, labels, centers)
    num_clusters = size(centers, 1);
    db_index = 0;
    for i = 1:num_clusters
        cluster_i = data(labels == i, :);
        s_i = mean(pdist2(cluster_i, centers(i, :)));
        max_ratio = 0;
        for j = 1:num_clusters
            if i ~= j
                cluster_j = data(labels == j, :);
                s_j = mean(pdist2(cluster_j, centers(j, :)));
                d_ij = norm(centers(i, :) - centers(j, :));
                ratio = (s_i + s_j) / d_ij;
                max_ratio = max(max_ratio, ratio);
            end
        end
        db_index = db_index + max_ratio;
    end
    db_index = db_index / num_clusters;
end

function sse = sum_of_squared_errors(data, labels, centers)
    sse = 0;
    for i = 1:size(centers, 1)
        cluster_data = data(labels == i, :);
        sse = sse + sum(sum((cluster_data - centers(i, :)).^2));
    end
end

function wcs = within_cluster_scatter(data, labels, centers)
    wcs = 0;
    for i = 1:size(centers, 1)
        cluster_data = data(labels == i, :);
        wcs = wcs + sum(sum((cluster_data - centers(i, :)).^2));
    end
end

function bcs = between_cluster_scatter(data, labels, centers)
    global_center = mean(data);
    bcs = 0;
    for i = 1:size(centers, 1)
        cluster_size = sum(labels == i);
        bcs = bcs + cluster_size * sum((centers(i, :) - global_center).^2);
    end
end

function ch_index = calinski_harabasz(data, labels, centers)
    num_clusters = size(centers, 1);
    num_points = size(data, 1);
    bcs = between_cluster_scatter(data, labels, centers);
    wcs = within_cluster_scatter(data, labels, centers);
    ch_index = (bcs / (num_clusters - 1)) / (wcs / (num_points - num_clusters));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{'Clusters Visualization with Lines to Cluster Centers', ...
%sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)}
function plot_clusters(data, labels, centers, title_text,topologyNum,dimensionNum)
    if size(data, 2) == 2
        gscatter(data(:,1), data(:,2), labels);
        hold on;
        %plot(centers(:,1), centers(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
        hold off;
        title({title_text, ...
            sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)});
    elseif size(data, 2) == 3
        scatter3(data(:,1), data(:,2), data(:,3), 15, labels, 'filled');
        hold on;
        %plot3(centers(:,1), centers(:,2), centers(:,3), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
        hold off;
        title({title_text, ...
            sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuzzy C-Means Initialization
function [centers, U] = fcm_initialization(data, num_clusters)
    options = [2, 100, 1e-5, 0];%[fuzziness coefficient,maximum number of iterations,minimum amount of improvement,display info]
    [centers, U] = fcm(data, num_clusters, options);
end

% Fuzzy Maximum Likelihood Clustering Method
function [centers, U] = fmcm(data, initial_centers, num_clusters)
    max_iter = 100;
    epsilon = 1e-5;
    centers = initial_centers;
    U = zeros(size(data, 1), num_clusters);
    for iter = 1:max_iter
        for i = 1:num_clusters
            for j = 1:size(data, 1)
                U(j, i) = exp(-norm(data(j, :) - centers(i, :))^2);
            end
        end
        U = U ./ sum(U, 2);
        new_centers = (U' * data) ./ sum(U, 1)';
        if norm(new_centers - centers) < epsilon
            break;
        end
        centers = new_centers;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dimensionNum=2:6 %choose topology dimension
    for topologyNum=1:4 %choose topology number
        filename=['gaussian_cluster_', num2str(topologyNum),'_',num2str(dimensionNum),'D', '.mat'];
        load(filename);
        % Vizualization for 2D and 3D
%{
        if m==3
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
            title({'Clusters Visualization with Lines to Cluster Centers', ...
       sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)});
            grid on; 
            view(3); % 3D vizualization
            xlim auto;
            ylim auto;
            zlim auto;
        end
        if m==2
            figure;
            hold on;
            for i = 1:c
                cluster_points = all_vectors(labels(:,1) == i, :);
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
            title({'Clusters Visualization with Lines to Cluster Centers', ...
       sprintf('gaussian cluster number %d %dD', topologyNum, dimensionNum)});
            grid on; 
            xlim auto;
            ylim auto;
        end
%}
        % Run FCM for initialization
        [initial_centers, ~] = fcm_initialization(all_vectors, c);
        % Run FMCM
        [centers, U] = fmcm(all_vectors, initial_centers, c);
        % Assign labels based on maximum membership
        [~, assigned_labels] = max(U, [], 2);
        % Evaluate clusters
        %evaluation_scores = evaluate_clusters(all_vectors, labels, centers);
        %disp(evaluation_scores);
    
        % Plot original clusters
        if dimensionNum==2 || dimensionNum==3
            figure;
            subplot(1, 2, 1);
            true_labels = labels(:, 1);
            plot_clusters(all_vectors, true_labels, initial_centers, 'Original Clusters',topologyNum,dimensionNum);%
        
            % Plot reconstructed clusters
            subplot(1, 2, 2);
            plot_clusters(all_vectors, assigned_labels, centers, 'Reconstructed Clusters',topologyNum,dimensionNum);
        end
         %%  Check the accuracy of result
        % Calculate centroids of original clusters
        true_labels = labels(:, 1); % Assuming true labels are stored in the first column of `labels`
        original_centroids = zeros(c, m);
        for i = 1:c
            original_centroids(i, :) = mean(all_vectors(true_labels == i, :), 1);
        end
        
        % Calculate centroids of estimated clusters
        estimated_centroids = zeros(c, m);
        for i = 1:c
            estimated_centroids(i, :) = mean(all_vectors(assigned_labels == i, :), 1);
        end
        
        % Assign each estimated cluster to the closest original cluster
        cluster_map = zeros(c, 1);
        for i = 1:c
            distances = zeros(c, 1);
            for j = 1:c
                distances(j) = norm(estimated_centroids(i, :) - original_centroids(j, :));
            end
            [~, cluster_map(i)] = min(distances);
        end
        
        % Compare points in matched clusters
        correct_points = 0;
        for i = 1:c
            true_cluster = original_centroids(cluster_map(i), :);
            estimated_cluster_points = all_vectors(assigned_labels == i, :);
            correct_points = correct_points + sum(true_labels(assigned_labels == i) == find(ismember(original_centroids, true_cluster, 'rows')));
        end
        
        accuracy = correct_points / size(all_vectors, 1);
        disp([filename,' accuracy for unknown mu: ', num2str(accuracy * 100), '%']);
        %%
    end
end