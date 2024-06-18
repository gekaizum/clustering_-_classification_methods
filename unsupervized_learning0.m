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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q 1.1 no mean value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_cluster=[];
[unique_values, idx, ~] = unique(labels(:,1), 'stable');
merged_labels = labels(idx, :); %matrix of parameters: cluster num/size/mean/variance

%% Function checks probability that point belongs to current cluster and returns probabilities
function res = pdf_of_cluster(x_val,y_val,z_val,mu_val,sigma_val)
    sigma_val=sigma_val^2;
    answer_x=(1/sqrt(2*pi*sigma_val)*exp((-(x_val-mu_val(1))^2)/(2*sigma_val)));
    if answer_x < 1e-200 % Set a minimum threshold
        answer_x = 1e-200;
    end
    answer_y=(1/sqrt(2*pi*sigma_val)*exp((-(y_val-mu_val(2))^2)/(2*sigma_val)));
    if answer_y < 1e-200 % Set a minimum threshold
        answer_y = 1e-200;
    end
    answer_z=(1/sqrt(2*pi*sigma_val)*exp((-(z_val-mu_val(3))^2)/(2*sigma_val)));
    if answer_z < 1e-200 % Set a minimum threshold
        answer_z = 1e-200;
    end
    res(1,:) = [answer_x, answer_y, answer_z];

end
%% Function checks probability that point belongs to current cluster and returns 1 or 0
function result = estimation_of_cluster(x_val,y_val,z_val,mu_val,sigma_val)
    sigma_val=sigma_val^2;
    answer_x=(1/sqrt(2*pi*sigma_val)*exp((-(x_val-mu_val)^2)/(2*sigma_val)));
    if answer_x < 1e-200 % Set a minimum threshold
        answer_x = 1e-200;
    end
    answer_y=(1/sqrt(2*pi*sigma_val)*exp((-(y_val-mu_val)^2)/(2*sigma_val)));
    if answer_y < 1e-200 % Set a minimum threshold
        answer_y = 1e-200;
    end
    answer_z=(1/sqrt(2*pi*sigma_val)*exp((-(z_val-mu_val)^2)/(2*sigma_val)));
    if answer_z < 1e-200 % Set a minimum threshold
        answer_z = 1e-200;
    end
    if (answer_x < (0.1) && answer_y < (0.1))  || (answer_x < (0.1) && answer_z < (0.1)) || (answer_y  < (0.1) && answer_z  < (0.1))
        result=0;
    else
        result=1;
    end
end
%%
aprior_prob=[merged_labels(:,1)];
aprior_prob=[aprior_prob,merged_labels(:,2)/sum(merged_labels(:,2))];

mean_vector = repmat(mean(all_vectors, 1),c,1); %init vector of estimated mean values
tolerance = 0.1;
% Tolerance matrix
random_offsets = tolerance * merged_labels(:,4) .* randn(c, size(all_vectors, 2));
mean_vector = mean_vector + random_offsets;

maxIter=100; % number of iterations of search for mean values

%% MLE algorithm for finding mean value
for iterNumber = 1:maxIter
    new_mean_vector=zeros(c,3);
    for i = 1:c
        mone=zeros(1,3);
        mahane=zeros(1,3);
        for j = 1: size(all_vectors,1)
            estimProb_i=pdf_of_cluster(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(i,:),merged_labels(i,4));
            estimProb_i=estimProb_i*aprior_prob(i,2);
            total_estimProb_i=zeros(1,3);
            for try_mean = 1:c
                temp=pdf_of_cluster(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(try_mean,:),merged_labels(try_mean,4));
                temp=temp*aprior_prob(try_mean,1);
                total_estimProb_i=total_estimProb_i+temp;
            end
            end_prob=estimProb_i./total_estimProb_i;
            mone=mone+(end_prob.*all_vectors(j,:));
            mahane=mahane+end_prob;
        end
        new_mean_vector(i,:)=mone./mahane;
    end
    mean_vector=new_mean_vector;
end
%% Estimation of clusters based on estimated mean values
vector_matrix_vizual=all_vectors;
index=0;
figure;
hold on;
current_cluster_centers=[];
list_of_estimated_clusters= struct('Matrix', cell(c, 1), 'ClusterMean', cell(c, 1));
for i = 1:c
    estimated_cluster=[];
        for j =1:size(vector_matrix_vizual,1) %loop must find vectors which are belong to same cluster i
            %for search we use function which checks PDF value for current x
            if estimation_of_cluster(vector_matrix_vizual(j,1),vector_matrix_vizual(j,2),vector_matrix_vizual(j,3),mean(mean_vector(i,:)),merged_labels(i,4)) == 1
                estimated_cluster=[estimated_cluster; vector_matrix_vizual(j,:)];
            end
        end
        if size(estimated_cluster)>0
            current_cluster_centers = [current_cluster_centers; mean(estimated_cluster, 1)];
            [~, rows_to_remove] = ismember(estimated_cluster, vector_matrix_vizual, 'rows');
            % Remove rows from A that are in B
            vector_matrix_vizual(rows_to_remove, :) = [];
        else
            current_cluster_centers = [current_cluster_centers; [-666, -666, -666]];
        end
        
        list_of_estimated_clusters(i).Matrix = estimated_cluster;
        list_of_estimated_clusters(i).ClusterCenter = mean_vector(i,:);
end


%% If some of points doesnt belongs to any cluster we will assign them here
for i = 1:size(vector_matrix_vizual)
    distances = zeros(c, 1);
    for j = 1:c
        distances(j) = norm(vector_matrix_vizual(i, :) - mean_vector(j,:));
    end
    [~, cluster_idx] = min(distances);
    list_of_estimated_clusters(cluster_idx).Matrix = vertcat(list_of_estimated_clusters(cluster_idx).Matrix, vector_matrix_vizual(i, :));
    current_cluster_centers(cluster_idx,:)= mean(list_of_estimated_clusters(cluster_idx).Matrix, 1);
end
%% Vizualization
for i = 1:c
    if size(list_of_estimated_clusters(i).Matrix)>0
        scatter3(list_of_estimated_clusters(i).Matrix(:, 1), list_of_estimated_clusters(i).Matrix(:, 2), list_of_estimated_clusters(i).Matrix(:, 3), 36, colors(i), 'filled');
        % Grid line from center to points in cluster
        for j = 1:size(list_of_estimated_clusters(i).Matrix, 1)
            plot3([list_of_estimated_clusters(i).Matrix(j, 1), current_cluster_centers(i-index, 1)], ...
                  [list_of_estimated_clusters(i).Matrix(j, 2), current_cluster_centers(i-index, 2)], ...
                  [list_of_estimated_clusters(i).Matrix(j, 3), current_cluster_centers(i-index, 3)], colors(i));
        end
    end
end

hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Estimated Clusters with estimated mean value');
grid on;
if m==3
    view(3); % 3D vizualization
end
xlim auto;
ylim auto;
zlim auto;

%%  Check the accuracy of result

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q 1.2 only amount of clusters known %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = pdf_of_cluster2(x_val,y_val,z_val,mu_val,sigma_val)
    answer_x=(1/sqrt(2*pi*sigma_val(1)^2)*exp((-(x_val-mu_val(1))^2)/(2*sigma_val(1)^2)));
    if answer_x < 1e-200 % Set a minimum threshold
        answer_x = 1e-200;
    end
    answer_y=(1/sqrt(2*pi*sigma_val(2)^2)*exp((-(y_val-mu_val(2))^2)/(2*sigma_val(2)^2)));
    if answer_y < 1e-200 % Set a minimum threshold
        answer_y = 1e-200;
    end
    answer_z=(1/sqrt(2*pi*sigma_val(3)^2)*exp((-(z_val-mu_val(3))^2)/(2*sigma_val(3)^2)));
    if answer_z < 1e-200 % Set a minimum threshold
        answer_z = 1e-200;
    end
    res(1,:) = [answer_x, answer_y, answer_z];

end

mean_vector = repmat(mean(all_vectors, 1),c,1); %init vector of estimated mean values
deviation_vector=repmat(var(all_vectors,1)/100,c,1); %init vector of estimated standart deviation values
estim_aprior_prob=repmat((sum(merged_labels(:,2))/c)/sum(merged_labels(:,2)),c,1);%init vector of estimated aprioric probabilities
new_aprior_prob=zeros(c,1);

tolerance = 0.1;

% Tolerance matrix
random_offsets = tolerance * deviation_vector .* randn(c, size(all_vectors, 2));
mean_vector = mean_vector + random_offsets;
deviation_vector=deviation_vector + random_offsets;

maxIter=200;
for iterNumber = 1:maxIter
    %% Aprioric probability estimation
    for i = 1:c
        end_prob=zeros(1,3);
        for j = 1: size(all_vectors,1)
            estimProb_i=pdf_of_cluster2(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(i,:),deviation_vector(i,:));
            estimProb_i=estimProb_i*estim_aprior_prob(i,1);
            total_estimProb_i=zeros(1,3);
            for try_mean = 1:c
                temp=pdf_of_cluster2(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(try_mean,:),deviation_vector(try_mean,:));
                temp=temp*estim_aprior_prob(try_mean,1);
                total_estimProb_i=total_estimProb_i+temp;
            end
            end_prob=end_prob+(estimProb_i./total_estimProb_i);
        end
        new_aprior_prob(i)=mean(end_prob/size(all_vectors,1));
    end
    total_prob = sum(new_aprior_prob);
    if total_prob > 1
        new_aprior_prob = new_aprior_prob / total_prob;
    end
    
    %% Mean value eatimation
    new_mean_vector=zeros(c,3);
    for i = 1:c
        end_prob=0;
        mone=zeros(1,3);
        mahane=zeros(1,3);
        for j = 1: size(all_vectors,1)
            estimProb_i=pdf_of_cluster2(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(i,:),deviation_vector(i,:));
            estimProb_i=estimProb_i*estim_aprior_prob(i,1);
            total_estimProb_i=zeros(1,3);
            for try_mean = 1:c
                temp=pdf_of_cluster2(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(try_mean,:),deviation_vector(try_mean,:));
                temp=temp*estim_aprior_prob(try_mean,1);
                total_estimProb_i=total_estimProb_i+temp;
            end
            end_prob=estimProb_i./total_estimProb_i;
            mone=mone+(end_prob.*all_vectors(j,:));
            mahane=mahane+end_prob;
        end
        new_mean_vector(i,:)=mone./mahane;
    end
    
    %% Variance value eatimation
    new_deviation_vector=zeros(c,3);
    for i = 1:c
        end_prob=0;
        mone=zeros(1,3);
        mahane=zeros(1,3);
        for j = 1: size(all_vectors,1)
            estimProb_i=pdf_of_cluster2(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(i,:),deviation_vector(i,:));
            estimProb_i=estimProb_i*estim_aprior_prob(i,1);
            total_estimProb_i=zeros(1,3);
            for try_mean = 1:c
                temp=pdf_of_cluster2(all_vectors(j,1),all_vectors(j,2),all_vectors(j,3),mean_vector(try_mean,:),deviation_vector(try_mean,:));
                temp=temp*estim_aprior_prob(try_mean,1);
                total_estimProb_i=total_estimProb_i+temp;
            end
            end_prob=estimProb_i./total_estimProb_i;
            mone=mone+(end_prob.*((all_vectors(j,:)-mean_vector(i,:)).^2));
            mahane=mahane+end_prob;
        end
        new_deviation_vector(i,:)=mone./mahane;
    end
    mean_vector=new_mean_vector;
    estim_aprior_prob=new_aprior_prob;
    deviation_vector=new_deviation_vector;
end

vector_matrix_vizual=all_vectors;
index=0;
fig=figure;
hold on;
current_cluster_centers=[];
list_of_estimated_clusters= struct('Matrix', cell(c, 1), 'ClusterMean', cell(c, 1));
for i = 1:c
    estimated_cluster=[];
        for j =1:size(vector_matrix_vizual,1) %loop must find vectors which are belong to same cluster i
            %for search we use function which check PDF value for current x
            if estimation_of_cluster(vector_matrix_vizual(j,1),vector_matrix_vizual(j,2),vector_matrix_vizual(j,3),mean(mean_vector(i,:)),mean(deviation_vector(i,:))) == 1
                estimated_cluster=[estimated_cluster; vector_matrix_vizual(j,:)];
            end
        end
        if size(estimated_cluster)>0
            current_cluster_centers = [current_cluster_centers; mean(estimated_cluster, 1)];
            [~, rows_to_remove] = ismember(estimated_cluster, vector_matrix_vizual, 'rows');
            % Remove rows from A that are in B
            vector_matrix_vizual(rows_to_remove, :) = [];
        else
            current_cluster_centers = [current_cluster_centers; [-666, -666, -666]];
        end
        
        list_of_estimated_clusters(i).Matrix = estimated_cluster;
        list_of_estimated_clusters(i).ClusterCenter = mean_vector(i,:);
end


% If some of points doesnt belongs to any cluster we will assign them here
for i = 1:size(vector_matrix_vizual)
    distances = zeros(c, 1);
    for j = 1:c
        distances(j) = norm(vector_matrix_vizual(i, :) - mean_vector(j,:));
    end
    [~, cluster_idx] = min(distances);
    list_of_estimated_clusters(cluster_idx).Matrix = vertcat(list_of_estimated_clusters(cluster_idx).Matrix, vector_matrix_vizual(i, :));
    current_cluster_centers(cluster_idx,:)= mean(list_of_estimated_clusters(cluster_idx).Matrix, 1);
end

for i = 1:c
    if size(list_of_estimated_clusters(i).Matrix)>0
        scatter3(list_of_estimated_clusters(i).Matrix(:, 1), list_of_estimated_clusters(i).Matrix(:, 2), list_of_estimated_clusters(i).Matrix(:, 3), 36, colors(i), 'filled');
        % Grid line from center to points in cluster
        for j = 1:size(list_of_estimated_clusters(i).Matrix, 1)
            plot3([list_of_estimated_clusters(i).Matrix(j, 1), current_cluster_centers(i-index, 1)], ...
                  [list_of_estimated_clusters(i).Matrix(j, 2), current_cluster_centers(i-index, 2)], ...
                  [list_of_estimated_clusters(i).Matrix(j, 3), current_cluster_centers(i-index, 3)], colors(i));
        end
    end
end

hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Estimated Clusters with estimated mean, variance and aprioric probability');
grid on;
if m==3
    view(3); % 3D vizualization
end
xlim auto;
ylim auto;
zlim auto;

%%  Check the accuracy of result
estimated_merged_labels=round(estim_aprior_prob*size(all_vectors,1));
estimated_merged_labels=[estimated_merged_labels,mean(mean_vector,2)];
estimated_merged_labels=[estimated_merged_labels,mean(deviation_vector,2)]
%%