load('gaussian_cluster_2_3D.mat');
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

%% Gaussian PDF for vectors
function res = pdf_of_cluster(cur_vector,mu_vector, sigma_val)
    cov_vector=zeros(length(cur_vector));
    cov_vector(1:length(cur_vector)+1:end) = sigma_val^2;
    answer=(1/sqrt(2*pi*sqrt(det(cov_vector)))*exp((-(1/2)*(cur_vector-mu_vector)*inv(cov_vector)*transpose(cur_vector-mu_vector))));
    if answer < 1e-200 % Set a minimum threshold
        answer = 1e-200;
    end
    res = answer;
end
%%
aprior_prob=[merged_labels(:,1)];
aprior_prob=[aprior_prob,merged_labels(:,2)/sum(merged_labels(:,2))];

mean_vector=repmat(mean(mean(all_vectors, 1)),c,m); %init vector of estimated mean values
tolerance = 0.1;
% Tolerance matrix
random_offsets = tolerance * merged_labels(:,4) .* randn(c, 1);
mean_vector = mean_vector + random_offsets;

maxIter=100; % number of iterations of search for mean values

%% MLE algorithm for finding mean value
for iterNumber = 1:maxIter
    new_mean_vector=zeros(c,3);
    for i = 1:c
        mone=zeros(1,3);
        mahane=zeros(1,3);
        for j = 1: size(all_vectors,1)
            estimProb_i=pdf_of_cluster(all_vectors(j,:),mean_vector(i,:), merged_labels(i, 4));
            estimProb_i=estimProb_i*aprior_prob(i,2);
            total_estimProb_i=zeros(1,3);
            for try_mean = 1:c
                temp=pdf_of_cluster(all_vectors(j,:),mean_vector(try_mean,:), merged_labels(i, 4));
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
%% Bayesian Decision Theory
assigned_labels = zeros(size(all_vectors, 1), 1);
for j = 1:size(all_vectors, 1)
    posterior_probs = zeros(c, 1);
    for i = 1:c
        posterior_probs(i) = prod(pdf_of_cluster(all_vectors(j,:), mean_vector(i, :), merged_labels(i, 4))) * aprior_prob(i, 2);
    end
    [~, assigned_labels(j)] = max(posterior_probs);
end

% Visualization of the result after Bayesian Decision Theory
figure;
hold on;
for i = 1:c
    cluster_points = all_vectors(assigned_labels == i, :);
    scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(i), 'filled');
    % Grid line from center to points in cluster
    for j = 1:size(cluster_points, 1)
        plot3([cluster_points(j, 1), mean_vector(i, 1)], ...
              [cluster_points(j, 2), mean_vector(i, 2)], ...
              [cluster_points(j, 3), mean_vector(i, 3)], colors(i));
    end
end
hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Clusters Visualization with Bayesian Decision Theory and Grid Lines');
grid on;
if m == 3
    view(3); % 3D visualization
end
xlim auto;
ylim auto;
zlim auto;

%%  Check the accuracy of result

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q 1.2 only amount of clusters known %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
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

%% Bayesian Decision Theory
assigned_labels = zeros(size(all_vectors, 1), 1);
for j = 1:size(all_vectors, 1)
    posterior_probs = zeros(c, 1);
    for i = 1:c
        posterior_probs(i) = prod(pdf_of_cluster(all_vectors(j, 1), all_vectors(j, 2), all_vectors(j, 3), mean_vector(i, :), merged_labels(i, 4))) * aprior_prob(i, 2);
    end
    [~, assigned_labels(j)] = max(posterior_probs);
end

% Visualization of the result after Bayesian Decision Theory
figure;
hold on;
for i = 1:c
    cluster_points = all_vectors(assigned_labels == i, :);
    scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(i), 'filled');
    % Grid line from center to points in cluster
    for j = 1:size(cluster_points, 1)
        plot3([cluster_points(j, 1), mean_vector(i, 1)], ...
              [cluster_points(j, 2), mean_vector(i, 2)], ...
              [cluster_points(j, 3), mean_vector(i, 3)], colors(i));
    end
end
hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Clusters Visualization with Bayesian Decision Theory and Grid Lines');
grid on;
if m == 3
    view(3); % 3D visualization
end
xlim auto;
ylim auto;
zlim auto;

%%  Check the accuracy of result
estimated_merged_labels=round(estim_aprior_prob*size(all_vectors,1));
estimated_merged_labels=[estimated_merged_labels,mean(mean_vector,2)];
estimated_merged_labels=[estimated_merged_labels,mean(deviation_vector,2)]
%%
%}