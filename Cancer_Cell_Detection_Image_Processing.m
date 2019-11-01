Detection of cancer cells over blood microscopy images based on shape anomaly

Samples of 20 images’ segmentation
ab_global=[]
for image_index=1:5:108
    file_name=dir(strcat('~/im/*.jpg'));
    input_image=imread(strcat('~/im/',file_name(image_index).name));
    %% L*a*b Color Space
    inputImLAB = rgb2lab(input_image);
    % Extract a* and b* channels and reshape
    ab = double(inputImLAB(:,:,2:3));
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,nrows*ncols,2);
    ab_global = [ab_global;ab];
end
%% Segmentation usign K-Means
nColors = 6;
[cluster_idx, cluster_center] = kmeans(ab_global,nColors) 

Image Segmentation
load computedClusterCenters.mat

stats_total = [];
file_name=dir(strcat('~/im/*.jpg'));
cancer_cell_coords=dir(strcat('~/edited_coords/*.xyc.txt'));
for main_loop_index=1:108
    input_image=imread(strcat('~/im/',file_name(main_loop_index).name));

%% Image segmentation
    size_image = size(input_image);
    n_cols = size_image(1);
    n_rows = size_image(2);
    inputImLAB = rgb2lab(input_image);
    
    %% for k=6 cluster matrix created for every a and b values 
    cluster1_ab = zeros([n_cols, n_rows, 2]);
    cluster1_ab(:,:,1) = cluster1_ab(:,:,1) + cluster_center(1,1);
    cluster1_ab(:,:,2) = cluster1_ab(:,:,2) + cluster_center(1,2);
    cluster2_ab = zeros([n_cols, n_rows, 2]);
    cluster2_ab(:,:,1) = cluster2_ab(:,:,1) + cluster_center(2,1);
    cluster2_ab(:,:,2) = cluster2_ab(:,:,2) + cluster_center(2,2);
    cluster3_ab = zeros([n_cols, n_rows, 2]);
    cluster3_ab(:,:,1) = cluster3_ab(:,:,1) + cluster_center(3,1);
    cluster3_ab(:,:,2) = cluster3_ab(:,:,2) + cluster_center(3,2);
    cluster4_ab = zeros([n_cols, n_rows, 2]);
    cluster4_ab(:,:,1) = cluster4_ab(:,:,1) + cluster_center(4,1);
    cluster4_ab(:,:,2) = cluster4_ab(:,:,2) + cluster_center(4,2);
    cluster5_ab = zeros([n_cols, n_rows, 2]);
    cluster5_ab(:,:,1) = cluster5_ab(:,:,1) + cluster_center(5,1);
    cluster5_ab(:,:,2) = cluster5_ab(:,:,2) + cluster_center(5,2);
    cluster6_ab = zeros([n_cols, n_rows, 2]);
    cluster6_ab(:,:,1) = cluster6_ab(:,:,1) + cluster_center(6,1);
    cluster6_ab(:,:,2) = cluster6_ab(:,:,2) + cluster_center(6,2);

    %% ab channel of input image
    ab = double(inputImLAB(:,:,2:3));
 
    %% Calculation Euclidian distance for every channel for 6 cluster
    diff_to_cluster_1 = ab - cluster1_ab;
    distance_cluster_1 = sqrt(diff_to_cluster_1(:,:,1).^2 + diff_to_cluster_1(:,:,2).^2);
    
    diff_to_cluster_2 = ab - cluster2_ab;
    distance_cluster_2 = sqrt(diff_to_cluster_2(:,:,1).^2 + diff_to_cluster_2(:,:,2).^2);
    
    diff_to_cluster_3 = ab - cluster3_ab;
    distance_cluster_3 = sqrt(diff_to_cluster_3(:,:,1).^2 + diff_to_cluster_3(:,:,2).^2);
    
    diff_to_cluster_4 = ab - cluster4_ab;
    distance_cluster_4 = sqrt(diff_to_cluster_4(:,:,1).^2 + diff_to_cluster_4(:,:,2).^2);
    
    diff_to_cluster_5 = ab - cluster5_ab;
    distance_cluster_5 = sqrt(diff_to_cluster_5(:,:,1).^2 + diff_to_cluster_5(:,:,2).^2);
    diff_to_cluster_6 = ab - cluster6_ab;
    distance_cluster_6 = sqrt(diff_to_cluster_6(:,:,1).^2 + diff_to_cluster_6(:,:,2).^2);

    %% Every distance converted to the matrix
    dist_3D = zeros([n_cols, n_rows 6]);
    dist_3D(:, :, 1) = distance_cluster_1;
    dist_3D(:, :, 2) = distance_cluster_2;
    dist_3D(:, :, 3) = distance_cluster_3;
    dist_3D(:, :, 4) = distance_cluster_4;
    dist_3D(:, :, 5) = distance_cluster_5;
    dist_3D(:, :, 6) = distance_cluster_6;
    
    %%Finding minimum distance of every matrix
    [min_data_3d, min_data_3d_ind] = min(dist_3D, [], 3);
    min_data_3d_ind(min_data_3d_ind ~= 4) = 0;
    min_data_3d_ind(min_data_3d_ind == 4) = 1;
 
Morphological Operations
    %% Morphological operations to remove small regions than 15 and convert binary image
    segmented_image = imbinarize(min_data_3d_ind);  
    se = strel('disk',15);

%% Watershed Algorithm implementation
    segmented_image = imopen(segmented_image,se);
    bw2=bwareaopen(~segmented_image, 10);
    D = -bwdist(~segmented_image,'euclidean');
    Ld=watershed(D);
    bw2 = segmented_image;
    bw2(Ld == 0) = 0;
    mask = imextendedmin(D,2);
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw3 = segmented_image;
    bw3(Ld2 == 0) = 0;

%%Connecting component labelling
    CC = bwconncomp(bw3,8);
    
    stats = regionprops(CC,'Area', 'BoundingBox', 'Centroid','ConvexArea', 'Eccentricity', 'EquivDiameter','EulerNumber','Extent','FilledArea','Perimeter', 'MajorAxisLength', 'MinorAxisLength', 'Image','Solidity');
    [stats(:).CancerCell] = deal(0);
    centers = reshape([stats.Centroid], 2, CC.NumObjects);
    majors = [stats.MajorAxisLength];
    minors = [stats.MinorAxisLength];
    diameters= [stats.EquivDiameter];
    radii = [diameters./2];

%%Normalization and detecting sub-images
    for i=1:CC.NumObjects
        %imshow(stats(i).Image);
        bb = stats(i).BoundingBox;
        region_img = rgb2gray(imcrop(input_image, bb));
        %% Normalization
        region_img = uint8(255*mat2gray(region_img));
    end
end
 
Feature Extraction
        offsets = [0 1; -1 1;-1 0;-1 -1];
        [glcms,SI] = graycomatrix(region_img,'Offset', offsets); 
        %% Entropy
        ent_result= entropy(region_img);
        %% Standard deviation of matrix elements
        std_dev = std2(region_img);
        %% Skewness and Kurtosis
        [pixelCount, GLs] = imhist(region_img);
        [skew kurtosis] = GetSkewAndKurtosis(GLs, pixelCount);
        %% mean of image
        meanval = mean2(region_img);
        out=graycoprops(glcms);
        % Add calculated features to stats
        stats(i).Contrast = out.Contrast;
        stats(i).Correlation = out.Correlation;
        stats(i).Energy = out.Energy;
        stats(i).Homogeneity = out.Homogeneity;
        stats(i).ent = ent_result;
        stats(i).std_dev  =std_dev ;
        stats(i).skew  = skew ;
        stats(i).kurtosis  =kurtosis ;
        stats(i).meanval = meanval;
        %%Checking labels of images
        if contains(file_name(main_loop_index).name, '_1')
            cells=tdfread(strcat('/~/edited_coords/'
));
            x1 = bb(1);
            y1 = bb(2);
            width = bb(3);
            height = bb(4);
            x2 = x1 + width;
            y2 = y1 - height;
            for j=1:length(cells.X)
                if x1 < cells.X(j) && cells.X(j) < x2 && y1 < cells.Y(j) && cells.Y(j) > y2
                    stats(i).CancerCell = 1;
                    
                end
            end
        end
        
    end
    
    stats_total = [stats; stats_total];
 
Classification
TEST_DATA_PERCENTAGE = 40;
ITERATION_COUNT = 5;
 
features = csvread('~/test_data_features_extracted.csv',1, 0);
 
overall_success = zeros(ITERATION_COUNT);
 
for i=1:ITERATION_COUNT
    shuffledArray = features(randperm(size(features,1)),:);
    % Calculate 40% of training data size as test data to calculate confidence of model
    total_training_set_size = size(shuffledArray,1);
    test_data_set_size = round(total_training_set_size/ TEST_DATA_PERCENTAGE);
    % Shuffle the all training data so we can pick %40 of rows distributed
    % randomly from all images
    training_data = shuffledArray(1:end-test_data_set_size + 1,:);
    % Cut %60 of total data to be used to train the model
    training_data_features = training_data(:,1:end-1);
    training_data_labels = training_data(:,end);
    % Cut 40% of total data to be used to test the model confidence
    test_data = shuffledArray(end-test_data_set_size:end,:);
    test_data_features = test_data(:,1:end-1);
    test_data_labels = test_data(:,end);
 
    % Train model
    trained_model = fitcsvm(training_data_features,training_data_labels,'Standardize',true,'KernelFunction','RBF',...
  'KernelScale', 'auto');
    % Test model
    test_result = predict(trained_model,test_data_features);
    confMat_svm = confusionmat(test_data_labels,test_result);
    confMat_svm = bsxfun(@rdivide,confMat_svm,sum(confMat_svm,2));
    disp(confMat_svm);
    fprintf('confmat-Completed %d iteration with %s success rate.\n',i ,num2str(mean(diag(confMat_svm)), '%.2f'));
end

%Multilayer Perceptron
training_data_features_neuralnetwork = zscore(features(:,1:end-1));
training_data_labels_neuralnetwork = features(:,end);
inputs = training_data_features_neuralnetwork';
targets = training_data_labels_neuralnetwork';
 
% Create a Feedforward with Levenberg-Marquardt backpropagation
net = feedforwardnet(6, 'trainlm');
% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 60/100;
net.divideParam.valRatio = 20/100;
net.divideParam.testRatio = 20/100;
 
% Train the Network
[net,tr] = train(net,inputs,targets);
% Test the Network
outputs = net(inputs); % predict()
