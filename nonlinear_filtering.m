% read in the image
im = imread('F:\PSU第一学期相关书籍及其工作\Digtal Image Processing II\Project 3\disk.gif');
im = mat2gray(im, double([min(min(im)) max(max(im))]));
figure, imshow(im)
title('Original')

%% Median filters
% 1 iteration
medfilt1_im = medfilt2(im, [5 5]);
figure, imshow(medfilt1_im)
title('Median Filter: 1 Iteration')

% histogram: 1 iteration
figure, imhist(medfilt1_im)
title('Histogram: Median filter, 1 iteration')

% statistics: 1 iteration
medfilt1_subset = medfilt1_im(100:150,50:100);
medfilt1_mean = mean(medfilt1_subset(:));
medfilt1_stdev = std(medfilt1_subset(:));

% 5 iterations
medfilt5_im = im;
for i = 1:5
    medfilt5_im = medfilt2(medfilt5_im);
end
figure, imshow(medfilt5_im)
title('Median Filter: 5 Iterations')

% histogram: 5 iterations
figure, imhist(medfilt5_im)
title('Histogram: Median filter, 5 iterations')

% statistics: 5 iterations
medfilt5_subset = medfilt5_im(100:150,50:100);
medfilt5_mean = mean(medfilt5_subset(:));
medfilt5_stdev = std(medfilt5_subset(:));


%% Alpha-trimmed mean
alpha = 0.25;
nhood = ones(5,5);
n = length(nhood(:));
coeff = 1/(n - (2*floor(alpha*n)));
lower_idx = floor(alpha*n) + 1;
upper_idx = n - floor(alpha*n);
idx_array = linspace(lower_idx, upper_idx, upper_idx - lower_idx + 1); 

% 1 iteration
orderfilt1_im = zeros(size(im));
for i = 1:length(idx_array)
    orderfilt1_im = orderfilt1_im + ordfilt2(im, idx_array(i), nhood);
end
orderfilt1_im = orderfilt1_im*coeff;

figure, imshow(orderfilt1_im)
title('Alpha - trimmed mean: 1 Iteration')

% histogram: 1 iteration
figure, imhist(orderfilt1_im)
title('Histogram: Alpha - trimmed mean, 1 iteration')

% statistics: 1 iteration
orderfilt1_subset = orderfilt1_im(100:150,50:100);
orderfilt1_mean = mean(orderfilt1_subset(:));
orderfilt1_stdev = std(orderfilt1_subset(:));

% 5 iterations
orderfilt5_im = im;
for i = 1:5
    orderfilt1_im = zeros(size(im));
    for j = 1:length(idx_array)
        orderfilt1_im = orderfilt1_im + ordfilt2(orderfilt5_im, idx_array(j), nhood);
    end
    orderfilt5_im = orderfilt1_im*coeff;
end

figure, imshow(orderfilt5_im)
title('Alpha - trimmed mean: 5 Iterations')

% histogram: 5 iterations
figure, imhist(orderfilt5_im)
title('Histogram: Alpha - trimmed mean, 5 iterations')

% statistics: 5 iterations
orderfilt5_subset = orderfilt5_im(100:150,50:100);
orderfilt5_mean = mean(orderfilt5_subset(:));
orderfilt5_stdev = std(orderfilt5_subset(:));


%% Sigma filter
sigma = 20;
[num_rows, num_cols] = size(im);

% 1 iteration
sigma1_im = im;
for row = 3:num_rows-3
    for col = 3:num_cols-3
        nhood = im(row-2:row+2, col-2:col+2);
        difference = abs(nhood - im(row,col));
        vals = find(difference <= 2*sigma);
        im_vals = nhood(vals);
        avg_vals = mean(im_vals(:));
        sigma1_im(row, col) = avg_vals;
    end
end

figure, imshow(sigma1_im)
title('Sigma filter: 1 Iteration')

% histogram: 1 iteration
figure, imhist(sigma1_im)
title('Histogram: Sigma filter, 1 iteration')

% statistics: 1 iteration
sigma1_subset = sigma1_im(100:150,50:100);
sigma1_mean = mean(sigma1_subset(:));
sigma1_stdev = std(sigma1_subset(:));

% 5 iterations
sigma5_im = im;
temp_im = zeros(size(im));
for i = 1:5
    for row = 3:num_rows-3
        for col = 3:num_cols-3
            nhood = sigma5_im(row-2:row+2, col-2:col+2);
            difference = abs(nhood - sigma5_im(row,col));
            vals = find(difference <= 2*sigma);
            im_vals = nhood(vals);
            avg_vals = mean(im_vals(:));
            temp_im(row, col) = avg_vals;
        end
    end
    sigma5_im = temp_im;
end

figure, imshow(sigma5_im)
title('Sigma filter: 5 Iterations')

% histogram: 5 iterations
figure, imhist(sigma5_im)
title('Histogram: Sigma filter, 5 iterations')

% statistics: 5 iterations
sigma5_subset = sigma5_im(100:150,50:100);
sigma5_mean = mean(sigma5_subset(:));
sigma5_stdev = std(sigma5_subset(:));


%% Symmetric nearest-neighbor mean
[num_rows, num_cols] = size(im);

near_neighbor_im1 = zeros(size(im));
for k = 3:num_rows - 3
    for l = 3:num_cols - 3
        nhood = im(k-2:k+2,l-2:l+2);
        center = nhood(3,3);
        min_vals = [];
        for i = 0:2
            for j = 0:2
                pair = [nhood(3-i,3-j), nhood(3+i,3+j)];
                diff = abs([pair - center]);
                [~, min_idx] = min(diff);
                min_vals = cat(1, min_vals, pair(min_idx));
            end
        end
        near_neighbor_im1(k,l) = mean(min_vals);
    end
end


figure, imshow(near_neighbor_im1)
title('Symmetric Nearest-Neighbor Mean: 1 Iteration')

% histogram: 1 iteration
figure, imhist(near_neighbor_im1)
title('Histogram: Symmetric Nearest-Neighbor Mean, 1 iteration')

% statistics: 1 iteration
nearest1_subset = near_neighbor_im1(100:150,50:100);
nearest1_mean = mean(nearest1_subset(:));
nearest1_stdev = std(nearest1_subset(:));                
                
% 5 iterations
near_neighbor_im5 = im;
for z = 1:5
    temp = zeros(size(im));
    for k = 3:num_rows - 3
        for l = 3:num_cols - 3
            nhood = near_neighbor_im5(k-2:k+2,l-2:l+2);
            center = nhood(3,3);
            min_vals = [];
            for i = 0:2
                for j = 0:2
                    pair = [nhood(3-i,3-j), nhood(3+i,3+j)];
                    diff = abs([pair - center]);
                    [~, min_idx] = min(diff);
                    min_vals = cat(1, min_vals, pair(min_idx));
                end
            end
            temp(k,l) = mean(min_vals);
        end
    end
    near_neighbor_im5 = temp;
end

figure, imshow(near_neighbor_im5)
title('Symmetric Nearest-Neighbor Mean: 5 Iterations')

% histogram: 1 iteration
figure, imhist(near_neighbor_im5)
title('Histogram: Symmetric Nearest-Neighbor Mean, 5 iterations')

% statistics: 1 iteration
nearest5_subset = near_neighbor_im5(100:150,50:100);
nearest5_mean = mean(nearest5_subset(:));
nearest5_stdev = std(nearest5_subset(:));



