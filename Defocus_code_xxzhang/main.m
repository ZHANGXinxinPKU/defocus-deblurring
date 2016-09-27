
% Zhang X, Wang R, Jiang X, et al.
% Spatially variant defocus blur map estimation and deblurring from a single image.
% Journal of Visual Communication and Image Representation, 2016, 35: 257-264

clc;
clear all;
close all;

addpath(genpath('superpixel'));
addpath(genpath('BM3D'));
% addpath(genpath('vlfeat-0.9.20-bin'));

%% parameter setting
% they can be adjusted by users
number = 5; % the number of superpixels
std=3.5;  %3
lambda=0.05;   %0.1; 0.3
maxBlur=3;
sigma = 0.005;  %0.01 BM3D restoration  The larger the sigma is, the restored image becomes more smooth.
prescale = 1;

%%
fn = 'image\out_of_focus0474.jpg';
kernel_size = 15;
y = imread(fn);
y = imresize(y,[size(y,1)*prescale,size(y,2)*prescale]);
[hei wid dimen] = size(y);

isselect = 0; % deblur an image or a patch
if isselect ==1
    figure, imshow(y);
    fprintf('Please choose an area for deblurring \n');
    h = imrect;
    position = wait(h);
    close;
    B_patch = imcrop(y,position);
    y = (B_patch);
end
if size(y,3)==3
    yg = im2double(rgb2gray(y));
else
    yg = im2double(y);
end
y = im2double(y);

%% defocus map estimation

eth=0.0051; % thershold for canny edge detector   0.1
edgeMap=edge(rgb2gray(y),'canny',eth,1);
%estimate the defocus map
[sparsemap, fullmap, sharp] = defocus(y,edgeMap,std,lambda,maxBlur);

figure,
imshow(sparsemap);title('sparse_map')

fullmap(fullmap<0)=0;
figure,
imshow(fullmap/max(fullmap(:)));title('smooth_fullmap')

J = fullmap/max(fullmap(:));


%% superpixel
% from Greg Mori's work.  "Guiding Model Search Using Segmentation"
% Greg Mori's method is slow, it can be replaced by SLIC superpixel method

% [labels, numlabels] = slicmex(im2uint8(J),100,0); % SLIC    faster

if ~exist('cncut')
    addpath('yu_imncut');
end
I(:,:,1) = J;
I(:,:,2) = J;
I(:,:,3) = J;
N = size(I,1);
M = size(I,2);

% Number of superpixels.
N_sp=200;
N_sp2=1000;
% Number of eigenvectors.
N_ev=number;   % pumpkinN_ev=30;   % faceN_ev=10;  % can be adjusted by users


% ncut parameters for superpixel computation
diag_length = sqrt(N*N + M*M);
par = imncut_sp;
par.int=0;
par.pb_ic=1;
par.sig_pb_ic=0.05;
par.sig_p=ceil(diag_length/50);
par.verbose=0;
par.nb_r=ceil(diag_length/60);
par.rep = -0.005;  % stability?  or proximity?
par.sample_rate=0.2;
par.nv = N_ev;
par.sp = N_sp;

% Intervening contour using mfm-pb
fprintf('running PB\n');
[emag,ephase] = pbWrapper(I,par.pb_timing);
emag = pbThicken(emag);
par.pb_emag = emag;
par.pb_ephase = ephase;
clear emag ephase;

st=clock;
fprintf('Ncutting...');
[Sp,Seg] = imncut_sp(I,par);
fprintf(' took %.2f minutes\n',etime(clock,st)/60);

st=clock;
fprintf('Fine scale superpixel computation...');

I_seg = segImage(I,Seg);
figure,
imshow(I_seg);
% imshow(labels,[0 numlabels]); title('superpixel');
labels = Seg;

%% deconvolution and remove ringing
ave = zeros(N_ev,1);
for la = 1:N_ev
    map = zeros(hei,wid);
    i = 0;
    for n = 1:hei
        for m = 1:wid
            if labels(n,m)==la
                map(n,m) = fullmap(n,m);
                i = i+1;
            end
        end
    end
    ave(la,1) = sum(map(:))/i;   % average value
    
    [qq ww ee]=find(map(:));
    sig(la,1) = min(ee);   % minimum value
    for n = 1:hei
        for m = 1:wid
            if labels(n,m)==la
                ffullmap(n,m) = ave(la,1);   % average value as the sigma
                fffullmap(n,m) = sig(la,1);   % minimum value as the sigma
            end
        end
    end
    
end

figure,
imshow(ffullmap,[0 max(ffullmap(:))]);title('superpixel_fDmap');  % sigma
figure,
imshow(fffullmap,[0 max(fffullmap(:))]);title('superpixel_fffDmap');  % sigma

num = N_ev;
min_sig = min(fullmap(:));
max_sig = max(fullmap(:));

b = zeros(kernel_size);
bhs = floor(size(b, 1)/2);
for j = 1:size(y,3)
    pad(:,:,j) = padarray(y(:,:,j), [1 1] * bhs, 'replicate', 'both');
end

Latent = zeros(hei,wid,3);

%% refocus in each superpixel
for j = 1:num
    Latent_1{j} = zeros(hei,wid,3);
end
%% multi-scale sigma

for j = 1:num   % each superpixel
    Latent_1{j} = zeros(hei,wid,3);
    blur = zeros(hei,wid,3);
    del_img1 = zeros(hei,wid);
    kernel = fspecial('gaussian',kernel_size,ave(j,1));  %  the average value of all pixels in a superpixel as the blur amount of the superpixel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BM3D
    
    for c = 1:size(pad,3)
        Latent2(:,:,c) = BM3DDEB(pad(:,:,c),kernel,sigma);
    end
    Latent1 = Latent2(bhs + 1 : end - bhs, bhs + 1 : end - bhs, :);
    
    % function (18)
    re_blur = conv2(rgb2gray(Latent1),kernel,'same');
    del_img = re_blur-rgb2gray(y);   % the difference of the current superpixel
    
    for n = 1:hei
        for m = 1:wid
            if labels(n,m)==j
                Latent_1{j}(n,m,:) = Latent1(n,m,:);     % the latent image of the current superpixel
                del_img1(n,m) = del_img(n,m);
            end
        end
    end
    %     figure,imshow(del_img1);
    
    del_abs = abs(del_img1);
    
    [hh ll] = find(abs(del_img1)>0.0001);
    [ch ku] = size(hh);
    kernel1 = fspecial('gaussian',kernel_size,sig(j,1));
    for c = 1:size(pad,3)
        Latent3(:,:,c) = BM3DDEB(pad(:,:,c),kernel1,sigma);
    end
    Latent4 = Latent3(bhs + 1 : end - bhs, bhs + 1 : end - bhs, :);
    for c = 1:ch
        Latent_1{j}(hh(c),ll(c),:) = Latent4(hh(c),ll(c),:);   % remove ringing  % Function£¨19£©in the paper
    end
    
    Latent = Latent+Latent_1{j};
    figure,imshow(Latent_1{j});
end

figure,imshow(Latent); title('Latent image');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write to file
imwrite(J,['result_set\' fn(7:end-4) '_JDmap4.png']);
imwrite(sparsemap,['result_set\' fn(7:end-4) '_sparsemap4.png']);
imwrite(ffullmap/max(ffullmap(:)),['result_set\' fn(7:end-4) '_ffullmap4.png']);
imwrite(Latent,['result_set\' fn(7:end-4) '_result4.png']);
imwrite(I_seg,['result_set\' fn(7:end-4) '_I_seg4.png']);
imwrite(fffullmap/max(fffullmap(:)),['result_set\' fn(7:end-4) '_fffullmap4.png']);

