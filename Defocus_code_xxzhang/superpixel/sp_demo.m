% sp_demo.m
%
% See instructions in README.

if ~exist('cncut')
    addpath('yu_imncut');
end

J = im2double(imread('pumpkin.png'));
I(:,:,1) = J;
I(:,:,2) = J;
I(:,:,3) = J;
N = size(I,1);
M = size(I,2);

% Number of superpixels coarse/fine.
N_sp=200;
N_sp2=1000;
% Number of eigenvectors.
N_ev=50;


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
% Sp2 = clusterLocations(Sp,ceil(N*M/N_sp2));
% fprintf(' took %.2f minutes\n',etime(clock,st)/60);


I_sp = segImage(I,Sp);
% I_sp2 = segImage(I,Sp2);
I_seg = segImage(I,Seg);
figure;
subplot(1,4,1);
imshow(I);
subplot(1,4,2);
imshow(I_seg);
subplot(1,4,3);
imshow(I_sp);
% subplot(1,4,4);
% imshow(I_sp2);


