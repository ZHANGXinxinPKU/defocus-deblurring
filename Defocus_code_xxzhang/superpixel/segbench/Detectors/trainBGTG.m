function trainBGTG(pres)

% get features and labels
[f,y] = sampleDetector(@detector,pres);
f=f'; y=y';
% normalize features to unit variance
fstd = std(f);
fstd = fstd + (fstd==0);
f = f ./ repmat(fstd,size(f,1),1);
% fit the model
fprintf(2,'Fitting model...\n');
beta = logist2(y,f);
% save the result
save(sprintf('beta_bgtg_%s.txt',pres),'fstd','beta','-ascii');

function [f] = detector(im)
[bg,tg] = detBGTG(im);
b = max(bg,[],3);
t = max(tg,[],3);
b = b(:);
t = t(:);
f = [ ones(size(b)) b t ]';
