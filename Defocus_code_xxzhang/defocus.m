function [sparseDMap, fullDMap, sharpness] = defocus(I,edgeMap,std,lambda,maxBlur)

[hei wid dimen] = size(I);
dx = [1 -1];
dy = [1 -1]';

for hh=1:hei 
    for ww=1:wid
        gI(hh,ww) = max(I(hh,ww,:));
    end
end
w1 = (2*ceil(2* std))+1;
g1 = fspecial('gaussian',w1,std);
I1 = conv2(gI,g1,'same');

I1x = conv2(I1,dx,'same');
I1y = conv2(I1,dy,'same');
Igx = conv2(gI,dx,'same');
Igy = conv2(gI,dy,'same');

Igg = sqrt(Igx.^2 + Igy.^2);
I11 = sqrt(I1x.^2 + I1y.^2);

sharpness = ((Igg-I11)./(Igg+10^(-6))); % Function (4)
% figure,imshow(sharpness);

% Generate the sparse defocus map 
idx=find(edgeMap~=0);
sparseDMap=zeros(size(gI));

% 2015/3/22
for jj=1:length(idx)
    if sharpness(idx(jj))<2&sharpness(idx(jj))>0
        sparseDMap(idx(jj))=sqrt(((1-sharpness(idx(jj)))^2*std^2)/(2*sharpness(idx(jj))-sharpness(idx(jj)).^2)); % Function (8)
    end
end

% remove noise
sparseDMap(sparseDMap>maxBlur)=maxBlur;
[sparseDMap]= guidedfilter(edgeMap,sparseDMap,5,1e-9);
sparseDMap(sparseDMap<1e-8)=0;
%% Generate the full defocus map
im = I;
image = sparseDMap;
level=30;    %5
[m n d]=size(im);
nn=[10;1];    %nn=[10;1];
[a b]=ind2sub([m n],1:m*n);
feature=[reshape(im,m*n,d)';[a;b]/sqrt(m*m+n*n)*level+rand(2,m*n)*1e-6];
% size(feature)
now=0;
for i=1:size(nn,1)
    ind=vl_kdtreequery(vl_kdtreebuild(feature),feature,feature,'NUMNEIGHBORS',nn(i),'MAXNUMCOMPARISONS',nn(i)*2);
    a=reshape(repmat(uint32(1:m*n),nn(i),1),[],1);
    b=reshape(ind,[],1);
    row(now+1:now+m*n*nn(i),:)=[min(a,b) max(a,b)];
    feature(d+1:d+2,:)=feature(d+1:d+2,:)/100;
    now=now+m*n*nn(i);
end

value=max(1-sum(abs(feature(1:d+2,row(:,1))-feature(1:d+2,row(:,2))))/(d+2),0);%1/(norm(x1-x2,1)+0.1);
A=sparse(double(row(:,1)),double(row(:,2)),value,m*n,m*n);
A=A+A';
D=spdiags(sum(A,2),0,n*m,n*m);
num = 1;
map = image>0.00001;
M=D-A+lambda*spdiags(map(:),0,m*n,m*n);
L=ichol(M);
i=1;
tot = 0;
while i<=num
    val = image(:);
    i=i+1;
    tot=tot+1;
    alpha(:,tot)=pcg(M,lambda*val(:),1e-10,2000,L,L');
    fullDMap=reshape(alpha(:,tot),m,n);
%     imshow(fullDMap/max(fullDMap(:)));
end


