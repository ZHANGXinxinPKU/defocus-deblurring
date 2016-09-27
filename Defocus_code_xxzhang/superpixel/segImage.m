function I_s = segImage(I,S)
[cx,cy] = gradient(S);
ccc = (abs(cx)+abs(cy))~=0;
I_s = I;
I_s(:,:,1) = max(I_s(:,:,1),ccc);
I_s(:,:,2) = min(I_s(:,:,2),~ccc);
I_s(:,:,3) = min(I_s(:,:,3),~ccc);
