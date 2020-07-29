X = imread('lake.gif');
i_1 = inputdlg('number of iterations');
k_1 = inputdlg('number of K');
g_1 = inputdlg('g():1 for exponential  2 for inverse quadratic') 
i = str2num(cell2mat(i_1));
k = str2num(cell2mat(k_1));
g = str2num(cell2mat(g_1));
[M,N] = size(X);
iterateImage = double(X);
FinalImage = zeros(M,N);1
iterateNum = i;
K = k;
g_mode = g;

while(iterateNum > 0)
    
    for x = 1 : M 
        for y = 1 : N
            if(x ~=1 &&x~=M &&y ~=1&&y ~=N)
            GradientN = iterateImage(x,y+1) - iterateImage(x,y);
            GradientS = iterateImage(x,y-1) - iterateImage(x,y);
            GradientE = iterateImage(x+1,y) - iterateImage(x,y);
            GradientW = iterateImage(x-1,y) - iterateImage(x,y);
            Cn = ConCoefficient(g_mode,GradientN,K);
            Cs = ConCoefficient(g_mode,GradientS,K);
            Ce = ConCoefficient(g_mode,GradientE,K);
            Cw = ConCoefficient(g_mode,GradientW,K);
            FinalImage(x,y) = iterateImage(x,y) +0.25*(GradientN * Cn + GradientS * Cs + GradientE * Ce + GradientW * Cw);
            else
                FinalImage(x,y) = iterateImage(x,y);
            end
            
        
        end
    end
    iterateNum = iterateNum -1;
    iterateImage = FinalImage;
end
gray128 = zeros(1,M)
for x = 1 : M 
    gray128(1,x) = iterateImage(x,128);
end
I = uint8(iterateImage) 
SpokeImage = iterateImage
for x = 1 : M 
       for y = 1 : N 
           if SpokeImage(x,y) < 89 || SpokeImage(x,y)>105
               SpokeImage(x,y) = 255;
           end
       end
end       
figure(1)
imshow(X)
title('original image')


figure(2)
imshow(I)
title('new image')

           
           
           