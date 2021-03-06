clear all; close all; clc;

%% C)

f = [0 0 0 0 0;
     0 0 0 0 0;
     5 5 0 0 0;
     5 5 0 0 0;
     5 5 0 0 0]

fx = conv2(f(2:end-1,:), [1 0 -1]/2, 'valid')
fy = conv2(f(:,2:end-1), [1 0 -1]'/2, 'valid')
 
for i=1:3
    for j=1:3
        J0{j,i} = [ fx(j,i)^2,       fx(j,i)*fy(j,i); 
                    fx(j,i)*fy(j,i), fy(j,i)^2 ];
        eigJ0{j,i} = eig(J0{j,i});
    end
end


for i=1:3
    for j=1:3
        J0{j,i} = [ fx(j,i)^2,       fx(j,i)*fy(j,i); 
                    fx(j,i)*fy(j,i), fy(j,i)^2 ];
        eigJ0{j,i} = eig(J0{j,i});
    end
end


%% D)

% before convolution
[EIGV,EIG] = eig(J0{2,2})    


gausk = [1/16 1/8 1/16;
          1/8 1/4 1/8;
          1/16 1/8 1/16]

fx2 = cellfun( @(f) f(1,1), J0 );
fy2 = cellfun( @(f) f(2,2), J0 );
fxy = cellfun( @(f) f(1,2), J0 );

fx2_filt = conv2(gausk,fx2,'same')
fy2_filt = conv2(gausk,fy2,'same')
fxy_filt = conv2(gausk,fxy,'same')

% J0 after convolution with Gaussian kernel
Jrho = [fx2_filt, fxy_filt; 
        fxy_filt, fy2_filt]
[EIGV,EIG] = eig(Jrho)



%% E) F) H)
img = zeros(3,3,2);
img(:,:,1) = [5 5 0; 5 5 0; 5 5 0];
img(:,:,2) = [0 0 5;0 0 5;0 0 5];
img(:,:,3) = zeros(3);

fx = imfilter(img,[-1 0 1]/2);
fy = imfilter(img,[-1 0 1]'/2);

dfR = [fx(2,2,1),fy(2,2,1)]'
dfG = [fx(2,2,2),fy(2,2,2)]'
dfB = [fx(2,2,3),fy(2,2,3)]'

% E)
norm( dfR + dfG + dfB )
% F)
norm( [norm(dfR) norm(dfG) norm(dfB) ] )
% H)
[EIGV,EIG] = eig( dfR*dfR' + dfG*dfG' + dfB*dfB' )