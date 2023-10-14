function [binaryIm, grad_angle] = detectEdge(N, sigma, lowThresh, highThresh, image)

I = image;

%Gaussian Kernel Calculation
W = zeros(N,N);
for i=1:size(W,1)
    for j=1:size(W,2)
        W(i,j) = exp(-(((i-((N+1)/2))^2) + ((j-((N+1)/2))^2))/2*sigma^2);
    end
end

%Discretizing and Normalization
Min = min(W);
Min = min(Min);
W_disc = ceil(W/Min); %discretize W kernel
W_norm = W_disc/(sum(sum(W_disc)));

gaussian_smoothing = conv2(I, W_norm,'same');
I = gaussian_smoothing;
%gaussian_smoothing = uint8(gaussian_smoothing*255);
%figure;imshow(gaussian_smoothing);

%gradient kernels
dx = [0 0 0; -1/2 0 1/2; 0 0 0];
dx = flip(dx, 2);
dy = [0 -1/2 0; 0 0 0; 0 1/2 0];
dy = flip(dy, 1);

dW_dx = conv2(W_norm, dx, 'same');
dW_dy = conv2(W_norm, dy, 'same');
dI_dx = conv2(I, dW_dx, 'valid');
dI_dy = conv2(I, dW_dy, 'valid');
Ixy(:,:) = sqrt(dI_dx.^2 + dI_dy.^2); %magnitude
%Ixy = dI_dx + dI_dy;
Ixy = uint8(Ixy*255);

%gradient_magnitude_image = Ixy;
%figure;imshow(gradient_magnitude_image);

%Hysterisis
low = lowThresh*255;
high = highThresh*255;
binaryIm = zeros(size(Ixy));

for i=1:size(Ixy,1)
    for j=1:size(Ixy,2)
        if ((Ixy(i,j) > low))
            
            %Check neighbors
            neighbors = [i-1 j-1; i-1 j; i-1 j+1; i j-1; i j+1; i+1 j-1; i+1 j; i+1 j+1];
            for r=1:size(neighbors,1)
                if ((neighbors(r,1) > 0) && ((neighbors(r,1) <= size(Ixy,1))) && ...
                    (neighbors(r,2) > 0) && ((neighbors(r,2) <= size(Ixy,2))))
                    
                    if (Ixy(neighbors(r,1),neighbors(r,2)) > high)
                        binaryIm(i,j) = 1;
                    end        
                end
            end
            
        elseif (Ixy(i,j) > high)
            binaryIm(i,j) = 1;
        end
    end
end

%histerisis = uint8(binaryIm*255);
%figure;imshow(histerisis);

grad_angle = zeros(size(binaryIm)); %gradient angle matrix
for i=1:size(grad_angle,1)
    for j=1:size(grad_angle,2)
        if (binaryIm(i,j) == 1)
            grad_angle(i,j) = atan(dI_dy(i,j)/dI_dx(i,j));
            %grad_angle(i,j) = atan2(dI_dy(i,j),dI_dx(i,j));
        end
    end
end

%{
%non-max suppression calculations
im = Ixy;
for i=1:size(grad_angle,1)
    for j=1:size(grad_angle,2)
        if (grad_angle(i,j) ~= 0)
            if ((grad_angle(i,j) >= (7*pi)/8) && (grad_angle(i,j) < pi/8))
                if ((j-1 > 0) && (j+1 <= size(im,2)))
                    left = im(i,j-1);
                    right = im(i,j+1);
                    if ((im(i,j)> left) && (im(i,j) > right))
                        im(i,j) = 255;
                    end
                end
            elseif ((grad_angle(i,j) <= pi/8) && (grad_angle(i,j) < (3*pi)/8))
                if ((i-1 > 0) && (i+1 < size(im,1)) && ...
                        (j-1 > 0) && (j+1 <= size(im,2)))
                    up_right = im(i-1,j+1);
                    down_left = im(i+1,j-1);
                    if ((im(i,j)> up_right) && (im(i,j) > down_left))
                        im(i,j) = 255;
                    end
                end
            elseif ((grad_angle(i,j) <= (3*pi)/8) && (grad_angle(i,j) < (5*pi)/8))
                if ((i-1 > 0) && (i+1 <= size(im,1)))
                    up = im(i-1,j);
                    down = im(i+1,j);
                    if ((im(i,j)> up) && (im(i,j) > down))
                        im(i,j) = 255;
                    end
                end
            elseif ((grad_angle(i,j) <= (5*pi)/8) && (grad_angle(i,j) < (7*pi)/8))
                if ((i-1 > 0) && (i+1 < size(im,1)) && ...
                        (j-1 > 0) && (j+1 <= size(im,2)))
                    up_left = im(i-1,j-1);
                    down_right = im(i+1,j+1);
                    if ((im(i,j)> up_left) && (im(i,j) > down_right))
                        im(i,j) = 255;
                    end
                end
            end
        end
    end
end

%figure;imshow(im);
%}

