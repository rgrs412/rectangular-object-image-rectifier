function I = getRectification(grayIm,newImageSize,intersections)

image = uint8(grayIm*255);
%figure;imshow(image);
%{
A = [1 1 1 0 0 0 0 0 0;
    0 0 0 1 1 1 0 0 0;
    343 1 1 0 0 0 0 0 0;
    0 0 0 343 1 1 0 0 0;
    1 471 1 0 0 0 0 0 0;
    0 0 0 1 471 1 0 0 0;
    343 471 1 0 0 0 0 0 0;
    0 0 0 343 471 1 0 0 0];

B = [0 0 0 0 0 0 108*1 108*1 108;
    0 0 0 0 0 0 70*1 70*1 70;
    0 0 0 0 0 0 316*343 316*1 316;
    0 0 0 0 0 0 78*343 78*1 78;
    0 0 0 0 0 0 18*1 18*471 18;
    0 0 0 0 0 0 326*1 326*471 326;
    0 0 0 0 0 0 308*343 308*471 308;
    0 0 0 0 0 0 373*343 373*471 373];
%}
x = newImageSize(1);
y = newImageSize(2);

A = [1 1 1 0 0 0 0 0 0;
    0 0 0 1 1 1 0 0 0;
    x 1 1 0 0 0 0 0 0;
    0 0 0 x 1 1 0 0 0;
    1 y 1 0 0 0 0 0 0;
    0 0 0 1 y 1 0 0 0;
    x y 1 0 0 0 0 0 0;
    0 0 0 x y 1 0 0 0];

x1 = intersections(1,1);
y1 = intersections(1,2);
x2 = intersections(2,1);
y2 = intersections(2,2);
x3 = intersections(3,1);
y3 = intersections(3,2);
x4 = intersections(4,1);
y4 = intersections(4,2);

B = [0 0 0 0 0 0 x1*1 x1*1 x1;
    0 0 0 0 0 0 y1*1 y1*1 y1;
    0 0 0 0 0 0 x2*x x2*1 x2;
    0 0 0 0 0 0 y2*x y2*1 y2;
    0 0 0 0 0 0 x3*1 x3*y x3;
    0 0 0 0 0 0 y3*1 y3*y y3;
    0 0 0 0 0 0 x4*x x4*y x4;
    0 0 0 0 0 0 y4*x y4*y y4];

M = ((A.')*A) + ((B.')*B) - ((B.')*A) - ((A.')*B);
[U,S,V] = svd(M);
h = reshape(V(:,9), 3, 3);
h = h.';

I = zeros(y,x);
a= [];
b= [];
for i=1:size(I,1)
    for j=1:size(I,2)
        w_p = (h(3,1)*j + h(3,2)*i + h(3,3));
        x_p = (h(1,1)*j + h(1,2)*i + h(1,3))/w_p;
        y_p = (h(2,1)*j + h(2,2)*i + h(2,3))/w_p;
        
        %x = (h(1,1)*x_p + h(1,2)*y_p + h(1,3))/w_p;
        %y = (h(2,1)*x_p + h(2,2)*y_p + h(2,3))/w_p;
        a(end+1,:) = x_p;
        b(end+1,:) = y_p;

        I(i,j) = image(round(y_p),round(x_p));
        %I(ceil(abs(x_p)),ceil(abs(y_p))) = image(i,j);
    end
end
I = uint8(I);
%figure;imshow(I);