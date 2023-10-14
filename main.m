clear all;

%Edge Detection
image = double(imread('image-to-rectify.png'))/255;

image = rgb2gray(image);
N = 11;            %Gaussian Kernel Size
sigma = 1;
lowThresh = 0; %Hysterisis Threshhold
highThresh = 0.05;
[edgeIm, grad_angle] = detectEdge(N, sigma, lowThresh, highThresh, image);

%Hough Transform
d = sqrt(size(image,1)^2 + size(image,2)^2); %The diagonal length of polar coordinate
H = zeros(180,ceil(2*d)); 

xyTheta =[]; %which x y coordinate has this Theta
xyR = []; % which x y coordinate has this Rho

%Generate the 2D Hough Histogram
for i=1:size(grad_angle,1)
    for j=1:size(grad_angle,2)
        if (edgeIm(i,j) ~= 0 && grad_angle(i,j) < pi/2 && grad_angle(i,j) > -1*(pi/2))
            angle = grad_angle(i,j);
            if (grad_angle(i,j) <= 0)
                angle = grad_angle(i,j) + (pi);
            end
            theta = ceil(((angle*180)/pi));
            r = sin(angle)*i + cos(angle)*j; %x = j, y = i
            r = ceil(r) + ceil(d);
            
            xyTheta(i,j) = theta; 
            xyR(i,j) = r;
            
            %{
            if (theta == 109 && r == 560)
                disp([i j]);
            end
            %}
            H(theta,r) = H(theta,r) + 1;
        end
    end
end
H = uint8(H); %Hough Transform

%s = max(H,[],2);

localMax = zeros(size(H,1),size(H,2)); %logical matrix showing which location in the 2D Hough has maxima

thresh = round(( 0.20)*max(max(H))); % threshold for filtering out small bins in the 2D Hough histogram

for i=1:size(H,1)
    for j=1:size(H,2)
        if (H(i,j) > thresh)
            localMax(i,j) = H(i,j);
        end
    end
end
localMax = imregionalmax(localMax); %find local maxima locations

localMaxLoc = []; %local max location coord

for i=1:size(localMax,1)
    for j=1:size(localMax,2)
        if (localMax(i,j) == 1)
            %localMaxVector(end+1,1) = H(i,j);
            localMaxLoc(end+1,1) = i;
            localMaxLoc(end,2) = j;
        end
    end
end

%Theta threshold. 2 lines are parallel if the absolute value of the difference
%   of their angles on the 2D hough is less than this. 
%The absolute value of the difference between the angles of
%   2 lines that are not parallel has to be more than this 
theta_thresh = 45;

%Rho threshold. To check that the Rho difference between 2 location on the 2D hough 
%   must be greater than this. 
%If 2 lines have similar angle and Rho then they should be parallel lines in the same cluster
%   of the 2D hough
%If 2 lines have similar angle and different Rho, then they should parallel
%   lines from different clusters of the 2D hough.
r_thresh = d/4; 

len_thresh = 0.2; 

H1= []; %Parallel line pairs horizontal lines, I believe.
H2 = []; %Parallel lines pairs for vertical lines.
%Each row is a pair formatted as following: [Theta_1 Rho_1 Theta_2 Rho_2] 

%Generate the H1 and H2 pairs
for i=1:size(localMaxLoc,1)
    for j=1:size(localMaxLoc,1)
        
        pair = [localMaxLoc(i,1) localMaxLoc(i,2) localMaxLoc(j,1) localMaxLoc(j,2)];

        %isMember = sum(ismember(pair,H1));
        %isMember2 = sum(ismember(pair,H2));
        
        if (isempty(H1))
            isMember = 0;
        else
            isMember = sum(ismember(pair,H1,'rows'));
        end
        
        if (isempty(H2))
            isMember2 = 0;
        else
            isMember2 = sum(ismember(pair,H2,'rows'));
        end
        
        %{
        isValidAngle1 = abs(localMaxLoc(i,1) - localMaxLoc(j,1))  < theta_thresh ...
            && abs(localMaxLoc(i,1) - max(localMaxLoc(:,1))) > theta_thresh ...
            && abs(localMaxLoc(i,1) - min(localMaxLoc(:,1))) > theta_thresh;
        %}
        isValidAngle1 = abs(localMaxLoc(i,1) - localMaxLoc(j,1))  < theta_thresh ...
            && localMaxLoc(i,1)  > ((min(localMaxLoc(:,1)) + max(localMaxLoc(:,1)))/2);
        
        %isValidAngle1 = abs(localMaxLoc(i,1) - localMaxLoc(j,1))  < theta_thresh;
        isValidAngle2 = abs(localMaxLoc(i,1) - localMaxLoc(j,1))  < theta_thresh;
        
        len1 = getPeakLength(xyTheta, xyR, localMaxLoc(i,1), localMaxLoc(i,2));
        len2 = getPeakLength(xyTheta, xyR, localMaxLoc(j,1), localMaxLoc(j,2));
        new_len_thresh = len_thresh*((len1 + len2)/2);
        %abs((len1 - len2)) < new_len_thresh
        
        %if ( i ~= j && isValidAngle1 && abs(localMaxLoc(i,2) - localMaxLoc(j,2)) > r_thresh)
        if ( i ~= j && isValidAngle1 && abs(localMaxLoc(i,2) - localMaxLoc(j,2)) > r_thresh && abs((len1 - len2)) < new_len_thresh)
                if  (isMember == 0 && isMember2 == 0 && localMaxLoc(i,2) < localMaxLoc(j,2))
                    H1(end+1,:) = [localMaxLoc(i,1) localMaxLoc(i,2) localMaxLoc(j,1) localMaxLoc(j,2)];
                elseif (isMember == 0 && isMember2 == 0)
                    H1(end+1,:) = [localMaxLoc(j,1) localMaxLoc(j,2) localMaxLoc(i,1) localMaxLoc(i,2)];
                end
        %elseif ( i ~= j && isValidAngle2 && abs(localMaxLoc(i,2) - localMaxLoc(j,2)) > r_thresh)
        elseif ( i ~= j && isValidAngle2 && abs(localMaxLoc(i,2) - localMaxLoc(j,2)) > r_thresh && abs((len1 - len2)) < new_len_thresh)    
                if  (isMember == 0 && isMember2 == 0 && localMaxLoc(i,2) < localMaxLoc(j,2))
                    H2(end+1,:) = [localMaxLoc(i,1) localMaxLoc(i,2) localMaxLoc(j,1) localMaxLoc(j,2)];
                elseif (isMember == 0 && isMember2 == 0)
                    H2(end+1,:) = [localMaxLoc(j,1) localMaxLoc(j,2) localMaxLoc(i,1) localMaxLoc(i,2)];
                end
        end
    end
end

%------------

%Calculate the areas of each quadrilateral formed by parallel lines
areas = [];
for i=1:size(H1,1)
    for j=1:size(H2,1)       
        
        a1 = cos((H1(i,1)*pi)/180);
        b1 = sin((H1(i,1)*pi)/180);
        c1 = H1(i,2) - d;
        
        a2 = cos((H1(i,3)*pi)/180);
        b2 = sin((H1(i,3)*pi)/180);
        c2 = H1(i,4) - d;

        a3 = cos((H2(j,1)*pi)/180);
        b3 = sin((H2(j,1)*pi)/180);
        c3 = H2(j,2) - d;

        a4 = cos((H2(j,3)*pi)/180);
        b4 = sin((H2(j,3)*pi)/180);
        c4 = H2(j,4) - d;
        
        intersections = [];
        intersections(1,:) = getIntersection([a1,b1,c1],[a3,b3,c3]);
        intersections(2,:) = getIntersection([a1,b1,c1],[a4,b4,c4]);
        intersections(3,:) = getIntersection([a2,b2,c2],[a3,b3,c3]);
        intersections(4,:) = getIntersection([a2,b2,c2],[a4,b4,c4]);
        
        %disp(intersections);
        
        x_axis = intersections(:,1);
        y_axis = intersections(:,2);
        noNegativeX = sum(sum(((x_axis < 1) | (x_axis > size(image,2)) )));
        noNegativeY = sum(sum(((y_axis < 1) | (y_axis > size(image,1)) )));
        if (noNegativeX == 0 && noNegativeY == 0)
            P1 = [intersections(1,1) intersections(2,1) intersections(4,1) intersections(3,1)];
            P2 = [intersections(1,2) intersections(2,2) intersections(4,2) intersections(3,2)];
            polyin = polyshape(P1,P2);
            %figure;plot(polyin)
            %axis equal

            areas(end+1,1) = area(polyin);
            areas(end,2) = i;
            areas(end,3) = j;
        end
        
        %{
        P1 = [intersections(1,1) intersections(2,1) intersections(4,1) intersections(3,1)];
        P2 = [intersections(1,2) intersections(2,2) intersections(4,2) intersections(3,2)];
        polyin = polyshape(P1,P2);
        %figure;plot(polyin)
        %axis equal

        areas(end+1,1) = area(polyin);
        areas(end,2) = i;
        areas(end,3) = j;
        %}
        %disp(area(polyin));
    end
end

%------------

%{
minH1 = [];
minH2 = [];

for i=1:size(H1,1)
    minH1(i) = abs(H1(i,2) - H1(i,4));
end

for i=1:size(H2,1)
    minH2(i) = abs(H2(i,2) - H2(i,4));
end

[M1,ind1] = min(minH1);
[M2,ind2] = min(minH2);
%}
[Min,ind] = max(areas(:,1)); 


ind1 = areas(ind,2); %Index of a parallel line pair of H1
ind2 = areas(ind,3); %Index of a parallel line pair of H2
%The the quadrilateral formed by ind1 and ind2 is the one that has the
%largest area among the lines detected through the thresholds.

figure;imshow(edgeIm); 
%plot the lines on the edge image
hold on

a1 = cos((H1(ind1,1)*pi)/180);
b1 = sin((H1(ind1,1)*pi)/180);
c1 = H1(ind1,2) - d;
x1 = -size(image,2):size(image,2);
y1 = (c1 - a1*x1)/b1;
plot(x1,y1,'red');

a2 = cos((H1(ind1,3)*pi)/180);
b2 = sin((H1(ind1,3)*pi)/180);
c2 = H1(ind1,4) - d;
x2 = -size(image,2):size(image,2);
y2 = (c2 - a2*x2)/b2;
plot(x2,y2,'red');

a3 = cos((H2(ind2,1)*pi)/180);
b3 = sin((H2(ind2,1)*pi)/180);
c3 = H2(ind2,2) - d;
x3 = -size(image,2):size(image,2);
y3 = (c3 - a3*x3)/b3;
plot(x3,y3,'red');

a4 = cos((H2(ind2,3)*pi)/180);
b4 = sin((H2(ind2,3)*pi)/180);
c4 = H2(ind2,4) - d;
x4 = -size(image,2):size(image,2);
y4 = (c4 - a4*x4)/b4;
plot(x4,y4,'red');

hold off

figure;imshow(H,[0,10]);

%plot or superimpose the location of the lines in the 2D hough histogram
hold on

plot(H1(ind1,2), H1(ind1,1),'x', 'MarkerSize', 5, 'color', 'red', 'LineWidth', 1);
plot(H1(ind1,4), H1(ind1,3),'x', 'MarkerSize', 5, 'color', 'red', 'LineWidth', 1);
plot(H2(ind2,2), H2(ind2,1),'x', 'MarkerSize', 5, 'color', 'red', 'LineWidth', 1);
plot(H2(ind2,4), H2(ind2,3),'x', 'MarkerSize', 5, 'color', 'red', 'LineWidth', 1);

hold off

%calculate the intersection points
intersections = [];
intersections(1,:) = getIntersection([a1,b1,c1],[a3,b3,c3]);
intersections(2,:) = getIntersection([a1,b1,c1],[a4,b4,c4]);
intersections(3,:) = getIntersection([a2,b2,c2],[a3,b3,c3]);
intersections(4,:) = getIntersection([a2,b2,c2],[a4,b4,c4]);

aspectRatio = (8.5/11);
newWidth = round(size(image,1) * aspectRatio);
newImageSize = [newWidth size(image,1)];

rectifiedImage = getRectification(image,newImageSize,intersections);
figure;imshow(rectifiedImage);