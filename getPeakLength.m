function len = getPeakLength(xyTheta, xyR, theta, r)

coords = [];

for i=1:size(xyTheta,1)
    for j=1:size(xyTheta,2)
        if (xyTheta(i,j) == theta && xyR(i,j) == r)
            coords(end+1,:) = [i j];
        end
    end
end

sortedCoords = sortrows(coords);
x1 = sortedCoords(1,1);
y1 = sortedCoords(1,2);
x2 = sortedCoords(end,1);
y2 = sortedCoords(end,2);

len = sqrt((x2 - x1)^2 + (y2 - y1)^2);