clc
clear
close all
set(0,'DefaultFigureWindowStyle','docked')
% convert tif output of tissue analyzer into polygons

%parameter list
reducePolyTolerance = [0.01,0.02]; 
%tolerance for reducePoly algorithm

imageFolder = "/Users/AndrewTon/Documents/YalePhD/projects/imageProcessing/";
%subPath = "woundHealingVideos/Drosophila/Wood_et_al_Nat_Cell_Bio_2002/mov1-WT/healingStack54/";
%subPath = "woundHealingVideos/testImagesForShapeCalculation/honeycomb_color_inverted/";
subPath = "woundHealingVideos/Drosophila/Tetley_et_al_Nat_Phys_2019/embryoHealedColored/";
filename = "handCorrection.tif";
fullPath = imageFolder + subPath + filename;

frameii = fullPath;
originalImage = imread(frameii);

figure(1)
imshow(originalImage)
set(gca, 'YDir', 'normal')
gray1 = im2gray(originalImage);

figure(2)
[BWmask, maskedImage1] = createMask_lighter(originalImage);
imshow(BWmask)

%% BWmask is the base image we're going to work with using blob counting
% get labels for connected objects in BWmask (polygons)
%figure(3)
labeledImage = bwlabel(BWmask);
%imshow(labeledImage, []);  % Show the gray scale image.
%set(gca, 'YDir', 'normal')

% Let's assign each blob a different color to visually show the user the distinct blobs.
%figure(4)
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');
%imshow(coloredLabels);
%set(gca, 'YDir', 'normal')
% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(labeledImage, BWmask, 'all');
numberOfBlobs = size(blobMeasurements, 1);

% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
figure(5)
%imshow(BWmask);

hold on;
% store box dimensions for boxEdgeMask, which throws out polygons too close
% to the edge
boxDims = [length(BWmask(1,:)) length(BWmask(:,1))];

boundaries = bwboundaries(BWmask);
numberOfBoundaries = size(boundaries, 1);
a_array = zeros(numberOfBoundaries,1);
pmeter_array = zeros(numberOfBoundaries,1);
centroid_array = zeros(numberOfBoundaries,2);
boxEdgeMask = zeros(numberOfBoundaries,1,'logical');


shapeArraysTolerance = zeros(length(reducePolyTolerance), 20);
% looping over tolerances, record shapes to get standard error
for ii = 1:length(reducePolyTolerance)
    vertexPositions = [];
    pixelBuffer = 1;
    for kk = 1 : numberOfBoundaries
        thisBoundary = boundaries{kk};
    
        P_reduced = reducepoly(thisBoundary,reducePolyTolerance(ii));
        %reducepoly to reduce vertex density and do some smoothing
    
        plot(P_reduced(:,2), P_reduced(:,1), 'r', 'linewidth', 1);
        poly_kk = polyshape(P_reduced(:,2),P_reduced(:,1));
        a_array(kk) = area(poly_kk);
        pmeter_array(kk) = perimeter(poly_kk);
        vertices = poly_kk.Vertices;
        centroid_array(kk,:) = mean(vertices);
        xMask = (vertices(:,1) < boxDims(1) - pixelBuffer) & (vertices(:,1) > pixelBuffer);
        yMask = (vertices(:,2) < boxDims(2) - pixelBuffer) & (vertices(:,2) > pixelBuffer);
        boxEdgeMask(kk) = all(xMask) & all(yMask);
    end
    
    %area mask thresholds polygons based on their areas
    sorted_areas = sort(a_array);
    [~,sort_order] = sort(pmeter_array);
    lowBoundAreas = sorted_areas(1);
    upBoundAreas = sorted_areas(end-2);
    %areaMask = (a_array > lowBoundAreas) & (a_array < upBoundAreas);
    areaMask = (a_array-median(a_array)<median(a_array));
    
    compositeMask = boxEdgeMask & areaMask;
    
    filteredPerimeter = pmeter_array(compositeMask);
    filteredArea = a_array(compositeMask);
    shapeParameters = filteredPerimeter.^2 ./ filteredArea / (4*pi);
    mean(shapeParameters)

    %save shapeParameters at this tolerance
    shapeArraysTolerance(ii,1:length(shapeParameters)) = shapeParameters;
    
    scatter(centroid_array(compositeMask,1), centroid_array(compositeMask,2),'k','x')
    
    % for each valid segmented cell, record vertex positions separated by nans
    % and a header, sorted increasing by cell perimeters
    figure(6)
    clf; hold on;
    set(gca, 'YDir', 'reverse')
    for kk = 1:numberOfBoundaries
        jj = sort_order(kk);
        thisBoundary = boundaries{jj};
        if (compositeMask(jj) == false)
            continue;
        end
        % fill vertexPositions with a header
        vertexPositions = [vertexPositions; nan nan];
        vertexPositions = [vertexPositions; kk length(thisBoundary(:,1))];
        for ll = 1:length(thisBoundary(:,1))
            %fill vertexPositions with lines of coordinates
            vertexPositions = [vertexPositions; thisBoundary(ll,1) thisBoundary(ll,2)];
        end
        
        % redo old plot to show new segmentation
        P_reduced = reducepoly(thisBoundary,reducePolyTolerance(ii));
        plot(P_reduced(:,2), P_reduced(:,1), 'r', 'linewidth', 1);
        daspect([1 1 1])
    end
    figure(7)
    scatter(vertexPositions(:,2), vertexPositions(:,1),'rx')
    daspect([1 1 1])
end

diffShapeArrays = diff(shapeArraysTolerance);
diffShapeArrays(diffShapeArrays == 0) = [];
std_error = sqrt(mean(diffShapeArrays.^2))
