clc
clear
close all

% convert tif output of tissue analyzer into polygons

imageFolder = "/Users/AndrewTon/Documents/YalePhD/projects/imageProcessing/";
subPath = "woundHealingVideos/Drosophila/Wood_et_al_Nat_Cell_Bio_2002/mov1-WT/healingStack54/";
filename = "handCorrection.tif";
fullPath = imageFolder + subPath + filename;

%first, use ffmpeg -i pcbi_mixed_modes_healing_murrellBanerjee.mp4  healingStack%02d.bmp
%%
numFrames = 1;
writerObj = VideoWriter('segmentedWoundHealingVideo.mp4', 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

for ii=1:numFrames
    clf
    set(0,'DefaultFigureWindowStyle','docked')
    
    %frameii = imgFolder + "healingStack"+sprintf('%02d',ii)+".bmp";
    frameii = fullPath;
    originalImage = imread(frameii);
    
    figure(1)
    imshow(originalImage)
    set(gca, 'YDir', 'normal')
    
    gray1 = im2gray(originalImage);
    
    figure(2)
    [BWmask, maskedImage1] = createMask_lighter(originalImage);
    
    imshow(BWmask)

    %structuring element for imdilate
    se1 = strel('line',6,0);
    se2 = strel('line',6,90);
    %imshow(imdilate(~BWmask,[se1 se2]))
    set(gca, 'YDir', 'normal')
    %BWmask = ~imdilate(~BWmask,[se1 se2]);
    
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
    set(gca, 'YDir', 'reverse')
    
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
    
    vertexPositions = [];
    pixelBuffer = 1;
    for kk = 1 : numberOfBoundaries
	    thisBoundary = boundaries{kk};
	    %plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
        %poly_kk = polyshape(thisBoundary(:,2),thisBoundary(:,1));
        P_reduced = reducepoly(thisBoundary,0.01);
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
    %hold off;
    
    %figure(6)

    %area mask thresholds polygons based on their areas
    sorted_areas = sort(a_array);
    [~,sort_order] = sort(pmeter_array);
    lowBoundAreas = sorted_areas(1);
    upBoundAreas = sorted_areas(end-2);
    areaMask = (a_array > lowBoundAreas) & (a_array < upBoundAreas);

    compositeMask = boxEdgeMask & areaMask;
    
    filteredPerimeter = pmeter_array(compositeMask);
    filteredArea = a_array(compositeMask);
    shapeParameters = filteredPerimeter.^2 ./ filteredArea / (4*pi);
    
    scatter(centroid_array(compositeMask,1), centroid_array(compositeMask,2),'k','x')
   
    currframe = getframe(gcf);
    writeVideo(writerObj,currframe);
end
close(writerObj)

% for each valid segmented cell, record vertex positions separated by nans
% and a header
for kk = 1:numberOfBoundaries
    ii = sort_order(kk);
    thisBoundary = boundaries{ii};
    if (compositeMask(ii) == false)
        continue;
    end
    % fill vertexPositions with a header
    vertexPositions = [vertexPositions; nan nan];
    vertexPositions = [vertexPositions; kk length(thisBoundary(:,1))];
    for jj = 1:length(thisBoundary(:,1))
        %fill vertexPositions with lines of coordinates
        vertexPositions = [vertexPositions; thisBoundary(jj,1) thisBoundary(jj,2)];
    end
end
figure(6)
scatter(vertexPositions(:,2), vertexPositions(:,1),'rx')
daspect([1 1 1])
set(gca, 'Ydir', 'reverse')
