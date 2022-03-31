clc
clear
close all
%first, use ffmpeg -i pcbi_mixed_modes_healing_murrellBanerjee.mp4  healingStack%02d.bmp
%%
numFrames = 92;
showShapeHistogram = 0;
healingTypeStr = "TetleyStaddonWT";
writerObj = VideoWriter(healingTypeStr+".mp4", 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

for ii=1:numFrames
    clf
    imgFolder = "videosForImageJ/";
    set(0,'DefaultFigureWindowStyle','docked')
    
    frameii = imgFolder + healingTypeStr+"_healing_"+sprintf('%02d',ii)+".bmp";
    originalImage = imread(frameii);
    
    %figure(1)
    %imshow(originalImage)
    %set(gca, 'YDir', 'normal')
    
    gray1 = im2gray(originalImage);
    
    % mask for murrellBanerjee
    %[BWmask, maskedImage1] = createMask_lighter(originalImage);

    %figure(2)

    % mask for tetley-staddon
    [BWmask, maskedImage1] = createMaskTetleyStaddon(originalImage);
    %BWmask = ~BWmask;
    %imshow(BWmask)
    %set(gca, 'YDir', 'normal')

    %structuring element for imdilate (use 6,0 and 6,90 for
    %murrellBanerjee)
    se1 = strel('line',2,0);
    se2 = strel('line',2,90);
    %imshow(imdilate(~BWmask,[se1 se2]))
    %set(gca, 'YDir', 'normal')
    BWmask = ~imdilate(~BWmask,[se1 se2]);
    
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
    
    pixelBuffer = 10;
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
    lowBoundAreas = 0.05*1e4;
    upBoundAreas = 1e4;
    areaMask = (a_array > lowBoundAreas) & (a_array < upBoundAreas);

    compositeMask = boxEdgeMask & areaMask;
    
    %scatter(pmeter_array(areaMask),a_array(areaMask),'filled','k')
    %ylabel('polygon area', 'fontsize',24)
    %xlabel('polygon perimeter', 'fontsize',24)
    %title('areaMask test')
    
    %figure(7)
    %scatter(pmeter_array(boxEdgeMask),a_array(boxEdgeMask),'filled','k')
    %ylabel('polygon area', 'fontsize',24)
    %xlabel('polygon perimeter', 'fontsize',24)
    %title('boxEdgeMask test')
    
    %figure(8)
    %scatter(pmeter_array(compositeMask),a_array(compositeMask),'filled','k')
    %ylabel('polygon area', 'fontsize',24)
    %xlabel('polygon perimeter', 'fontsize',24)
    %title('compositeMask test')
    
    filteredPerimeter = pmeter_array(compositeMask);
    filteredArea = a_array(compositeMask);
    shapeParameters = filteredPerimeter.^2 ./ filteredArea / (4*pi);
    
    scatter(centroid_array(compositeMask,1), centroid_array(compositeMask,2),'k','x')
    
    %figure(9)
    %histEdges = 1.0:0.04:3.0;
    %histogram(shapeParameters,histEdges);
    %xticks([1.0:0.4:3.0])
    %xlabel('$\mathcal{A}$','Interpreter','LaTeX','Fontsize',16)
    %ylabel('$\mathcal{A}$','Interpreter','LaTeX','Fontsize',16)
    %ylim([0 10])

    if (showShapeHistogram)
        figure(5)
        histAx = axes('Position', [.7 .7   .25 .25]);
        box on
        histEdges = 1.0:0.05:4.0;
        histogram(shapeParameters,histEdges);
        xticks(1.0:0.5:4.0);
        xlabel('$\mathcal{A}$','Interpreter','LaTeX','Fontsize',16)
        ylabel('$\mathcal{A}$','Interpreter','LaTeX','Fontsize',16)
        ylim([0 20])
        annotation('textbox',[.75 .8 .1 .1], 'edgecolor', 'none', 'string', "mean="+num2str(mean(shapeParameters)))
        annotation('textbox',[.75 .77 .1 .1], 'edgecolor', 'none', 'string', "std="+num2str(std(shapeParameters)))
    end
    currframe = getframe(gcf);
    writeVideo(writerObj,currframe);
end
close(writerObj)
