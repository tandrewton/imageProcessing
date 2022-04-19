clc
clear
close all
set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultAxesFontSize', 24);
% convert tif output of tissue analyzer into polygons

%parameter list
% 0 tolerance means no transformation is applied
reducePolyTolerance = linspace(0,0.1,50); 
makeAMovie = true;
%tolerance for reducePoly algorithm

imageFolder = "/Users/AndrewTon/Documents/YalePhD/projects/imageProcessing/";
outputFolder = imageFolder+"segmentationShapeCalculationTests/";
subPath_Wood_1 = "woundHealingVideos/Drosophila/Wood_et_al_Nat_Cell_Bio_2002/mov1-WT/healingStack54/";
subPath_honeycomb = "woundHealingVideos/testImagesForShapeCalculation/honeycomb_color_inverted/";
subPath_Tetley = "woundHealingVideos/Drosophila/Tetley_et_al_Nat_Phys_2019/embryoHealedColored/";
subPath_Campos = "woundHealingVideos/Drosophila/Campos_et_al_Genetics_2010/Fig3_cropped_final/";
subPath_Garcia_Fernandez = "woundHealingVideos/Drosophila/Garcia_Fernandez_et_al_Int_Dev_Bio_2009/Fig1_final/";
subPath_Abreu_Blanco = "woundHealingVideos/Drosophila/Abreu_Blanco_et_al_Cell_Science_2012/healingStack42_cropped/";
subPath_Stramer = "woundHealingVideos/Drosophila/Stramer_et_al_Cell_Bio_2005/Stramer_Fig_1A_cropped_PreProcessed/";
subPath_Verboon = "woundHealingVideos/Drosophila/Verboon_et_al_Small_GTPases_2015/Verboon_Fig1A_cropped_healed-01/";
subPath_Wood_2 = "woundHealingVideos/Drosophila/Wood_et_al_Nat_Cell_Bio_2002/WoodAdaptedFromParkhurstFigure_healed-01/";

allWoundedPaths = [subPath_Wood_1, subPath_Tetley, subPath_Campos, ...
    subPath_Garcia_Fernandez, subPath_Abreu_Blanco, subPath_Stramer,...
    subPath_Verboon, subPath_Wood_2];
%allWoundedPaths = [subPath_Verboon];

subPath_unwounded_Abreu_Blanco = "woundHealingVideos/Drosophila/Abreu_Blanco_et_al_Cell_Science_2012/healingStack01_cropped/";
subPath_unwounded_Verboon = "woundHealingVideos/Drosophila/Verboon_et_al_Small_GTPases_2015/Verboon_Fig_1A_unwounded_cropped-01/";
%subPath_unwounded_Fernandez_Gonzalez = "woundHealingVideos/Drosophila/Fernandez_Gonzalez_et_al_MBOC_2013/healingStack01/";
subPath_unwounded_Hunter_2015_1 = "woundHealingVideos/Drosophila/Hunter_et_al_JCB_2015/Fig1_cropped_unwounded-01/";
subPath_unwounded_Hunter_2015_7 = "woundHealingVideos/Drosophila/Hunter_et_al_JCB_2015/Fig7_cropped_unwounded-01/";

allUnwoundedPaths = [subPath_unwounded_Abreu_Blanco, ...
    subPath_unwounded_Verboon,... % subPath_unwounded_Fernandez_Gonzalez,...
    subPath_unwounded_Hunter_2015_1, subPath_unwounded_Hunter_2015_7];
%allUnwoundedPaths = [];


allPaths = [allWoundedPaths, allUnwoundedPaths];

% for all wounds and unwoundeds, create arrays to hold mean and error
% values
woundMeanShapes = zeros(length(allWoundedPaths), 1);
woundStdShapes = zeros(length(allWoundedPaths),1);

unwoundedMeanShapes = zeros(length(allUnwoundedPaths),1);
unwoundedStdShapes = zeros(length(allUnwoundedPaths),1);

for woundIt=1:length(allPaths)
    filename = "handCorrection.tif";
    %fullPath = imageFolder + subPath_unwounded_Hunter_2015_7 + filename;
    fullPath = imageFolder + allPaths(woundIt) + filename;
    frameii = fullPath;

    if makeAMovie == 1
        moviestr = outputFolder+num2str(woundIt)+"segmentationTimeLapse.mp4";
        vobj = VideoWriter(moviestr, 'MPEG-4');
        vobj.Quality = 100;
            
        vobj.FrameRate = 5;
        open(vobj);
    end


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
    labeledImage = bwlabel(BWmask);
    
    % Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
    blobMeasurements = regionprops(labeledImage, BWmask, 'all');
    numberOfBlobs = size(blobMeasurements, 1);
    
    % bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
    figure(5)
    
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
    areaArraysTolerance = zeros(length(reducePolyTolerance), 20);
    % looping over tolerances, record shapes and areas to get standard error
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
            % discard nan edges due to segmentation
            vertices(isnan(vertices)) = [];
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
        areaMask = (a_array-median(a_array)<3.0*median(a_array));
        
        compositeMask = boxEdgeMask & areaMask;
        
        filteredPerimeter = pmeter_array(compositeMask);
        filteredArea = a_array(compositeMask);
        shapeParameters = filteredPerimeter.^2 ./ filteredArea / (4*pi);
        mean(shapeParameters)
    
        %save shapeParameters at this tolerance
        shapeArraysTolerance(ii,1:length(shapeParameters)) = shapeParameters;
        areaArraysTolerance(ii,1:length(shapeParameters)) = filteredArea;

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
        if makeAMovie == 1
            text(0.25,0.25, "\epsilon = "+num2str(reducePolyTolerance(ii)),'Units','Normalized');
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end
        %figure(7)
        %scatter(vertexPositions(:,2), vertexPositions(:,1),'rx')
        %daspect([1 1 1])
    end
    
    shapeArraysTolerance(:,~any(shapeArraysTolerance,1)) = [];
    areaArraysTolerance(:, ~any(areaArraysTolerance,1)) = [];
    shapeError = shapeArraysTolerance.*0;
    areaError = areaArraysTolerance.*0;

    % for each cell, compute error down the rows (1 row per tolerance, 1
    % col per cell)
    for cellIt=1:length(shapeArraysTolerance(1,:))
        shapeError(:,cellIt) = (shapeArraysTolerance(:,cellIt) - shapeArraysTolerance(1,cellIt))./shapeArraysTolerance(:,cellIt);
        areaError(:,cellIt) = (areaArraysTolerance(:,cellIt) - areaArraysTolerance(1,cellIt))./areaArraysTolerance(:,cellIt);
    end

    figure(8); clf; hold on
    plot(reducePolyTolerance,mean(areaError,2),"k", "DisplayName","% area change",'linewidth', 6)
    plot(reducePolyTolerance,mean(shapeError,2),"r","DisplayName","% shape change",'linewidth', 6)
    ylabel("Error")
    xlabel("$\epsilon$","interpreter","latex")
    legend('Location','best')
    saveas(gcf, outputFolder+num2str(woundIt)+"shapeSegmentedErrors.eps", 'epsc')

    [stdAreas,meanAreas] = std(areaArraysTolerance,0,2);
    [stdShapes,meanShapes] = std(shapeArraysTolerance,0,2);

    figure(9); clf; hold on
    %errorbar(reducePolyTolerance,meanAreas(:),stdAreas(:), "k", "DisplayName","area",'linewidth', 6)
    plot(reducePolyTolerance,meanAreas(:), "k", "DisplayName","area",'linewidth', 6)
    ylabel("area")
    yyaxis right
    %errorbar(reducePolyTolerance,meanShapes(:),stdShapes(:),"r","DisplayName","shape",'linewidth', 6)
    plot(reducePolyTolerance,meanShapes(:),"r","DisplayName","shape",'linewidth', 6)
    ylabel("shape")
    ylim([1 5])
    xlabel("$\epsilon$", "interpreter", "latex")
    legend
    saveas(gcf, outputFolder+num2str(woundIt)+"shapeSegmentedChanges.eps", 'epsc')

    figure(11); clf; hold on
    %errorbar(reducePolyTolerance,meanAreas(:),stdAreas(:), "k", "DisplayName","area",'linewidth', 6)
    plot(reducePolyTolerance,meanAreas(:), "k", "DisplayName","area",'linewidth', 6)
    ylabel("area")
    yyaxis right
    %errorbar(reducePolyTolerance,meanShapes(:),stdShapes(:),"r","DisplayName","shape",'linewidth', 6)
    plot(reducePolyTolerance,meanShapes(:),"r","DisplayName","shape",'linewidth', 6)
    for cellIt=1:length(shapeArraysTolerance(1,:))
        scatter(reducePolyTolerance, shapeArraysTolerance(:,cellIt),'r.','HandleVisibility','off')
    end
    ylabel("shape")
    ylim([1 5])
    xlabel("$\epsilon$", "interpreter", "latex")
    legend
    saveas(gcf, outputFolder+num2str(woundIt)+"shapeSegmentedChangesError1.eps", 'epsc')


    % note to self, current std is just the std among the cells, not really
    % a useful metric unless I'm interest in sample std rather than
    % measurement error, which I'm currently not.

    if makeAMovie == 1
        close(vobj);
    end
   
end