% read in a segmentation and try to fit polygons to each cell
% fixed points of the polygons are determined by junctions between 
% 3 or more cells

close all
clear
set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultAxesFontSize', 24);

rng(1)
addpath '/Users/AndrewTon/Documents/YalePhD/projects/imageProcessing/';
addpath '/Users/AndrewTon/Documents/YalePhD/projects/imageProcessing/testImagesForShapeCalibration'

% convert tif output of tissue analyzer into polygons
% note - install David Legland's matGeom toolbox to do polygon
% interpolation calculations

%parameter list
outputInterpolatedConfigurations = 0;
outputFileID = fopen('testDPMConfigurationFile.txt','w');
%tolerance for reducePoly algorithm

boxEdgeMaskOff = 1; % if true, disable boxEdgeMask

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

%allWoundedPaths = [subPath_Wood_1, subPath_Tetley, subPath_Campos, ...
%    subPath_Garcia_Fernandez, subPath_Abreu_Blanco, subPath_Stramer,...
%    subPath_Verboon, subPath_Wood_2];

%allWoundedPaths = [subPath_Campos];
allWoundedPaths = [];

subPath_unwounded_Abreu_Blanco = "woundHealingVideos/Drosophila/Abreu_Blanco_et_al_Cell_Science_2012/healingStack01_cropped/";
subPath_unwounded_Verboon = "woundHealingVideos/Drosophila/Verboon_et_al_Small_GTPases_2015/Verboon_Fig_1A_unwounded_cropped-01/";
%subPath_unwounded_Fernandez_Gonzalez = "woundHealingVideos/Drosophila/Fernandez_Gonzalez_et_al_MBOC_2013/healingStack01/";
subPath_unwounded_Hunter_2015_1 = "woundHealingVideos/Drosophila/Hunter_et_al_JCB_2015/Fig1_cropped_unwounded-01/";
subPath_unwounded_Hunter_2015_7 = "woundHealingVideos/Drosophila/Hunter_et_al_JCB_2015/Fig7_cropped_unwounded-01/";
subPath_after_wounding_Verboon = "woundHealingVideos/Drosophila/Verboon_et_al_Small_GTPases_2015/Verboon_Fig_1A_wounded_cropped-01/";


%allUnwoundedPaths = [subPath_unwounded_Abreu_Blanco, ...
%    subPath_unwounded_Verboon,... % subPath_unwounded_Fernandez_Gonzalez,...
%    subPath_unwounded_Hunter_2015_1, subPath_unwounded_Hunter_2015_7];
%allUnwoundedPaths = [subPath_after_wounding_Verboon];
allUnwoundedPaths = [];

subPath_duringWounding_Fernandez_Gonzalez = "woundHealingVideos/Drosophila/" + ...
    "Fernandez_Gonzalez_et_al_MBOC_2013/Fernandez_Gonzalez_et_al_MBOC_20130035/";
allDuringWoundingPaths = [subPath_duringWounding_Fernandez_Gonzalez];
%allDuringWoundingPaths = [];

allPaths = [allWoundedPaths, allUnwoundedPaths, allDuringWoundingPaths];

% for all wounds and unwoundeds, create arrays to hold mean and error
% values
woundMeanShapes = zeros(length(allWoundedPaths), 1);
woundStdShapes = zeros(length(allWoundedPaths),1);

unwoundedMeanShapes = zeros(length(allUnwoundedPaths),1);
unwoundedStdShapes = zeros(length(allUnwoundedPaths),1);

for woundIt=1:length(allPaths)
    filename = "handCorrection.tif";
    fullPath = imageFolder + allPaths(woundIt) + filename;
    frameii = fullPath;

    originalImage = imread(frameii);
    
    %figure(1)
    %imshow(originalImage)
    %set(gca, 'YDir', 'normal')
    gray1 = im2gray(originalImage);
    
    figure(2); hold on
    [BWmask, maskedImage1] = createMask_lighter(originalImage);
    %BWmask(1,:) = 0;
    %BWmask(:,1) = 0;
    %BWmask(:,end) = 0;
    %BWmask(end,:) = 0;
    imshow(BWmask)
    h=gca;
    h.Visible = 'On';

    %imshow(bwmorph(BWmask,'branchpoints'))
    branchPoints = bwmorph(BWmask,'branchpoints');
    branchCoords = [];
    for branch_ii=1:length(branchPoints(:,1))
        for branch_jj=1:length(branchPoints(1,:))
            if (branchPoints(branch_jj,branch_ii) == 1)
                branchCoords = [branchCoords ; branch_jj branch_ii];
            end
        end
    end
    %scatter(branchPoints(:,2),branchPoints(:,1),'ro','filled')
    scatter(branchCoords(:,2),branchCoords(:,1),'ro','filled')

    %% BWmask is the base image to work with using blob counting
    % get labels for 4-connected objects in BWmask (polygons/blobs), and each
    % boundary for those blobs
    [B,labeledImage] = bwboundaries(~BWmask,4);
    
    % Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
    blobMeasurements = regionprops(labeledImage, BWmask, 'all');
    numberOfBlobs = size(blobMeasurements, 1);
    blobCentroids = zeros(numberOfBlobs,2);
    for ii=1:numberOfBlobs
        blobCentroids(ii,:) = [blobMeasurements(ii).Centroid];
    end

    %plot blob centroids
    scatter(blobCentroids(:,1),blobCentroids(:,2),'bx')

    % plot boundary of a cellID
    % plot(boundaries{cellID}(:,2),boundaries{40}(:,1),'yellow')
    
    % plot each blob and assign each blob the coordinates of
    % the branch points on its boundary
    branchMap = zeros(size(branchPoints));
    numBranchPoints = sum(branchPoints,'all');
    I = sub2ind(size(branchMap),branchCoords(:,1),branchCoords(:,2));
    branchMap(I) = 100;
    figure(4)
    labeledImageWithBranchPoints = labeledImage + branchMap;
    heatmap(labeledImageWithBranchPoints);
    colormap(jet)

    fixedPoints = cell(numberOfBlobs,1);
    for ii=1:numberOfBlobs
        % boundary_ii_thickened is a fuzzy version of B, which we'll use to
        % detect when branchCoords are near B
        boundary_ii_thickened = [B{ii}+[1 1];B{ii}+[1 0];B{ii}+[1 -1]...
            ;B{ii}+[0 1];B{ii};B{ii}+[0 -1]...
            ;B{ii}+[-1 1];B{ii}+[-1 0];B{ii} + [-1 -1]];
        for jj=1:numBranchPoints
            % keep branchCoords if it borders B, since branchCoords are
            % exterior points and B are interior points
            if (ismember(branchCoords(jj,:), boundary_ii_thickened,'rows'))
                fixedPoints{ii} = [fixedPoints{ii} ;branchCoords(jj,:)];
            end
        end
    end

    % if cellID has fixed points on the edge of the image, ignore that cell
    % if cellID has more than 1 region, then it'll be hard to sort the
    % vertices, so ignore that cell too.
    cellsToIgnore = zeros(numberOfBlobs,1);
    for ii=1:numberOfBlobs
        xBounds = fixedPoints{ii}(:,1);
        yBounds = fixedPoints{ii}(:,2);
        if (max(ismember(xBounds, 1))) % if any xBounds are 1, then ignore this cell
            cellsToIgnore(ii) = 1;
        elseif (max(ismember(yBounds,1)))
            cellsToIgnore(ii) = 1;
        elseif (max(ismember(xBounds,length(BWmask(1,:))))) % not sure if i should swap this and the next one
            cellsToIgnore(ii) = 1;
        elseif (max(ismember(yBounds,length(BWmask(:,1)))))
            cellsToIgnore(ii) = 1;
        elseif (polyshape(B{ii}).NumRegions ~= 1)
            cellsToIgnore(ii) = 1;
        end
    end

    % bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
    figure(5)
    
    hold on;
    % store box dimensions for boxEdgeMask, which throws out polygons too close
    % to the edge
    boxDims = [length(BWmask(1,:)) length(BWmask(:,1))];
    
    boundaries = B;
    numberOfBoundaries = size(boundaries, 1);
    a_array = zeros(numberOfBoundaries,1);
    pmeter_array = zeros(numberOfBoundaries,1);
    centroid_array = zeros(numberOfBoundaries,2);
    isValidCell = ones(numberOfBoundaries,1,'logical');
 
 %  record shapes and areas to get standard error
    vertexPositions = [];
    pixelBuffer = 1;
    for kk = 1 : numberOfBoundaries
        if (cellsToIgnore(kk))
            isValidCell(kk) = 0;
            continue
        end
        thisBoundary = boundaries{kk};
    
        %poly_kk = polyshape(P_reduced(:,2),P_reduced(:,1));
        poly_kk = polyshape(thisBoundary);
        plot(poly_kk,'EdgeColor',[rand rand rand])
        a_array(kk) = area(poly_kk);
        pmeter_array(kk) = perimeter(poly_kk);
        vertices = poly_kk.Vertices;
        % discard nan edges due to segmentation
        vertices(isnan(vertices)) = [];
        centroid_array(kk,:) = mean(vertices);
    end
    
    filteredPerimeter = pmeter_array(isValidCell);
    filteredArea = a_array(isValidCell);
    shapeParameters = pmeter_array.^2 ./ a_array / (4*pi);

    % centroid array marks the centers of each polygon visually
    scatter(centroid_array(isValidCell,1), centroid_array(isValidCell,2),'k','x')

    figure(6)
    % show a cell, its centroid, and its fixed points
    [~,validCellID] = max(isValidCell); 
    %validCellID=42;
    scatter(fixedPoints{validCellID}(:,1),fixedPoints{validCellID}(:,2),'rx')
    hold on
    scatter(centroid_array(validCellID,1),centroid_array(validCellID,2),'bo')
    plot(polyshape(boundaries{validCellID}).Vertices(:,1),polyshape(boundaries{validCellID}).Vertices(:,2),'k')
    daspect([1 1 1])

    % adjust fixed points so they coincide with the nearest boundary pixel
    adjustedFixedPoints = fixedPoints;
    fixedPointIndices = cell(numberOfBoundaries,1);
    orderedFixedPointIndices = fixedPointIndices;

    for ii=1:numberOfBoundaries
        if (cellsToIgnore(ii))
            continue
        end
        this_poly = polyshape(boundaries{ii});
        for jj=1:length(fixedPoints{ii}(:,1))
            pixelDiffs = this_poly.Vertices - fixedPoints{ii}(jj,:);
            pixelDistSq = pixelDiffs(:,1).^2 + pixelDiffs(:,2).^2;
            [~,argmin] = min(pixelDistSq);
            nearestPixel = this_poly.Vertices(argmin,:);
            adjustedFixedPoints{ii}(jj,:) = nearestPixel;
            fixedPointIndices{ii} = [fixedPointIndices{ii} argmin];
        end
        % put the indices of fixed points in order so they can be
        % correlated with pixel locations in polyshape.Vertices
        orderedFixedPointIndices{ii} = sort(unique(fixedPointIndices{ii}));
    end
    
    scatter(adjustedFixedPoints{validCellID}(:,1),adjustedFixedPoints{validCellID}(:,2),'bx')

    % next step: feed adjustedFixedPoints{ii} along with boundaries{ii}
    % into polygonCurveFit.m
    % custom function polygonCurveFit, computes shape using curve
    % fitting of degree polynomialDegree between each specified fixed 
    % point index in startingIndices, with optional plotting bool arg
    polynomialDegree = 10;
    shape = polygonCurveFit(polyshape(boundaries{validCellID}),orderedFixedPointIndices{validCellID},polynomialDegree,true);
    figure(7); hold on
    for ii=1:numberOfBoundaries
        if (cellsToIgnore(ii))
            continue
        end
        this_poly = polyshape(boundaries{ii});
        plot(this_poly)
        scatter(adjustedFixedPoints{ii}(:,1),adjustedFixedPoints{ii}(:,2),'bx')
        shape = polygonCurveFit(polyshape(boundaries{ii}),orderedFixedPointIndices{ii},polynomialDegree,true);
    end

%     % for each valid segmented cell, record vertex positions separated by nans
%     % and a header, sorted increasing by cell perimeters
%     figure(6)
%     clf; hold on;
%     set(gca, 'YDir', 'reverse')
% 
%     % minPerimeter will be split up into ~20 vertices, and the rest
%     % will use the corresponding particle radius from this splitting
%     minPerimeter = min(filteredPerimeter);
%     vertDiameter = minPerimeter / 20.0;
% 
%     
%     numDudCells = 0;
%     for kk = 1:numberOfBoundaries
%         jj = sort_order(kk);
%         thisBoundary = boundaries{jj};
%         if (isValidCell(jj) == false)
%             numDudCells = numDudCells + 1;
%             continue;
%         end
%         % fill vertexPositions with a header
%         vertexPositions = [vertexPositions; nan nan];
%         vertexPositions = [vertexPositions; kk length(thisBoundary(:,1))];
%         for ll = 1:length(thisBoundary(:,1))
%             %fill vertexPositions with lines of coordinates
%             vertexPositions = [vertexPositions; thisBoundary(ll,1) thisBoundary(ll,2)];
%         end
%     end
end