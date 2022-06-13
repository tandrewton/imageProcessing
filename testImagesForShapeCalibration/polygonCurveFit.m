% input: polygon pgon and list M with around 6 vertex indices 
%           which represent the turning points, and we'll fit curves
%           to all other N-6 points
% output: estimation of shape parameter = perimeter^2/(4*pi*area)
function shape = polygonCurveFit(polygon, fixedPointIndices, polynomialDegree, boolPlot)
    arc_length = 0;
    n = length(fixedPointIndices);
    % for each startingIndex to each endingIndex, fit a polynomial with fixed
    % endpoints at startingIndex and endingIndex. Plot them on top of the
    % original noisy polygon.
    sortedPoints = polygon.Vertices;
    fixedPoints = sortedPoints(fixedPointIndices,:);
    for ii=1:n
        appendedPoints = [sortedPoints ; sortedPoints]; 
        % append double sortedPoints for indexing purposes
        startingIndex = fixedPointIndices(ii);
        shiftedStartingIndex = circshift(fixedPointIndices,-1);
        endingIndex = shiftedStartingIndex(ii);
        if (startingIndex - endingIndex > 0)
            endingIndex = endingIndex + length(sortedPoints);
        end

        appendedPoints(startingIndex,1)
        appendedPoints(endingIndex, 1)
        
        % decide whether to fit X(Y) or Y(X) 
        B = appendedPoints(startingIndex:endingIndex,:);
        numOverlapsInX = length(B(:,1))-length(unique(B(:,1)))
        numOverlapsInY = length(B(:,2))-length(unique(B(:,2)))

        if (numOverlapsInY >= numOverlapsInX) % fit Y(X)
            [p,~] = polyfix(appendedPoints(startingIndex:endingIndex,1),...
                appendedPoints(startingIndex:endingIndex,2),polynomialDegree,... 
                [fixedPoints(ii,1) fixedPoints(mod(ii,n)+1,1)], ...
                [fixedPoints(ii,2) fixedPoints(mod(ii,n)+1,2)]);
    
            polygonLegX = sort(appendedPoints(startingIndex:endingIndex,1));
            % extract X values to evaluate the polynomial. discard points above and
            % below the endpoints
            minXValue = min(appendedPoints(startingIndex,1),appendedPoints(endingIndex,1));
            maxXValue = max(appendedPoints(startingIndex,1),appendedPoints(endingIndex,1));
            polygonLegX = polygonLegX(polygonLegX >= minXValue);
            polygonLegX = polygonLegX(polygonLegX <= maxXValue);
        elseif (numOverlapsInY < numOverlapsInX) % fit X(Y)
            [p,~] = polyfix(appendedPoints(startingIndex:endingIndex,2),...
                appendedPoints(startingIndex:endingIndex,1),polynomialDegree,... 
                [fixedPoints(ii,2) fixedPoints(mod(ii,n)+1,2)], ...
                [fixedPoints(ii,1) fixedPoints(mod(ii,n)+1,1)]);
    
            polygonLegX = sort(appendedPoints(startingIndex:endingIndex,2));
            % extract X values to evaluate the polynomial. discard points above and
            % below the endpoints
            minXValue = min(appendedPoints(startingIndex,2),appendedPoints(endingIndex,2));
            maxXValue = max(appendedPoints(startingIndex,2),appendedPoints(endingIndex,2));
            polygonLegX = polygonLegX(polygonLegX >= minXValue);
            polygonLegX = polygonLegX(polygonLegX <= maxXValue);
        end

        if (boolPlot == 1 && numOverlapsInY >= numOverlapsInX)
            plot(polygonLegX, polyval(p,polygonLegX), 'r')
        elseif (boolPlot == 1 && numOverlapsInY < numOverlapsInX)
            plot(polyval(p,polygonLegX), polygonLegX, 'r')
        end

        p_derivative = p(1:end-1).*(length(p)-1:-1:1); %dp/dx
        dydx_sq = polyval(p_derivative, polygonLegX).^2;
        %p(1:end-1)
        %(length(p)-1:-1:1)
        %[polygonLegX polyval(p,polygonLegX)]

        I = trapz(polygonLegX, sqrt(1+dydx_sq));
        arc_length = arc_length + I;
    end 
    shape = arc_length^2/ (4 * pi * area(polygon));
end
