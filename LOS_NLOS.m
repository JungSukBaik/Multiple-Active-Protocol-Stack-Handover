function [ NLOSflag ] = LOS_NLOS(testpoints,XYCoordinates,klast,count,XYcenter,MD)
%%
% Author: Mandar N. Kulkarni, Graduate Student, UT Austin
% Contact: mandar.kulkarni@utexas.edu

%In order to use this code, please cite the following paper:
%M. N. Kulkarni, S. Singh and J. G. Andrews, "Coverage and Rate Trends in Dense Urban mmWave Cellular Networks", IEEE Globecom, Austin, Dec. 2014. 
%OR
%S. Singh, M. N. Kulkarni, A. Ghosh, J. G. Andrews, "Tractable Model for Rate in Self-backhauled mmWave Cellular Networks", submitted to IEEE JSAC, July 2014. [Available Online] http://arxiv.org/abs/1407.5537


%   Checks if the line joining two points is intersected by 
% � any of the buildings.

%   testpoints is (2 x 2) matrix where
%   testpoints(:,1) is the first endpoint and
%   testpoints(:,2) is the the second endpoint of the line segment.

% XYCoordinates: Has coordinates of all buildings. It is a 3D array.
%    The indices are as follows BuildingXYCoord[building number][vertex
%    number][X or Y coordinate]
% klast: This array stores the number of vertices of the building+1 (one of the vertex is counted twice)
% count: Total number of buildings in the area
% XYcenter: Contains the centroid of each building
% MD: maximum distance of any vertex of any building from its 
%    centroid

% Helpful links: 
% http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
%
%
%   First we shift our origin to the first point
FirstPointNOTSHIFTED = testpoints(:,1);
SecondPointSHIFTED = testpoints(:,2)-FirstPointNOTSHIFTED;
NLOSflag=0;
%Instead of searching over all buildings just search over nearby buildings
%that could actually intersect the line segment based on the building
%statistics (maximum distance between the centroid of a building and the line segment)
%
total_count = 0;
for j=1:count 
    lengthseg = norm(SecondPointSHIFTED);
    DistanceCentroidFromSeg =  abs(det([SecondPointSHIFTED,XYcenter(:,j)-FirstPointNOTSHIFTED]))/lengthseg; 
    DistanceCentroidFromSecondEndpoint = norm(XYcenter(:,j)-FirstPointNOTSHIFTED-SecondPointSHIFTED);
    DistanceCentroidFromFirstEndpoint = norm(XYcenter(:,j)-FirstPointNOTSHIFTED);
    if DistanceCentroidFromSeg<=MD && (DistanceCentroidFromSecondEndpoint<=lengthseg+MD && DistanceCentroidFromFirstEndpoint<=lengthseg+MD)
        total_count= total_count+1;
        cnt_opt(total_count) = j;
    end
end

for j=1:total_count
        P = [XYCoordinates(cnt_opt(j),1:(klast(cnt_opt(j))),1)' XYCoordinates(cnt_opt(j),1:(klast(cnt_opt(j))),2)']';
        M = bsxfun(@minus, P, FirstPointNOTSHIFTED);
        tester = seg2poly(FirstPointNOTSHIFTED,SecondPointSHIFTED,M);
        if isempty(tester)~=1 
            NLOSflag = 1; break;
        end
end
end