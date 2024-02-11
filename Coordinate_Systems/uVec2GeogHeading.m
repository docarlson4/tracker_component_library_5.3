function [geoEastOfNorth,angUpFromLevel]=uVec2GeogHeading(point,u)
%%UVEC2GEOGHEADING Obtain the geographic headings in radians East of true
%                 North as well as the elevations above the local tangent
%                 plane according to a particular reference ellipsoid that
%                 correspond to directions specified by unit vectors in
%                 ECEF coordinates. Direction vectors that are perfectly
%                 orthogonal to the local tangent plane will result in 0
%                 being returned for the heading.
%
%INPUTS: point The location of the point in geodetic latitude and longitude
%              in radians at which the headings are taken. Point can be
%              [latitude;longitude] or [latitude;longitude;height]. The
%              height component is ignored if included since it does not
%              change the result.
%            u A 3XN matrix of unit vectors pointing in the directions that
%              are to be converted to heading in radians East of North at
%              the given point.
%
%OUTPUTS: geoEastOfNorth An NX1 array of N geographic headings in radians
%                clockwise from North (East of North) that correspond to
%                the unit vectors given in u at the point point. For u
%                values that are perfectly orthogonal to the local tangent
%                plane, 0 values are returned.
%        angUpFromLevel An NX1 array of N elevations in radians above the
%                local tangent plane.
%
%The components of the vector u in East and North directions can be
%obtained by projecting u onto the East and North unit vectors, which are
%obtained using the function getENUAxes. Given the projection East and
%North (two sides of a triangle), the four-quadrant inverse tangent
%function can be used to determine the geographic heading. Similarly, the
%sin of the elevation is given by the dot product of the given unit
%vector and the local up vector.
%
%EXAMPLE:
%Here, we show that the transformation is consistent with tis inverse:
% plhPoint=[0.1;-0.3;4000];
% geoEastOfNorth=35*(pi/180);
% angUpFromLevel=15*(pi/180);
% 
% uVec=geogHeading2uVec(plhPoint,geoEastOfNorth,angUpFromLevel);
% [geoEastOfNorthBack,angUpFromLevelBack]=uVec2GeogHeading(plhPoint,uVec);
% RelErrGeo=(geoEastOfNorth-geoEastOfNorthBack)/geoEastOfNorth
% RelErrEl=(angUpFromLevel-angUpFromLevelBack)/angUpFromLevel
%Both relative errors should indicate more than 15 digits of precision,
%which is consistent with finite precision limits.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

uLocal=getENUAxes(point,false);
uEast=uLocal(:,1);
uNorth=uLocal(:,2);
uUp=uLocal(:,3);

%Perform dot products over all of the vectors.
magEast=sum(bsxfun(@times,u,uEast),1)';
magNorth=sum(bsxfun(@times,u,uNorth),1)';
geoEastOfNorth=atan2(magEast,magNorth);
angUpFromLevel=asin(sum(bsxfun(@times,u,uUp),1))';

end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
