classdef GPS < handle
       
    % Class properties and variables
    properties
        
        % General parameters
        reference_latitude = 40.8623031194444;
        reference_longitude = -77.8362636138889;
        reference_altitude = 333.817;
        
    end
        
    properties (Constant)
        
        A_EARTH = 6378137;
        flattening = 1/298.257223563;
        NAV_E2 = (2 - GPS.flattening) * GPS.flattening; % also e^2
        
    end
    
    methods
        
        function obj = GPS(reference_latitude, reference_longitude, reference_altitude)
           
            if nargin == 3
                obj.reference_latitude = reference_latitude;
                obj.reference_longitude = reference_longitude;
                obj.reference_altitude = reference_altitude;
            else
                fprintf('Using default base station coordinates at LTI test track\n')
            end
            
        end
        
        function ENU = WGSLLA2ENU(obj, latitude, longitude, altitude, reference_latitude, reference_longitude, reference_altitude)
            
            if nargin == 4
               
                reference_latitude = obj.reference_latitude;
                reference_longitude = obj.reference_longitude;
                reference_altitude = obj.reference_altitude;
                
            end
            
            ENU = zeros(3,length(latitude));
            for i = 1:length(latitude)
            
                XYZ = WGSLLA2XYZ(obj, latitude(i), longitude(i), altitude(i));

                ENU(:,i) = WGSXYZ2ENU(obj, XYZ, reference_latitude, reference_longitude, reference_altitude);
                
            end
            
        end
        
        % Returns the equivalent WGS84 XYZ coordinates (in meters) for a given geodetic latitude (degrees),
        % longitude (degrees), and altitude above the WGS84 ellipsoid in meters.  Note: N latitude is positive,
        % S latitude is negative, E longitude is positive, W longitude is negative.
        %
        % Ref: Decker, B. L., World Geodetic System 1984,
        % Defense Mapping Agency Aerospace Center. 
        function XYZ = WGSLLA2XYZ(obj, latitude, longitude, altitude)
            
            slat = sin(deg2rad(latitude));
            clat = cos(deg2rad(latitude));
            r_n = GPS.A_EARTH / sqrt(1 - GPS.NAV_E2 * slat * slat);
            XYZ = [ (r_n + altitude) * clat * cos(deg2rad(longitude));
                    (r_n + altitude) * clat * sin(deg2rad(longitude));
                    (r_n * (1 - GPS.NAV_E2) + altitude) * slat ];

            if ((latitude < -90.0) || (latitude > 90.0) || (longitude < -180.0) || (longitude > 360.0))
                error('WGS lat or WGS lon out of range');
            end
            
        end
        
        function ENU = WGSXYZ2ENU(obj, XYZ, reference_latitude, reference_longitude, reference_altitude)
            
            if nargin == 2
               
                reference_latitude = obj.reference_latitude;
                reference_longitude = obj.reference_longitude;
                reference_altitude = obj.reference_altitude;
                
            end
           
            % First, calculate the xyz of reflat, reflon, refalt

            refXYZ = WGSLLA2XYZ(obj, reference_latitude, reference_longitude, reference_altitude);

            % Difference xyz from reference point

            diffXYZ = XYZ - refXYZ;

            % Now rotate the (often short) diffxyz vector to enu frame

            % Rotate about Z-axis
            R1 = eye(3);
            c = cos(deg2rad(90 + reference_longitude));
            s = sin(deg2rad(90 + reference_longitude));
            R1(1,1) = c;
            R1(2,2) = c;
            R1(2,1) = -s;
            R1(1,2) = s;
    
            % Rotate about X-axis
            R2 = eye(3);
            c = cos(deg2rad(90 - reference_latitude));
            s = sin(deg2rad(90 - reference_latitude));
            R2(2,2) = c;
            R2(3,3) = c;
            R2(2,3) = s;
            R2(3,2) = -s;
            
            R = R2 * R1;

            ENU = R * diffXYZ;
            
        end
        
        function XYZ = ENU2WGSXYZ(obj, ENU, reference_latitude, reference_longitude, reference_altitude)
            
            if nargin == 2
               
                reference_latitude = obj.reference_latitude;
                reference_longitude = obj.reference_longitude;
                reference_altitude = obj.reference_altitude;
                
            end
            
            % Rotate the ENU vector to XYZ frame
            
            % Rotate about Z-axis
            R1 = eye(3);
            c = cos(deg2rad(90 + reference_longitude));
            s = sin(deg2rad(90 + reference_longitude));
            R1(1,1) = c;
            R1(2,2) = c;
            R1(2,1) = -s;
            R1(1,2) = s;
    
            % Rotate about X-axis
            R2 = eye(3);
            c = cos(deg2rad(90 - reference_latitude));
            s = sin(deg2rad(90 - reference_latitude));
            R2(2,2) = c;
            R2(3,3) = c;
            R2(2,3) = s;
            R2(3,2) = -s;
            
            R = R2 * R1;
            
            % Rotate
            diffXYZ = R \ ENU;
            
            % Calculate the XYZ of the reference latitude, longitude, altitude
            refXYZ = WGSLLA2XYZ(obj, reference_latitude, reference_longitude, reference_altitude);
            
            % Add the difference to the reference
            XYZ = diffXYZ + refXYZ;
            
        end
        
        function LLA = WGSXYZ2LLA(obj, XYZ)
           
            if ((XYZ(1) == 0.0) && (XYZ(2) == 0.0))
                longitude = 0.0;
            else
                longitude = rad2deg(atan2(XYZ(2), XYZ(1)));
            end
            
            if ((XYZ(1) == 0.0) && (XYZ(2) == 0.0) && (XYZ(3) == 0.0))
                error('WGS XYZ located at the center of the earth')
            else
                
                % Make initial latitude and altitude guesses based on
                % spherical earth.
                
                rhosqrd = XYZ(1)^2 + XYZ(2)^2;
                rho = sqrt(rhosqrd);
                templat = atan2(XYZ(3), rho);
                tempalt = sqrt(rhosqrd + XYZ(3)^2) - GPS.A_EARTH;
                rhoerror = 1000.0;
                zerror   = 1000.0;
                
                %  Newton's method iteration on templat and tempalt makes
                %	the residuals on rho and z progressively smaller.  Loop
                %	is implemented as a 'while' instead of a 'do' to simplify
                %	porting to MATLAB

                while ((abs(rhoerror) > 1e-6) || (abs(zerror) > 1e-6)) 
                    slat = sin(templat);
                    clat = cos(templat);
                    q = 1 - GPS.NAV_E2 * slat * slat;
                    r_n = GPS.A_EARTH / sqrt(q);
                    drdl = r_n * GPS.NAV_E2* slat * clat / q; % d(r_n)/d(latitutde) 

                    rhoerror = (r_n + tempalt) * clat - rho;
                    zerror   = (r_n * (1 - GPS.NAV_E2) + tempalt) * slat - XYZ(3);

                    %			  --                               -- --      --
                    %			  |  drhoerror/dlat  drhoerror/dalt | |  a  b  |
                                % Find Jacobian           |		       		    |=|        |
                    %			  |   dzerror/dlat    dzerror/dalt  | |  c  d  |
                    %			  --                               -- --      -- 

                    aa = drdl * clat - (r_n + tempalt) * slat;
                    bb = clat;
                    cc = (1 - GPS.NAV_E2) * (drdl * slat + r_n * clat);
                    dd = slat;

                    %  Apply correction = inv(Jacobian)*errorvector

                    invdet = 1.0 / (aa * dd - bb * cc);
                    templat = templat - invdet * (dd*rhoerror - bb*zerror);
                    tempalt = tempalt - invdet * (-cc*rhoerror + aa*zerror);
                end

                latitude = rad2deg(templat);
                altitude = tempalt;
                
            end
            
            LLA = [latitude; longitude; altitude];
            
        end

        function LLA = ENU2WGSLLA(obj, ENU, reference_latitude, reference_longitude, reference_altitude)
            
            if nargin == 2
               
                reference_latitude = obj.reference_latitude;
                reference_longitude = obj.reference_longitude;
                reference_altitude = obj.reference_altitude;
                
            end
            
            LLA = zeros(3,size(ENU,2));
            for i = 1:size(ENU,2)
               
                XYZ = ENU2WGSXYZ(obj, ENU(:,i), reference_latitude, reference_longitude, reference_altitude);
                
                LLA(:,i) = WGSXYZ2LLA(obj, XYZ);
                
            end
            
        end
        
    end
    
end