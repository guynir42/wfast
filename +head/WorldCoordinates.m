classdef WorldCoordinates < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        WCSAXES = 2;
        RADECSYS = 'ICRS';
        LONPOLE = NaN;
        LATPOLE = NaN;
        EQUINOX = NaN;
        CTYPE1 = 'RA---TPV';
        CTYPE2 = 'DEC--TPV';
        CUNIT1 = 'deg';
        CUNIT2 = 'deg';
        CRPIX = []; % pixel position of transformation anchor 
        CRVAL = []; % RA/Dec position of transformation anchor 
        CDELT = []; % shift between anchors
        CD = []; % rotation matrix (and pixel scale)
        PV = []; % values of the polynomials needed for TPV transformation. Column 1 is for X, column 2 is for Y
        
        RA_deg_center;
        DE_deg_center;
        
    end
    
    properties(Dependent=true)
        
        rotation;
        
    end
    
    properties % switches/controls
        
        use_tpv = 1;
        max_pv_indices = 40;
        
        debug_bit = 0;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = WorldCoordinates(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.WorldCoordinates')
                if obj.debug_bit>1, fprintf('WorldCoordinates copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isstruct(varargin{1}) && isfield(varargin{1}, 'CRVAL') % can add more checks
                if obj.debug_bit>1, fprintf('WorldCoordinates struct based constructor v%4.2f\n', obj.version); end
                obj.input(varargin{1});
            elseif ~isempty(varargin) && (isa(varargin{1}, 'SIM') || isa(varargin{1},'ClassWCS'))
                if obj.debug_bit>1, fprintf('WorldCoordinates object based constructor v%4.2f\n', obj.version); end
                obj.input(varargin{1});
            else
                if obj.debug_bit, fprintf('WorldCoordinates constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function val = get.rotation(obj)
            
            if isempty(obj.PV) || size(obj.PV,1)<3
                val = [];
            else
                val = acotd(obj.PV(2,1)./obj.PV(3,1));
            end
            
        end
        
        function val = SCALE_DEG(obj) % plate scale in degrees
            
            val = obj.CD(1); % for different scales on x and y we will also need CD(4)
            
        end
        
%         function val = get.RA_deg_center(obj)
%             
%             if isempty(obj.PV) || isempty(obj.CRVAL)
%                 val = [];
%             else
%                 val = obj.CRVAL(1) - obj.PV(1,1);
%             end
%             
%         end
%         
%         function val = get.DE_deg_center(obj)
%             
%             if isempty(obj.PV) || isempty(obj.CRVAL)
%                 val = [];
%             else
%                 val = obj.CRVAL(2) - obj.PV(1,2);
%             end
%             
%         end
%         
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, w)
            
            if isa(w, 'ClassWCS')
                w = w.WCS;
            elseif isa(w, 'SIM')
                w = w.WCS.WCS;
            elseif isstruct(w)
                % pass
            end
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isfield(w, list{ii})
                    obj.(list{ii}) = w.(list{ii});
                end
                
            end
            
            if ~isempty(w.CTYPE) && iscell(w.CTYPE)
                obj.CTYPE1 = w.CTYPE{1};
                obj.CTYPE2 = w.CTYPE{2};
            end
            
            if ~isempty(w.CUNIT) && iscell(w.CUNIT)
                obj.CUNIT1 = w.CUNIT{1};
                obj.CUNIT2 = w.CUNIT{2};
            end
            
            if isfield(w, 'tpv')
                obj.parsePVstruct(w.tpv);
            elseif isfield(w, 'PV')
                obj.parsePVstruct(w.PV);
            end
            
        end
        
        function parsePVstruct(obj, pv)
            
            obj.PV = NaN(obj.max_pv_indices, 2);
            
            for ii = 1:2
            
                for jj = 1:obj.max_pv_indices

                    idx = find(pv.Ind(:,1)==ii & pv.Ind(:,2)==jj-1);
                    
                    if ~isempty(idx)
                        obj.PV(jj,ii) = pv.KeyVal{idx(1)};
                    end
                    
                end
            
            end
            
            vec = obj.xy2coo(obj.CRPIX(1), obj.CRPIX(2));
            
            obj.RA_deg_center = vec(1);
            obj.DE_deg_center = vec(2);
            
        end
        
        function RA_Dec = xy2coo(obj, x, y) % find the RA and Dec for the image xy positions
        % Usage: RA_Dec = xy2coo(obj, x, y)
        % Find the RA and Dec for the image xy positions, using the TPV 
        % coefficients in the WCS. 
        
        if nargin==0, help('head.WorldCoordinates.xy2coo'); return; end
        
            if nargin<3 || isempty(y)
                if size(x,2)==2
                    y = x(:,2);
                    x = x(:,1);
                else
                    error('Must input x and y!');
                end
            end
            
            x2 = x - obj.CRPIX(1);
            y2 = y - obj.CRPIX(2);
            
            vec = obj.CD*[x2;y2]; % rotation matrix, also stretches the pixels to units of degrees
            x3 = vec(1); % units of degrees
            y3 = vec(2); % units of degrees
            
            if obj.use_tpv
            
                if ~isempty(strfind(obj.CTYPE1, 'TPV')) && ~isempty(strfind(obj.CTYPE2, 'TPV'))

                    R  = sqrt(x3.^2 + y3.^2); % units of degrees
            
                    [Xpowers, Ypowers] = head.WorldCoordinates.tpv_polydef;
            
                    % these lines assume PV is a 2x40 coefficients matrix with
                    % NaNs or zeros where there is no contribution to that polynomial
                    Xout = sum(obj.PV(:,1).*x3.^Xpowers(:,1).*y3.^Xpowers(:,2).*R.^Xpowers(:,3), 1, 'omitnan');
                    Yout = sum(obj.PV(:,2).*x3.^Ypowers(:,1).*y3.^Ypowers(:,2).*R.^Ypowers(:,3), 1, 'omitnan');

                else
                    error('Unknown CTYPE1 "%s". Use tpv instead...', obj.CTYPE1);
                end

            end
            
            [RA, Dec] = head.WorldCoordinates.pr_ignomonic(Xout, Yout, obj.CRVAL);

            RA_Dec(1) = RA;
            RA_Dec(2) = Dec;

        end
        
        function XY = coo2xy(obj, RA, Dec) % invertion of xy2coo by calling fminsearch and calculating the point nearest to RA/Dec
        % Usage: xy = coo2xy(obj, RA, Dec)
        % Use the GAIA match to transform the given RA/Dec into xy on the 
        % image plane. 
        % Specify coordinates as hexagesimal strings or numeric degrees, 
        % DO NOT GIVE RA AS NUMERIC HOURS, if you give it as numeric value, 
        % it must be in degrees!
        
            if ischar(RA)
                RA = head.Ephemeris.hour2deg(RA);
            end
            
            if ischar(Dec)
                Dec = head.Ephemeris.sex2deg(Dec);
            end
            
            func = @(xy) sum(([RA, Dec] - obj.xy2coo(xy(1), xy(2)) ).^2); % minimization function
            
            XY = fminsearch(func, obj.CRPIX); % initial guess is middle of field 
            
        end
        
    end
    
    methods % utilities
        
        function struct2obj(obj, s) % must first construct an object and then use the struct to update its properties
            
            list = properties(obj);
            
            for ii = 1:length(list)
               
                if isfield(s, list{ii})
                    
                    if isobject(obj.(list{ii}))
                        
                        if isprop(obj.(list{ii}), 'struct2obj')
                            struct2obj(obj.(list{ii}), s.list{ii}); % call the sub-object method if it exists
                        else
                            
                            list2 = properties(obj.(list{ii}));
                            
                            for jj = 1:length(list2)
                            
                                if isfield(s.(list{ii}), list2{jj})
                                    
                                    try 
                                        obj.(list{ii}).(list2{jj}) = s.(list{ii}).(list2{jj}); 
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    else % regular data types
                        try % we can run into many problems trying to set unsetable properties
                            obj.(list{ii}) = s.(list{ii}); 
                        end
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods(Static=true) % some methods stolen from Eran's ClassWCS and celetial.proj
        
        function [XiTab,EtaTab]=tpv_polydef
            % Return the TPV transformation polynomial orders definitions
            % Stolen shamelessly from Eran's Package: @ClassWCS
            % Description: Return the TPV transformation polynomial orders
            %              definitions.
            % Input  : null
            % Output : - A three column matrix of the xi' polynomial order
            %            in xi, eta, r,
            %            where PV index 0 refers to the first table line.
            %          - A three column matrix of the eta' polynomial order
            %            in xi, eta, r,
            %            where PV index 0 refers to the first table line.
            % Example: [XiTab,EtaTab]=ClassWCS.tpv_polydef
            
            XiTab  = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 0 2 0; 3 0 0; 2 1 0; 1 2 0; 0 3 0; 0 0 3; 4 0 0; 3 1 0; 2 2 0; 1 3 0; 0 4 0; 5 0 0; 4 1 0; 3 2 0; 2 3 0; 1 4 0; 0 5 0; 0 0 5; 6 0 0; 5 1 0; 4 2 0; 3 3 0; 2 4 0; 1 5 0; 0 6 0; 7 0 0; 6 1 0; 5 2 0; 4 3 0; 3 4 0; 2 5 0; 1 6 0; 0 7 0; 0 0 7];
            EtaTab = [0 0 0; 0 1 0; 1 0 0; 0 0 1; 0 2 0; 1 1 0; 2 0 0; 0 3 0; 1 2 0; 2 1 0; 3 0 0; 0 0 3; 0 4 0; 1 3 0; 2 2 0; 3 1 0; 4 0 0; 0 5 0; 1 4 0; 2 3 0; 3 2 0; 4 1 0; 5 0 0; 0 0 5; 0 6 0; 1 5 0; 2 4 0; 3 3 0; 4 2 0; 5 1 0; 6 0 0; 0 7 0; 1 6 0; 2 5 0; 3 4 0; 4 3 0; 5 2 0; 6 1 0; 7 0 0; 0 0 7];
        
        end
        
        function [Long,Lat]=pr_ignomonic(X,Y,CenterVec)
            % Project coordinates using the inverse Gnomonic non conformal projection
            % Stolen shamelessly from Eran's Package: celestial.proj
            % Description: Project coordinates (X and Y) using the
            %              inverse Gnomonic non conformal projection,
            % Input  : - Vector of X, in degrees (not radians!).
            %          - Vector of Y, in degrees (not radians!).
            %          - Central coordinate vector [Long_center,Lat_center],
            % Output : - Vector of longitude in degrees (not radians!).
            %          - Vector of latitude in degrees (not radians!).
            % See also: pr_gnomonic.m
            % Tested : Matlab 5.3
            %     By : Eran O. Ofek                    Jun 2005
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [Long,Lat]=celestial.proj.pr_ignomonic(1,1,[0 0])
            % Reliable: 2
            %------------------------------------------------------------------------------
            
            X = X*pi/180;
            Y = Y*pi/180;
            
            Rho   = sqrt(X.^2+Y.^2);
            C     = atan(Rho);
            
            Long1 = CenterVec(1)*pi/180;
            Lat1  = CenterVec(2)*pi/180;
            
            if (Rho==0)
                Lat = Lat1;
                Long = Long1;
            else
                Lat   = asin(cos(C).*sin(Lat1) + Y.*sin(C).*cos(Lat1)./Rho);
                Long  = Long1 + atan(X.*sin(C)./(Rho.*cos(Lat1).*cos(C) - Y.*sin(Lat1).*sin(C)));
            end
            
            Lat = Lat*180/pi;
            Long = Long*180/pi;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
    properties(Dependent=true, Hidden=true) % bullshit properties so the header can be read off from txt/fits files
        
        CD1_1;
        CD1_2;
        CD2_1;
        CD2_2;
        
        PV1_0;
        PV1_1;
        PV1_2;
        PV1_3;
        PV1_4;
        PV1_5;
        PV1_6;
        PV1_7;
        PV1_8;
        PV1_9;
        PV1_10;
        PV1_11;
        PV1_12;
        PV1_13;
        PV1_14;
        PV1_15;
        PV1_16;
        PV1_17;
        PV1_18;
        PV1_19;
        PV1_20;
        PV1_21;
        PV1_22;
        PV1_23;
        PV1_24;
        PV1_25;
        PV1_26;
        PV1_27;
        PV1_28;
        PV1_29;
        PV1_30;
        PV1_31;
        PV1_32;
        PV1_33;
        PV1_34;
        PV1_35;
        PV1_36;
        PV1_37;
        PV1_38;
        PV1_39;

        PV2_0;
        PV2_1;
        PV2_2;
        PV2_3;
        PV2_4;
        PV2_5;
        PV2_6;
        PV2_7;
        PV2_8;
        PV2_9;
        PV2_10;
        PV2_11;
        PV2_12;
        PV2_13;
        PV2_14;
        PV2_15;
        PV2_16;
        PV2_17;
        PV2_18;
        PV2_19;
        PV2_20;
        PV2_21;
        PV2_22;
        PV2_23;
        PV2_24;
        PV2_25;
        PV2_26;
        PV2_27;
        PV2_28;
        PV2_29;
        PV2_30;
        PV2_31;
        PV2_32;
        PV2_33;
        PV2_34;
        PV2_35;
        PV2_36;
        PV2_37;
        PV2_38;
        PV2_39;
        
    end
    
    methods % getters for bullshit dependent properties
        
        function val = get.CD1_1(obj), if ~isempty(obj.CD), val = obj.CD(1,1); else, val = []; end; end
        function val = get.CD1_2(obj), if size(obj.CD,1)>1, val = obj.CD(2,1); else, val = []; end; end
        function val = get.CD2_1(obj), if size(obj.CD,2)>1, val = obj.CD(1,2); else, val = []; end; end
        function val = get.CD2_2(obj), if numel(obj.CD)>3, val = obj.CD(2,2); else, val = []; end; end
        
        function val = get.PV1_0(obj), if size(obj.PV,1)>0, val = obj.PV(1,1); else, val = []; end; end
        function val = get.PV1_1(obj), if size(obj.PV,1)>1, val = obj.PV(2,1); else, val = []; end; end
        function val = get.PV1_2(obj), if size(obj.PV,1)>2, val = obj.PV(3,1); else, val = []; end; end
        function val = get.PV1_3(obj), if size(obj.PV,1)>3, val = obj.PV(4,1); else, val = []; end; end
        function val = get.PV1_4(obj), if size(obj.PV,1)>4, val = obj.PV(5,1); else, val = []; end; end
        function val = get.PV1_5(obj), if size(obj.PV,1)>5, val = obj.PV(6,1); else, val = []; end; end
        function val = get.PV1_6(obj), if size(obj.PV,1)>6, val = obj.PV(7,1); else, val = []; end; end
        function val = get.PV1_7(obj), if size(obj.PV,1)>7, val = obj.PV(8,1); else, val = []; end; end
        function val = get.PV1_8(obj), if size(obj.PV,1)>8, val = obj.PV(9,1); else, val = []; end; end
        function val = get.PV1_9(obj), if size(obj.PV,1)>9, val = obj.PV(10,1); else, val = []; end; end
        function val = get.PV1_10(obj), if size(obj.PV,1)>10, val = obj.PV(11,1); else, val = []; end; end
        function val = get.PV1_11(obj), if size(obj.PV,1)>11, val = obj.PV(12,1); else, val = []; end; end
        function val = get.PV1_12(obj), if size(obj.PV,1)>12, val = obj.PV(13,1); else, val = []; end; end
        function val = get.PV1_13(obj), if size(obj.PV,1)>13, val = obj.PV(14,1); else, val = []; end; end
        function val = get.PV1_14(obj), if size(obj.PV,1)>14, val = obj.PV(15,1); else, val = []; end; end
        function val = get.PV1_15(obj), if size(obj.PV,1)>15, val = obj.PV(16,1); else, val = []; end; end
        function val = get.PV1_16(obj), if size(obj.PV,1)>16, val = obj.PV(17,1); else, val = []; end; end
        function val = get.PV1_17(obj), if size(obj.PV,1)>17, val = obj.PV(18,1); else, val = []; end; end
        function val = get.PV1_18(obj), if size(obj.PV,1)>18, val = obj.PV(19,1); else, val = []; end; end
        function val = get.PV1_19(obj), if size(obj.PV,1)>19, val = obj.PV(20,1); else, val = []; end; end
        function val = get.PV1_20(obj), if size(obj.PV,1)>20, val = obj.PV(21,1); else, val = []; end; end
        function val = get.PV1_21(obj), if size(obj.PV,1)>21, val = obj.PV(22,1); else, val = []; end; end
        function val = get.PV1_22(obj), if size(obj.PV,1)>22, val = obj.PV(23,1); else, val = []; end; end
        function val = get.PV1_23(obj), if size(obj.PV,1)>23, val = obj.PV(24,1); else, val = []; end; end
        function val = get.PV1_24(obj), if size(obj.PV,1)>24, val = obj.PV(25,1); else, val = []; end; end
        function val = get.PV1_25(obj), if size(obj.PV,1)>25, val = obj.PV(26,1); else, val = []; end; end
        function val = get.PV1_26(obj), if size(obj.PV,1)>26, val = obj.PV(27,1); else, val = []; end; end
        function val = get.PV1_27(obj), if size(obj.PV,1)>27, val = obj.PV(28,1); else, val = []; end; end
        function val = get.PV1_28(obj), if size(obj.PV,1)>28, val = obj.PV(29,1); else, val = []; end; end
        function val = get.PV1_29(obj), if size(obj.PV,1)>29, val = obj.PV(30,1); else, val = []; end; end
        function val = get.PV1_30(obj), if size(obj.PV,1)>30, val = obj.PV(31,1); else, val = []; end; end
        function val = get.PV1_31(obj), if size(obj.PV,1)>31, val = obj.PV(32,1); else, val = []; end; end
        function val = get.PV1_32(obj), if size(obj.PV,1)>32, val = obj.PV(33,1); else, val = []; end; end
        function val = get.PV1_33(obj), if size(obj.PV,1)>33, val = obj.PV(34,1); else, val = []; end; end
        function val = get.PV1_34(obj), if size(obj.PV,1)>34, val = obj.PV(35,1); else, val = []; end; end
        function val = get.PV1_35(obj), if size(obj.PV,1)>35, val = obj.PV(36,1); else, val = []; end; end
        function val = get.PV1_36(obj), if size(obj.PV,1)>36, val = obj.PV(37,1); else, val = []; end; end
        function val = get.PV1_37(obj), if size(obj.PV,1)>37, val = obj.PV(38,1); else, val = []; end; end
        function val = get.PV1_38(obj), if size(obj.PV,1)>38, val = obj.PV(39,1); else, val = []; end; end
        function val = get.PV1_39(obj), if size(obj.PV,1)>39, val = obj.PV(40,1); else, val = []; end; end

        function val = get.PV2_0(obj), if size(obj.PV,1)>0, val = obj.PV(1,2); else, val = []; end; end
        function val = get.PV2_1(obj), if size(obj.PV,1)>1, val = obj.PV(2,2); else, val = []; end; end
        function val = get.PV2_2(obj), if size(obj.PV,1)>2, val = obj.PV(3,2); else, val = []; end; end
        function val = get.PV2_3(obj), if size(obj.PV,1)>3, val = obj.PV(4,2); else, val = []; end; end
        function val = get.PV2_4(obj), if size(obj.PV,1)>4, val = obj.PV(5,2); else, val = []; end; end
        function val = get.PV2_5(obj), if size(obj.PV,1)>5, val = obj.PV(6,2); else, val = []; end; end
        function val = get.PV2_6(obj), if size(obj.PV,1)>6, val = obj.PV(7,2); else, val = []; end; end
        function val = get.PV2_7(obj), if size(obj.PV,1)>7, val = obj.PV(8,2); else, val = []; end; end
        function val = get.PV2_8(obj), if size(obj.PV,1)>8, val = obj.PV(9,2); else, val = []; end; end
        function val = get.PV2_9(obj), if size(obj.PV,1)>9, val = obj.PV(10,2); else, val = []; end; end
        function val = get.PV2_10(obj), if size(obj.PV,1)>10, val = obj.PV(11,2); else, val = []; end; end
        function val = get.PV2_11(obj), if size(obj.PV,1)>11, val = obj.PV(12,2); else, val = []; end; end
        function val = get.PV2_12(obj), if size(obj.PV,1)>12, val = obj.PV(13,2); else, val = []; end; end
        function val = get.PV2_13(obj), if size(obj.PV,1)>13, val = obj.PV(14,2); else, val = []; end; end
        function val = get.PV2_14(obj), if size(obj.PV,1)>14, val = obj.PV(15,2); else, val = []; end; end
        function val = get.PV2_15(obj), if size(obj.PV,1)>15, val = obj.PV(16,2); else, val = []; end; end
        function val = get.PV2_16(obj), if size(obj.PV,1)>16, val = obj.PV(17,2); else, val = []; end; end
        function val = get.PV2_17(obj), if size(obj.PV,1)>17, val = obj.PV(18,2); else, val = []; end; end
        function val = get.PV2_18(obj), if size(obj.PV,1)>18, val = obj.PV(19,2); else, val = []; end; end
        function val = get.PV2_19(obj), if size(obj.PV,1)>19, val = obj.PV(20,2); else, val = []; end; end
        function val = get.PV2_20(obj), if size(obj.PV,1)>20, val = obj.PV(21,2); else, val = []; end; end
        function val = get.PV2_21(obj), if size(obj.PV,1)>21, val = obj.PV(22,2); else, val = []; end; end
        function val = get.PV2_22(obj), if size(obj.PV,1)>22, val = obj.PV(23,2); else, val = []; end; end
        function val = get.PV2_23(obj), if size(obj.PV,1)>23, val = obj.PV(24,2); else, val = []; end; end
        function val = get.PV2_24(obj), if size(obj.PV,1)>24, val = obj.PV(25,2); else, val = []; end; end
        function val = get.PV2_25(obj), if size(obj.PV,1)>25, val = obj.PV(26,2); else, val = []; end; end
        function val = get.PV2_26(obj), if size(obj.PV,1)>26, val = obj.PV(27,2); else, val = []; end; end
        function val = get.PV2_27(obj), if size(obj.PV,1)>27, val = obj.PV(28,2); else, val = []; end; end
        function val = get.PV2_28(obj), if size(obj.PV,1)>28, val = obj.PV(29,2); else, val = []; end; end
        function val = get.PV2_29(obj), if size(obj.PV,1)>29, val = obj.PV(30,2); else, val = []; end; end
        function val = get.PV2_30(obj), if size(obj.PV,1)>30, val = obj.PV(31,2); else, val = []; end; end
        function val = get.PV2_31(obj), if size(obj.PV,1)>31, val = obj.PV(32,2); else, val = []; end; end
        function val = get.PV2_32(obj), if size(obj.PV,1)>32, val = obj.PV(33,2); else, val = []; end; end
        function val = get.PV2_33(obj), if size(obj.PV,1)>33, val = obj.PV(34,2); else, val = []; end; end
        function val = get.PV2_34(obj), if size(obj.PV,1)>34, val = obj.PV(35,2); else, val = []; end; end
        function val = get.PV2_35(obj), if size(obj.PV,1)>35, val = obj.PV(36,2); else, val = []; end; end
        function val = get.PV2_36(obj), if size(obj.PV,1)>36, val = obj.PV(37,2); else, val = []; end; end
        function val = get.PV2_37(obj), if size(obj.PV,1)>37, val = obj.PV(38,2); else, val = []; end; end
        function val = get.PV2_38(obj), if size(obj.PV,1)>38, val = obj.PV(39,2); else, val = []; end; end
        function val = get.PV2_39(obj), if size(obj.PV,1)>39, val = obj.PV(40,2); else, val = []; end; end
                     
    end
    
end

