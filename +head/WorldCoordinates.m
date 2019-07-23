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
        CTYPE = {'RA---TPV'  'DEC--TPV'};
        CUNIT = {'deg'  'deg'};
        CRPIX = []; % pixel position of transformation anchor 
        CRVAL = []; % RA/Dec position of transformation anchor 
        CDELT = []; % shift between anchors
        CD = []; % rotation matrix
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
                if obj.debug_bit, fprintf('WorldCoordinates copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isstruct(varargin{1}) && isfield(varargin{1}, 'CRVAL') % can add more checks
                if obj.debug_bit, fprintf('WorldCoordinates struct based constructor v%4.2f\n', obj.version); end
                obj.input(varargin{1});
            elseif ~isempty(varargin) && (isa(varargin{1}, 'SIM') || isa(varargin{1},'ClassWCS'))
                if obj.debug_bit, fprintf('WorldCoordinates object based constructor v%4.2f\n', obj.version); end
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
            
            if isempty(obj.PV)
                val = [];
            else
                val = acotd(obj.PV(2,1)./obj.PV(3,1));
            end
            
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
        
        function RA_Dec = xy2coo(obj, x, y)
            
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
            Xout = vec(1); % units of degrees
            Yout = vec(2); % units of degrees
            
            if obj.use_tpv
            
                if all(~cellfun(@isempty, strfind(obj.CTYPE, 'TPV')))

                    R  = sqrt(Xout.^2 + Yout.^2); % units of degrees
            
                    [Xpowers, Ypowers] = head.WorldCoordinates.tpv_polydef;
            
                    % these lines assume PV is a 2x40 coefficients matrix with
                    % NaNs or zeros where there is no contribution to that polynomial
                    Xout = sum(obj.PV(:,1).*Xout.^Xpowers(:,1).*Yout.^Xpowers(:,2).*R.^Xpowers(:,3), 1, 'omitnan');
                    Yout = sum(obj.PV(:,2).*Xout.^Ypowers(:,1).*Yout.^Ypowers(:,2).*R.^Ypowers(:,3), 1, 'omitnan');

                else
                    error('Unknown CTYPE "%s". Use tpv instead...', obj.CTYPE{1});
                end

            end
            
            [RA, Dec] = head.WorldCoordinates.pr_ignomonic(Xout, Yout, obj.CRVAL);

            RA_Dec(1) = RA;
            RA_Dec(2) = Dec;

        end
        
        function XY = coo2xy(obj, RA, Dec)
            
            error('Not yet implemented!');
            
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
        PV1_40;

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
        PV2_40;
        
    end
    
    methods % getters for bullshit dependent properties
        
        function val = get.CD1_1(obj), val = obj.CD(1,1); end
        function val = get.CD1_2(obj), val = obj.CD(2,1); end
        function val = get.CD2_1(obj), val = obj.CD(1,2); end
        function val = get.CD2_2(obj), val = obj.CD(2,2); end
        
        
        function val = get.PV1_0(obj), val = obj.PV(1,1); end
        function val = get.PV1_1(obj), val = obj.PV(2,1); end
        function val = get.PV1_2(obj), val = obj.PV(3,1); end
        function val = get.PV1_3(obj), val = obj.PV(4,1); end
        function val = get.PV1_4(obj), val = obj.PV(5,1); end
        function val = get.PV1_5(obj), val = obj.PV(6,1); end
        function val = get.PV1_6(obj), val = obj.PV(7,1); end
        function val = get.PV1_7(obj), val = obj.PV(8,1); end
        function val = get.PV1_8(obj), val = obj.PV(9,1); end
        function val = get.PV1_9(obj), val = obj.PV(10,1); end
        function val = get.PV1_10(obj), val = obj.PV(11,1); end
        function val = get.PV1_11(obj), val = obj.PV(12,1); end
        function val = get.PV1_12(obj), val = obj.PV(13,1); end
        function val = get.PV1_13(obj), val = obj.PV(14,1); end
        function val = get.PV1_14(obj), val = obj.PV(15,1); end
        function val = get.PV1_15(obj), val = obj.PV(16,1); end
        function val = get.PV1_16(obj), val = obj.PV(17,1); end
        function val = get.PV1_17(obj), val = obj.PV(18,1); end
        function val = get.PV1_18(obj), val = obj.PV(19,1); end
        function val = get.PV1_19(obj), val = obj.PV(20,1); end
        function val = get.PV1_20(obj), val = obj.PV(21,1); end
        function val = get.PV1_21(obj), val = obj.PV(22,1); end
        function val = get.PV1_22(obj), val = obj.PV(23,1); end
        function val = get.PV1_23(obj), val = obj.PV(24,1); end
        function val = get.PV1_24(obj), val = obj.PV(25,1); end
        function val = get.PV1_25(obj), val = obj.PV(26,1); end
        function val = get.PV1_26(obj), val = obj.PV(27,1); end
        function val = get.PV1_27(obj), val = obj.PV(28,1); end
        function val = get.PV1_28(obj), val = obj.PV(29,1); end
        function val = get.PV1_29(obj), val = obj.PV(30,1); end
        function val = get.PV1_30(obj), val = obj.PV(31,1); end
        function val = get.PV1_31(obj), val = obj.PV(32,1); end
        function val = get.PV1_32(obj), val = obj.PV(33,1); end
        function val = get.PV1_33(obj), val = obj.PV(34,1); end
        function val = get.PV1_34(obj), val = obj.PV(35,1); end
        function val = get.PV1_35(obj), val = obj.PV(36,1); end
        function val = get.PV1_36(obj), val = obj.PV(37,1); end
        function val = get.PV1_37(obj), val = obj.PV(38,1); end
        function val = get.PV1_38(obj), val = obj.PV(39,1); end
        function val = get.PV1_39(obj), val = obj.PV(40,1); end
        function val = get.PV1_40(obj), val = obj.PV(40,1); end

        function val = get.PV2_0(obj), val = obj.PV(1,2); end
        function val = get.PV2_1(obj), val = obj.PV(2,2); end
        function val = get.PV2_2(obj), val = obj.PV(3,2); end
        function val = get.PV2_3(obj), val = obj.PV(4,2); end
        function val = get.PV2_4(obj), val = obj.PV(5,2); end
        function val = get.PV2_5(obj), val = obj.PV(6,2); end
        function val = get.PV2_6(obj), val = obj.PV(7,2); end
        function val = get.PV2_7(obj), val = obj.PV(8,2); end
        function val = get.PV2_8(obj), val = obj.PV(9,2); end
        function val = get.PV2_9(obj), val = obj.PV(10,2); end
        function val = get.PV2_10(obj), val = obj.PV(11,2); end
        function val = get.PV2_11(obj), val = obj.PV(12,2); end
        function val = get.PV2_12(obj), val = obj.PV(13,2); end
        function val = get.PV2_13(obj), val = obj.PV(14,2); end
        function val = get.PV2_14(obj), val = obj.PV(15,2); end
        function val = get.PV2_15(obj), val = obj.PV(16,2); end
        function val = get.PV2_16(obj), val = obj.PV(17,2); end
        function val = get.PV2_17(obj), val = obj.PV(18,2); end
        function val = get.PV2_18(obj), val = obj.PV(19,2); end
        function val = get.PV2_19(obj), val = obj.PV(20,2); end
        function val = get.PV2_20(obj), val = obj.PV(21,2); end
        function val = get.PV2_21(obj), val = obj.PV(22,2); end
        function val = get.PV2_22(obj), val = obj.PV(23,2); end
        function val = get.PV2_23(obj), val = obj.PV(24,2); end
        function val = get.PV2_24(obj), val = obj.PV(25,2); end
        function val = get.PV2_25(obj), val = obj.PV(26,2); end
        function val = get.PV2_26(obj), val = obj.PV(27,2); end
        function val = get.PV2_27(obj), val = obj.PV(28,2); end
        function val = get.PV2_28(obj), val = obj.PV(29,2); end
        function val = get.PV2_29(obj), val = obj.PV(30,2); end
        function val = get.PV2_30(obj), val = obj.PV(31,2); end
        function val = get.PV2_31(obj), val = obj.PV(32,2); end
        function val = get.PV2_32(obj), val = obj.PV(33,2); end
        function val = get.PV2_33(obj), val = obj.PV(34,2); end
        function val = get.PV2_34(obj), val = obj.PV(35,2); end
        function val = get.PV2_35(obj), val = obj.PV(36,2); end
        function val = get.PV2_36(obj), val = obj.PV(37,2); end
        function val = get.PV2_37(obj), val = obj.PV(38,2); end
        function val = get.PV2_38(obj), val = obj.PV(39,2); end
        function val = get.PV2_39(obj), val = obj.PV(40,2); end
                
    end
    
end

