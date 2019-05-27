classdef Filter < handle
   
    properties
        
        name = '';
        
        wavelength; % nm
        bandwidth; % nm
        zero_point; % ergs/s/cm^2/Hz
        
        lambda; % wavelength axis for transmission curve
        trans; % transmission curve (fractions of 1) for each wavelength
        
    end
    
    properties(Hidden=true)
        
        version = 1.00; 

    end
    
    methods % constructor
        
        function obj = Filter(varargin)
            
            import util.text.cs;
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Filter')
                
                obj.name = varargin{1}.name;
                obj.wavelength = varargin{1}.wavelength;
                obj.bandwidth = varargin{1}.bandwidth;
                obj.zero_point = varargin{1}.zero_point;
                obj.lambda = varargin{1}.lambda;
                obj.trans = varargin{1}.trans;
                                
                return;
                
            end
            
            if isempty(varargin) || isempty(varargin{1})
                filter_type = 'V';
            elseif ischar(varargin{1})
                filter_type = varargin{1};
            else
                filter_type = '';
                % do something smart with the inputs...
            end
            
            obj.name = filter_type;
            
            if cs(filter_type, 'clear')
                obj.wavelength = 550;
                obj.bandwidth = 300;
                % should we have different zero points for each filter??
            elseif cs(filter_type, 'U')
                obj.wavelength = 365;
                obj.bandwidth = 66;
                obj.name = 'U';
            elseif cs(filter_type, 'B')
                obj.wavelength = 445;
                obj.bandwidth = 94;
                obj.name = 'B';
            elseif cs(filter_type, 'V')
                obj.wavelength = 551;
                obj.bandwidth = 88;
                obj.name = 'V';
            elseif cs(filter_type, 'R')
                obj.wavelength = 658;
                obj.bandwidth = 138;
                obj.name = 'R';
            elseif cs(filter_type, 'I')
                obj.wavelength = 806;
                obj.bandwidth = 149;
                obj.name = 'I';
            elseif cs(filter_type, 'K')
                obj.wavelength = 2190;
                obj.bandwidth = 390;
                obj.name = 'K';
            elseif cs(filter_type, 'luminance')
                obj.wavelength = 550;
                obj.bandwidth = 300;
                obj.name = 'luminance';
                % add more filters! 
            else
                error(['Unknown filter: ' filter_type ' try U,B,V, etc.']);
            end
            
        end
        
    end
    
    methods % getters
        
        function F = getFlux(obj, mag)
        % usage: getFlux(mag)
            if nargin==0
                help('head.Filter.getFlux')
                return;
            end
            
            if isempty(obj.zero_point)
                zp = 3.64e-20; % ergs/cm^2/sec/Hz
            else
                zp = obj.zero_point;
            end
            
            bandwidth_Hz = (1./(obj.wavelength-obj.bandwidth/2) - 1./(obj.wavelength+obj.bandwidth/2))*3e17; % speed of light in nm/second
            f0 = zp*bandwidth_Hz; % zero flux. ergs/cm^2/sec
                        
            F = f0.*10.^(-0.4*mag);
 
        end
        
        function mag = getMagFromFlux(obj, flux)
        % usage: getMagFromFlux(flux_erg_cm2_sec)

            if nargin==0
                help('head.Filter.getMagFromFlux');
                return;
            end
            
            if isempty(obj.zero_point)
                zp = 3.64e-20;
            else
                zp = obj.zero_point;
            end
            
            bandwidth_Hz = (1./(obj.wavelength-obj.bandwidth/2) - 1./(obj.wavelength+obj.bandwidth/2))*3e17; % speed of light in nm/second
            f0 = zp*bandwidth_Hz; % zero flux. ergs/cm^2/sec
                        
            mag = -2.5*log10(flux/f0);
            
        end
        
        function C = getCount(obj, mag, area_cm2, expT)
            % usage: getCount(mag, area_cm2=1, expT=1)
            
            if nargin==0
                help('head.Filter.getCount');
                return;
            end
            
            if nargin<3 || isempty(area_cm2)
                area_cm2 = 1;
            end
            
            if nargin<4 || isempty(expT)
                expT = 1;
            end
            
            flux = obj.getFlux(mag);
            
            h = 6.626e-27; % erg*sec
            c = 3e17; % nm/sec
        
            C = flux./(h.*c./obj.wavelength)*area_cm2*expT;
                
        end
        
        function mag = getMagFromCount(obj, count, area_cm2, expT)
            % usage: getMagFromCount(count, area_cm2=1, expT=1)
            
            if nargin==0
                help('head.Filter.getMagFromCount');
                return;
            end
            
            if nargin<3 || isempty(area_cm2)
                area_cm2 = 1;
            end
            
            if nargin<4 || isempty(expT)
                expT = 1;
            end
                        
            h = 6.626e-27; % erg*sec
            c = 3e17; % nm/sec
         
            flux = count*h*c/obj.wavelength/area_cm2/expT;
            
            mag = obj.getMagFromFlux(flux);
            
        end
        
        function disp(obj)
            
            if isempty(obj)
                disp('empty filter');
                return;
            end
            
            if isempty(obj.name)
                str = '';
            else
                str = ['"' obj.name '" '];
            end
            
            str = [str '(' num2str(obj.wavelength) 'nm / ' num2str(obj.bandwidth) 'nm)'];
            
            disp(str);
            
        end
        
        function val = double(obj)
           
            val = obj.wavelength;
            
        end
        
    end
    
end