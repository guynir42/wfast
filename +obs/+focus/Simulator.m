classdef Simulator < handle
    
    properties(Transient=true)
        
        gui@obs.focus.gui.SimGUI;
        
    end
    
    properties
        
        pos;
        best_pos = 3000;
        
        defocus_parameter = 0.01;
        inner_annulus = 0.25;
        
    end
    
    methods
        
        function obj = Simulator(varargin)
            
            obj.pos = round(normrnd(obj.best_pos, 1./obj.defocus_parameter));
            
        end
        
        function psf = makePSF(obj)
                
            df = obj.pos-obj.best_pos;
            
            if df==0
                psf = 1;
            else

                r2 = (obj.defocus_parameter.*df).^2;
                r1 = r2.*obj.inner_annulus;

                psf = util.img.annulusMask(2*r2, 'r_min', r1,'r_max', r2);

                if util.stat.sum2(psf)==0
                    psf = 1;
                else
                    psf = psf./util.stat.sum2(psf);
                end
                
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.focus.gui.SimGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
end