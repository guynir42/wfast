classdef VirtualSensor < dynamicprops 
    
    properties
        
        owner@obs.SensorChecker; 
        
        status = 1;
        use_this = 1; % is true get the data from this sensor. If not, ignore it and move on
        id = '';
        
        data_path = {}; 
        subsref_struct;
        data_struct = []; 
        
        jd = NaN; % latest date when data was successfully updated
        
        data_names = {'skyAmbientTemp', 'ambientTemp', 'windSpeed', 'windDir', 'humidity', 'rainFlag'...
            'outsideTemp', 'outsideHumidity', 'barometer', 'rainRate'};
        
    end
    
    methods
        
        function connect(obj)
            
            obj.status = 0;
            
            try

                if ~isempty(obj.owner) && ~isempty(obj.owner.wise_data_struct)

                    obj.subsref_struct = struct;

                    for ii = 1:length(obj.data_path)

                        if isnumeric(obj.data_path{ii})
                            obj.subsref_struct(ii).type = '()';
                            obj.subsref_struct(ii).subs = obj.data_path(ii); 
                        elseif ischar(obj.data_path{ii})
                            obj.subsref_struct(ii).type = '.';
                            obj.subsref_struct(ii).subs = obj.data_path{ii}; 
                        end
                        
                    end
                    
                    obj.data_struct = subsref(obj.owner.wise_data_struct, obj.subsref_struct);
                    
                    if ~isempty(obj.data_struct) % if we succeed in getting data
                    
                        for ii = 1:length(obj.data_names)
                            
                            new_prop = '';
                            
                            if isfield(obj.data_struct, obj.data_names{ii})
                                new_prop = obj.data_names{ii};
                            elseif isfield(obj.data_struct, obj.owner.capital(obj.data_names{ii}))
                                new_prop = obj.owner.capital(obj.data_names{ii});
                            end
                            
                            if ~isempty(new_prop)
                                
                                if ~isprop(obj, new_prop)
                                    addprop(obj, new_prop); 
                                end

                                obj.(new_prop) = obj.data_struct.(new_prop); 
                                
                            end
                            
                        end
                        
                        obj.status = 1;
                        
                    end
                    
                end

            catch ME
                fprintf('Cannot connect to virtual sensor %d\n', obj.id); 
                warning(ME.getReport); 
            end
            
        end
        
        function update(obj)
            
            import util.text.cs;
            
            obj.status = 0;
            
            obj.data_struct = subsref(obj.owner.wise_data_struct, obj.subsref_struct);
            
            if ~isempty(obj.data_struct)
                                
                for ii = 1:length(obj.data_names)
                            
                    field_name = '';

                    if isfield(obj.data_struct, obj.data_names{ii})
                        field_name = obj.data_names{ii};
                    elseif isfield(obj.data_struct, obj.owner.capital(obj.data_names{ii}))
                        field_name = obj.owner.capital(obj.data_names{ii});
                    end

                    if ~isempty(field_name)

                        obj.(field_name) = obj.data_struct.(field_name); 

                        % quality checks
                        if cs(field_name, 'clouds', 'sky ambient temperature')
                            if obj.(field_name)<-100
                                obj.(field_name) = NaN;
                            end
                        end
                        
                    end
                    
                end

                obj.status = 1;
                
            end
            
        end
        
    end
    
end