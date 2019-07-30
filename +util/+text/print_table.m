function txt_out = print_table(T, varargin)

    if nargin==0, help('util.text.print_table'); return; end
    
    input = util.text.InputVars;
    input.input_var('indent', 0);
    input.input_var('delim', '|');
    input.input_var('num_rows', []); 
    input.scan_vars(varargin{:});
    
    txt = '';
    
    if isa(T,'table')
    
        if ~isempty(input.num_rows) && input.num_rows<height(T)
            C = table2cell(T(1:input.num_rows,:));
            ellipses = 1;
        else
            C = table2cell(T);
            ellipses = 0;
        end
        
       
        
        C = [T.Properties.VariableNames; C];
        
        C = format_cell(C);
        
        new_line = repmat(' ', [1 input.indent]);
            
        for jj = 1:size(C,2)

            new_line = sprintf('%s %s %s', new_line, C{1,jj}, input.delim);

        end

        txt = sprintf('%s\n%s', txt, new_line);
        
        txt = sprintf('%s\n%s%s', txt, repmat(' ', [1 input.indent]), repmat('-', [1 length(new_line)-input.indent]));
        
        for ii = 2:size(C,1)
            
            new_line = repmat(' ', [1 input.indent]);
            
            for jj = 1:size(C,2)
                
                new_line = sprintf('%s %s %s', new_line, C{ii,jj}, input.delim);
                
            end
            
            txt = sprintf('%s\n%s', txt, new_line);
        
        end
        
        if ellipses
            
            new_line = repmat(' ', [1 input.indent]);
            
            for ii = 1:size(C,2)
                new_line = sprintf(['%s %' num2str(length(C{end,ii})) 's %s'], new_line, '... ', repmat(' ', [1, length(input.delim)]));
            end
            
            txt = sprintf('%s\n%s', txt, new_line);
            
        else
            txt = sprintf('%s\n%s%s', txt, repmat(' ', [1 input.indent]), repmat('-', [1 length(new_line)-input.indent]));
        end
        
        txt = sprintf('%s\n\n', txt);
        
    else
        error('Unsupprted class "%s" for first input. Try inputting a table...', class(T));
    end
    
    if nargout>0
        txt_out = txt;
    end

end

function C = format_cell(C)

    % first turn all cells into strings
    for ii = 1:size(C,1)
        
        for jj = 1:size(C,2)
            
            val = C{ii,jj};
            
            if isnumeric(val) && isscalar(val)
                C{ii,jj} = num2str(val);
            elseif ischar(val)
                C{ii,jj} = sprintf('''%s''', val);
            elseif isa(val, 'datetime')
                C{ii,jj} = util.text.time2str(val);
            else
                C{ii,jj} = sprintf('[%s %s]', util.text.print_vec(size(val)), class(val));
            end
            
        end
        
    end
    
    % now find the widest string in each column
     
    widths = max(max(cellfun(@length, C)),4); % minimal width is 4 characters anyway
    
    for ii = 1:size(C,1)
        
        for jj = 1:size(C,2)
            C{ii,jj} = sprintf(['%' num2str(widths(jj)) 's'], C{ii,jj});
        end
        
    end
    

end
