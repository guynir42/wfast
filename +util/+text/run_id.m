function val = run_id(folder)
% Usage: val = run_id(folder)
% Extract the run identifier from the given folder. 
% Finds the last place in the path where a date is followed by anohter 
% sub-folder (e.g., /path/to/data/2020-06-01/ecliptic_run1 will only get
% the last two folders). 
% This is useful for getting the path of a file in other computers where 
% the data folder is different but the unique run identifier is the same.  

    if nargin==0, help('util.text.run_id'); return; end
    
    idx = regexp(folder, '(\d{4}-\d{2}-\d{2}[\\/][a-zA-Z]+)');
    
    if isempty(idx)
        val = '';
    else
        
         val = regexp(folder(idx(end):end), '\d{4}-\d{2}-\d{2}[\\/].+', 'match'); 
         if iscell(val)
             val = val{1};
         end
         
    end

end