%% this script scans the header / text files and gets the telescope altitude vs. focus point for all runs in a given time range

start_date = datetime('2020-06-15');
end_date = datetime('2020-07-07'); 

d = util.sys.WorkingDirectory('E:\data_backup'); 

list = d.dir; 

data = []; 

for ii = 1:length(list)
    
    t = datetime(list{ii}(1:10)); 
    
    if t>=start_date && t<=end_date
        
%         disp(list{ii});

        d.cd(list{ii}); 
        
        folders = d.dir;
        
        for jj = 1:length(folders)
            
            if util.text.cs(folders{jj}, 'dark', 'flat')
                continue;
            end
            
            filename = fullfile(d.pwd, folders{jj}, 'A_README.txt'); 
            
            fid = fopen(filename); 
            
            if fid<0
                continue;
            end
            
            in_focus = 0; % this tells you when you are reading the focuser data
            st = struct;
            st.filename = filename;
            st.pos = NaN;
            st.tip = NaN;
            st.tilt = NaN;
            st.Alt = NaN;
            st.RA = NaN;
            st.Dec = NaN;
            st.LST = NaN;
            
            for kk = 1:1e4
                
                line = fgetl(fid); 
                if isnumeric(line)
                    break;
                end
                
                line = strtrim(line);
                
                if isempty(line)
                    continue;
                end
                
                if length(line)>=25 && strcmp(line(1:25), 'OBJECT LOG FOR: "focuser"')
                    in_focus = 1;
                elseif length(line)>=25 && strcmp(line(1:15), 'OBJECT LOG FOR:')
                    in_focus = 0;
                end
                
                if in_focus
                    [~,idx] = regexp(line, '^pos:\s+'); 

                    if ~isempty(idx)
                        st.pos = util.text.parse_value(line(idx+1:end)); 
                    end
                    
                    [~,idx] = regexp(line, '^tip:\s+'); 

                    if ~isempty(idx)
                        st.tip = util.text.parse_value(line(idx+1:end)); 
                    end
                    
                    [~,idx] = regexp(line, '^tilt:\s+'); 

                    if ~isempty(idx)
                        st.tilt = util.text.parse_value(line(idx+1:end)); 
                    end
                    
                end
                
                [~,idx] = regexp(line, '^ALT:\s+'); 

                if ~isempty(idx) && isnan(st.Alt)
                    st.Alt = util.text.parse_value(line(idx+1:end));
                    if isempty(st.Alt) || ischar(st.Alt), st.Alt = NaN; end
                end
                
                [~,idx] = regexp(line, '^LST_deg:\s+'); 

                if ~isempty(idx)
                    st.LST = util.text.parse_value(line(idx+1:end));
                    if isempty(st.LST) || ischar(st.LST), st.LST = NaN; end
                end
                
                [~,idx] = regexpi(line, '^TELRA_DEG:\s+'); 

                if ~isempty(idx)
                    st.RA = util.text.parse_value(line(idx+1:end));
                    if isempty(st.RA) || ischar(st.RA), st.RA = NaN; end
                end
                
                [~,idx] = regexpi(line, '^TELDEC_DEG:\s+'); 

                if ~isempty(idx)
                    st.Dec = util.text.parse_value(line(idx+1:end));
                    if isempty(st.Dec) || ischar(st.Dec), st.Dec = NaN; end
                end
                
                if isempty(st.RA) || isnan(st.RA)
                    [~,idx] = regexp(line, '^RA_DEG:\s+'); 

                    if ~isempty(idx)
                        st.RA = util.text.parse_value(line(idx+1:end));
                        if isempty(st.RA) || ischar(st.RA), st.RA = NaN; end
                    end
                end
                
                if isempty(st.Dec) || isnan(st.Dec)
                    [~,idx] = regexp(line, '^DEC_DEG:\s+'); 

                    if ~isempty(idx)
                        st.Dec = util.text.parse_value(line(idx+1:end));
                        if isempty(st.Dec) || ischar(st.Dec), st.Dec = NaN; end
                    end
                end
            end
            
            if isempty(data)
                data = st;
            else
                data(end+1) = st;
            end
            
            fclose(fid); 
            
        end
        
        d.up;
        
    end
    
end

%% plot the results

east_indices = find([data.LST]>[data.RA]-5); 
west_indices = find([data.LST]<[data.RA]-5); 

f1 = util.plot.FigHandler('altitude'); 
f1.clear;

ax = axes('Parent', f1.fig); 

plot(ax, [data(east_indices).Alt], [data(east_indices).pos], 'pb', 'MarkerSize', 15); 

hold(ax, 'on'); 

plot(ax, [data(west_indices).Alt], [data(west_indices).pos], 'pr', 'MarkerSize', 15); 

fr = util.fit.polyfit([data.Alt], [data.pos], 'order', 2); 

plot(ax, sort(fr.x), fr.func(sort(fr.x)), '-', 'LineWidth', 3);

util.plot.inner_title(sprintf('pos= %4.2f + %7.5f*A - %12.10f*A^2', fr.coeffs(1), fr.coeffs(2), abs(fr.coeffs(3))), 'ax', ax, 'Position', 'bottom', 'margin', 0.2); 

hold(ax, 'off'); 

box(ax, 'on'); 

xlabel(ax, 'Altitude above horizon [degrees]'); 
ylabel(ax, 'Focuser position [mm]'); 
ax.FontSize = 20;

%% save the plot

util.sys.print(fullfile(getenv('WFAST'), '/scripts/plots/focus_altitude')); 










