people = {'Avishay Gal-Yam', 'Eran Ofek', 'Sagi Ben-Ami', 'Guy Nir', 'Barak Zackay', 'Ilan Manulis', ...
    'David Polishook', 'Maayan Soumagnac', 'Noam Segev', 'Baptiste Callendret', 'Samantha Goldwasser',...
    'Betina Reyna', 'Rachel Bruch', 'Arka Ghosh', 'Jonathan Mushkin', 'Ofir Hershko', 'Oz Diner', ...
    'Ido Irani', 'Amir Sharon', 'Nora-Linn Strotjohann', 'Ofer Yaron', 'Danny Khazov', 'Tali Engel', ...
    'Nikola Knezevic', 'Sladjana Knezevic', 'Adam Rubin', 'Ira Bar', 'Ilan Sagiv',...
    'Michael Rappaport', 'Haim Sade', 'Sammy Ben-Gigi', 'Shai Kaspi', 'Yigal Shachar','Benjamin Pasmantirer'};

authors = {'Guy Nir', 'Eran Ofek', 'Avishay Gal-Yam', 'Sagi Ben-Ami', 'Noam Segev', 'David Polishook', 'Ilan Manulis', 'Ofir Hershko', 'Oz Diner', 'Barak Zackay', 'Ofer Yaron'}; 

T = table;

for ii = 1:length(people)
    
    if ~ismember(people{ii}, authors)
        
        name = strsplit(people{ii});
        
        T{end+1,:} = name; 
                
    end
    
end

T.Properties.VariableNames = {'first_name', 'last_name'}; 

T = sortrows(T, 2); % sort alphabetically

for ii = 1:height(T)
    
    fprintf('%s.~%s,\n', T{ii,1}{1}(1), T{ii,2}{1}); 
    
end

fprintf('\n\n'); 
