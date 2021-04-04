function show_efficiency(ev, ax)
    
    if nargin<2 || isempty(ax)
        ax = gca;
    end
    
    w = 0.1; 
    
    r_passed = [ev([ev.passed]).r]; 
    r_total = [ev.r]; 
    
    [N_total,E] = histcounts(r_total, 'BinWidth', w);     
    [N_passed] = histcounts(r_passed, 'BinEdges', E); 
    
    h_total = histogram(ax, r_total, 'BinEdges', E); 
    hold(ax, 'on'); 
    h_total.FaceColor = 'blue'; 
    h_passed = histogram(ax, r_passed, 'BinEdges', E); 
    hold(ax, 'off'); 
    h_passed.FaceColor = 'red'; 
        
    yyaxis(ax, 'right');
    
    plot(ax, E(1:end-1)+w/2, N_passed./N_total.*100); 
    
    ytickformat('%d%%'); 
    
    yyaxis(ax, 'left');
    
    ax.FontSize = 16; 
    ax.YScale = 'log'; 
    
end