function analysis_dv(sdata)

    % do the analysis over... (see 'psych' variable)
    % 1 = setsize
    % 2 = variance
    % 3 = constraint
    % 4 = duration
    vbs = [1,2,4];
    nb_vbs = length(vbs);

    % plots
    % 1 = y-scale
    % 2 = variance
    % 3 = x-mean
    toplots = [1];
    
    for toplot = toplots
        
        % figure
        figure('color',[1 1 1]);

        % change figure size
        f_xy = get(gcf,'position');
        set(gcf,'position',[f_xy(1),f_xy(2),nb_vbs*f_xy(3),f_xy(4)]);
        
        for i_vb = 1:nb_vbs
            vb = vbs(i_vb);

            % subplot
            subplot(1,nb_vbs,i_vb);
            
            % tools plot
            tools_plotdv(sdata,vb,toplot);
        end
    end
end