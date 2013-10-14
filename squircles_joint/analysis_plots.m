function analysis_plots(sdata)

    %load data
    load_data;

    % do the analysis over... (see 'psych' variable)
    % 1 = setsize
    % 2 = variance
    % 3 = constraint
    % 4 = duration
    vbs = [1,2,4];

    % plot ...
    % 1 = performance
    % 2 = error
    % 3 = absolute error
    % 4 = response
    % 5 = reaction times
    toplots = [4];

    for toplot = toplots
        % select variable to plot
        fit = 0;
        switch toplot
            case 1
                mvalues = mcor;
            case 2
                mvalues = merr;
            case 3
                mvalues = maerr;
            case 4
                if sdata.id == 4
                    mvalues = mcat;
                else
                    mvalues = .5*(mcat+1);
                end
                fit = 1;
            case 5
                mvalues = mrt;
            otherwise
                error(['analysis_plots: toplot=',num2str(toplot),'is not valid.']);
        end
        tools_psychplot(sdata,mvalues,vbs,fit,toplot);
    end
end