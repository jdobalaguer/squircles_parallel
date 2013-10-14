function [b,X,yest] = tools_logregression(y,X)
    %{
    INPUTS
        X = [vx1' ; vx2'; ... ]
        y = [y1 y2 y3 ...]'
    OUTPUT
        b = [b1 b2 b3 ...]'
    %}

    % check inputs
    if size(X,1)~=length(y)
        error('bad inputs');
    end
    
    % VARIABLES
    % step
    mu = .1;
    % weights
    b = rand(size(X,2),1);
    % error
    e = [];

    % ITERATION
    steps = 1;
    me1 = inf;
    looksteps = 10000;
    while steps<(10*looksteps) || me1<me2
        steps = steps+1;
        % estimation
        g = exp(X*b);
        yest = g./(g+1);
        % error
        e_last = mean(power(y - yest,2));
        %   windowed means
        if steps > 2*looksteps+1
            me2 = me1;
            me1 = mean(e((end-looksteps+1):end));
            me2 = mean(e((end-2*looksteps+1):(end-looksteps+1)));
        %   print error periodically
            if ~mod(steps,looksteps)
                fprintf(['logregression: error is ',num2str(me1),'\n']);    
            end
        end
        e(end+1) = e_last;
        % update from a random sample
        i_y = randi(length(y));
        dbx = X(i_y,:)*g(i_y)/power(1+g(i_y),2);
        dby = y(i_y) - yest(i_y);
        b = b + mu.*(dbx.*dby)';
    end
    fprintf(['logregression: error is ',num2str(e_last),'\n']);
    % PLOT ERROR
    figure;
    plot(e);
    drawnow;
end