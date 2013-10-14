   function [theta] = tools_plotsquircles(w,circle,stimulus,nstims,radi,S,C)
        rect = stimulus.rect;
        center = .5*[rect(3) , rect(4)];
   
        if circle
            % positions of the stimuli on the screen
            if stimulus.rndtheta && nstims<stimulus.rndnstim
                sort_theta = [1,1];
                while any(diff(sort_theta) < stimulus.mindtheta)
                    theta = 360*rand(1,nstims);
                    sort_theta = [sort(theta),0];
                    sort_theta(end) = sort_theta(1)+360;
                end
            else
                theta = linspace(0,360,nstims+1);
                theta = theta(1:end-1);
            end
            % [x,y] coordinates
            x        = center(1) + radi*sin(theta*pi/180);
            y        = center(2) + radi*cos(theta*pi/180);
        else
            nstims = stimulus.n_spectrum;
            x = stimulus.x_spectrum;
            y = stimulus.y_spectrum;
            S = stimulus.s_spectrum;
            C = stimulus.c_spectrum;
            theta = [];
        end
        
        % graphical construction
        Harea = stimulus.harea;          % area is kept constant for all shapes
        Ht = 0:(2*pi/180):(2*pi);   % shapes are defined by parametric curves of Ht (angle)
        Hr = 45*pi/180;             % the rotation of the hyperellipse
        for n = 1:nstims
            Ha = sqrt(abs(Harea./4.*gamma(1.+2./S(n))./(gamma(1.+1./S(n))).^2)); % define the parameter a from the expected area + curvature S(n)
            Hx = abs(cos(Ht)).^(2./S(n)).*Ha.*sign(cos(Ht));                     % x coordinates
            Hy = abs(sin(Ht)).^(2./S(n)).*Ha.*sign(sin(Ht));                     % y coordinates
            Hc = [C(n) 0 1-C(n)];                                                % color
            Hpoly = [Hx*cos(Hr)-Hy*sin(Hr) + x(n) ; Hx*sin(Hr)+Hy*cos(Hr) + y(n)]';
            Screen('FillPoly',w,RGBcor(Hc,stimulus.lumiIM),Hpoly,1);
        end
    end