% squircles_paris

    % notes
    %{
        single stimulus feature (shape)
        staircase momentum (mean)
            - convergence for 80% performance
            - intermingled conditions/staircases
        staircases defined by conditions
            - independent momentum (variance)
            - set size
            - constrained/unconstrained sampling
            - time exposition
            - starting from the harder case (d'=0)
        answers only after stimuli exposition
    %}
    
    %% initialization =====================================================
    
    % if a configuration file, run it!
    if exist('config.m','file')
        config();
    end
    
    %% variables
    if ~exist('variables','var')
        variables = struct();
        % conditions
        variables.constraints   = 1;
        variables.setsizes      = [2,4,8];
        variables.timings       = .2;
        variables.variances     = [.1 .2];
        variables.means         = 0:.05:.2;
        variables.counters      = 1:2;
        if min(variables.means) > .3 - max(variables.variances)
            error('squircles_paris: error: some means too big');
        end
        variables.mc            = 0;
        variables.vc            = 0.2;
    end
    max_setsizes = max(variables.setsizes);

    %% experiment
    if ~exist('experiment','var')
        experiment = struct();
        experiment.task = 'shape';                                         % task
        experiment.varying = 'ms';                                         % staircase value is shape mean
        experiment.nb_blocks = 3;                                         % number of blocks
        experiment.audio = 1;                                              % feedback audio
        experiment.text = 0;                                               % feedback text
        experiment.SR = round(rand);                                       % mouse button
        experiment.fixdur = .3;                                            % fixation duration (seconds)
        experiment.maskdur = .3;                                           % mask duration
        experiment.maskvar = 100;                                          % mask color variance
        experiment.pausedur = 0;
        experiment.date_start = datestr(now,'yyyymmddTHHMMSS');            % start time
        experiment.startexp = [];
        experiment.endexp = [];
        experiment.ttmin = 0;
        experiment.tth = 0;
        % response keys
        switch(experiment.SR)
            case 0
                experiment.keyleft = 1;
                experiment.keyright = 3;
            case 1
                experiment.keyleft=3;
                experiment.keyright=1;
        end
    end
    tic;                                                                   % timing
    rand('state',sum(100*clock));                                          % random seed
    
    %% conditions matrix
    nb_means      = length(variables.means);
    nb_blocks     = experiment.nb_blocks;
    nb_blocks_ss  = nb_blocks/length(variables.setsizes);
    nb_conditions = length(variables.constraints) * length(variables.timings) * length(variables.variances) * length(variables.means) * length(variables.counters) * length(variables.setsizes);
    l_blocks      = nb_conditions / nb_blocks;
    if (l_blocks-round(l_blocks)) || mod(nb_blocks,length(variables.setsizes))
        error('squircles_paris: error: nb_blocks incongruent');
    end
    % conditions matrix
    % [i_condition,setsize,variance,constraint,timing,i_mean,i_psych,counter]
    if ~exist('conditions','var')
        % create condition cell (by setsizes)
        i_psych     = 0;
        i_condition = 0;
        conditions_cell = {};
        conditions_matrix = [];
        for setsize = variables.setsizes
            for variance = variables.variances
                for constraint = variables.constraints
                    for timing = variables.timings
                        i_psych = i_psych + 1;
                        for i_mean = 1:length(variables.means)
                            for counter = variables.counters
                                % condition
                                i_condition = i_condition + 1;
                                % condition matrix
                                condition = [i_condition,setsize,variance,constraint,timing,i_mean,i_psych,counter];
                                % trials matrix
                                conditions_matrix = [conditions_matrix ; condition];
                            end
                        end
                    end
                end
            end
            conditions_cell{end+1} = conditions_matrix(randperm(size(conditions_matrix,1)),:);
            conditions_matrix = [];
        end
        % create matrix condition
        conditions = [];
        for i_block = 1:nb_blocks
            ib_block = (mod(i_block-1,nb_blocks_ss)  )*l_blocks + 1;
            ie_block = (mod(i_block-1,nb_blocks_ss)+1)*l_blocks;
            conditions = [conditions ; conditions_cell{ceil(i_block/nb_blocks_ss)}(ib_block:ie_block,:)];
        end
        clear condition
        clear conditions_cell
        clear conditions_matrix
    end
    
    %% participant
    if ~exist('participant','var')
        participant = struct();
        participant.name = input('Nom: ','s');
        participant.age = input('Age: ','s');
        participant.sex = input('Sexe: ','s');
        participant.hand = input('Droitier: ','s');
        participant.subject = Randi(1000);
        participant.filename = {'data_training',sprintf('squirclespsych_%s',participant.name),1};
    end
    
    %% stimulus
    if ~exist('stimulus','var')
        stimulus = struct();
        stimulus.tiny      = 4;                                            % fixation size
        stimulus.radi      = 200;                                          % radius
        stimulus.rndtheta  = 0;                                            % random locations
        stimulus.rndnstim  = 10;                                           % max nstims to use rand locations
        stimulus.mindtheta = 22.5;                                         % minimum angle distance for rand loc stimuli
        stimulus.harea     = 2500;                                         % area of stimulus (is constant)
        stimulus.whitecol  = RGBcor(1,1);                                  % white color
        stimulus.lumiBG    = .2;                                           % background luminance
        stimulus.bgcol     = RGBcor(1,stimulus.lumiBG);                    % background color
        stimulus.lumiIM    = .09;                                          % luminance of each item
        stimulus.crange    = [0.05 0.95];                                  % c=0 blue    -> c=1 red
        stimulus.srange    = [1.05 1.95];                                  % s=1 square  -> s=2 circle
        %stimulus.rect      = [0 0 1280 1024];
    end
    
    %% save
    while exist([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'file')
        participant.filename{3} = participant.filename{3} + 1;
    end
    save([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat']);
    
    %% task ===============================================================
    try
        %% open window
        if ~isfield(stimulus,'rect')
            [w, rect] = Screen('OpenWindow', 0, 0,[],32,2);
            stimulus.rect = rect;
        else
            [w, rect] = Screen('OpenWindow', 0, 0,stimulus.rect,32,2);
        end
        Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        HideCursor;
        Screen('TextFont', w, 'Helvetica');
        Screen('TextSize', w , 16 );
        center = [rect(3)/2 rect(4)/2];                                    % screen center
        tinyrect = [center-stimulus.tiny, center+stimulus.tiny];           % fixation rect
        KbName('UnifyKeyNames')

        %% instructions
        Screen(w,'FillRect',stimulus.bgcol);
        switch(experiment.SR)
            case 0
                Screen(w,'DrawText','SQUARE (left button) or CIRCLE (right button)',      300,350,            stimulus.whitecol);
            case 1
                Screen(w,'DrawText','CIRCLE (left button) or SQUARE (right button)',      300,350,            stimulus.whitecol);
        end
        xcoords = ((1:11) - 5.5)*50 +  rect(3)/2;
        for k = 1:11
            switch(experiment.SR)
                case 0
                    s = 1+(k-1)/10;
                case 1
                    s = 1+(11-k)/10;
            end
            Ht = 0 : (2*pi/180) : (2*pi);                                  % shapes are defined by parametric curves of Ht (angle)
            Hr = 45*pi/180;                                                % the rotation of the hyperellipse
            Ha = sqrt(abs(1500./4.*gamma(1.+2./s)./(gamma(1.+1./s)).^2));  % define the parameter a from the expected area + curvature S(n)
            Hx = abs(cos(Ht)).^(2./s).*Ha.*sign(cos(Ht));                  % x coordinates
            Hy = abs(sin(Ht)).^(2./s).*Ha.*sign(sin(Ht));                  % y coordinates
            Hc = [0.5 0 0.5];                                              % color
            Hpoly = [Hx*cos(Hr)-Hy*sin(Hr) + xcoords(k) ; Hx*sin(Hr)+Hy*cos(Hr) + 425]';
            RGBcolor = RGBcor(Hc,stimulus.lumiIM);
            Screen('FillPoly', w, RGBcolor, Hpoly, 1);
        end
        Screen(w,'DrawText','Hello and thanks for your collaboration!',                             200,100,   stimulus.whitecol);
        Screen(w,'DrawText','In this task you''ll see some coloured squircles.',                    200,150,   stimulus.whitecol);
        Screen(w,'DrawText','You will have to decide as fast as possible the shape',                200,200,   stimulus.whitecol);
        Screen(w,'DrawText','and to say if it looks more like a ...',                               200,250,   stimulus.whitecol);
        Screen(w,'DrawText','You will hear a high bip when the answer is correct',                  200,500,   stimulus.whitecol);
        Screen(w,'DrawText','and a low bip when it is not.',                                        200,550,   stimulus.whitecol);
        Screen(w,'DrawText','You might want to use the feedback to improve your responses.',        200,600,   stimulus.whitecol);
        Screen(w,'DrawText','Good luck!',                                                           200,700,   stimulus.whitecol);
        Screen(w,'Flip');       
        buttons = 0;
        while ~any(buttons)
            [x_mouse,y_mouse,buttons]=GetMouse;
        end
        WaitSecs(0.5);

        %% create mask texture
        stim_r = sqrt(.3183 * stimulus.harea);
        [circle_x,circle_y] = meshgrid(rect(1):rect(3),rect(2):rect(4));
        circle_d = sqrt((circle_x-center(1)).^2 + (circle_y-center(2)).^2); % distance matrix
        r_min = stimulus.radi - stim_r;                                     % annulus min radius
        r_max = stimulus.radi + stim_r;                                     % annulus max radius
        circle_a = 255*(circle_d<r_min | circle_d>r_max);                   % alpha (circular) layer
        
        %% blank wait
        Screen(w,'Flip');       
        buttons = 0;
        while ~any(buttons)
            [x_mouse,y_mouse,buttons]=GetMouse;
        end
        WaitSecs(0.5);
        
        %% create data struct
        if ~exist('data','var')
            data = struct();
            data.sub          = [];
            data.i_trial      = 0;
            data.i_block      = 0;
            data.RT           = [];
            data.cor          = [];
            data.cors         = cell(nb_conditions);
            data.err          = [];
            data.Mc           = [];
            data.Ms           = [];
            data.Vc           = [];
            data.Vs           = [];
            data.tasknum      = [];
            data.setsize      = [];
            data.C            = [];
            data.S            = [];
            data.Ccat         = [];
            data.Scat         = [];
            data.SR           = [];
            data.theta        = [];
            data.fix_terror   = [];
            data.stim_terror  = [];
            data.mask_terror  = [];
            data.pause_terror = [];
            data.stimcat      = [];
            data.respcat      = [];
            data.respcode     = [];
            data.condition    = [];
            data.psych        = [];
        end
        
        %% here we go
        experiment.startexp(end+1) = GetSecs;
        getout = 0;
        i_block = ceil((data.i_trial+1)/l_blocks)-1;
        for i_trial = (data.i_trial+1):size(conditions,1)
            
            % new block
            if ceil(i_trial/l_blocks) > i_block
                % new block
                i_block = i_block + 1;
                data.i_block = i_block;
                % block screen 1
                Screen(w,'DrawText',['Block ',num2str(i_block),' / ',num2str(nb_blocks)],    center(1)-100,  center(2)     );
                Screen(w,'Flip');
                buttons = 0;
                while ~any(buttons)
                    [x_mouse,y_mouse,buttons] = GetMouse;
                end
                WaitSecs(0.5);
                % block screen 1
                Screen(w,'DrawText',['Block ',num2str(i_block),' / ',num2str(nb_blocks)],    center(1)-100,  center(2)     );
                Screen(w,'DrawText', 'Click to continue'  ,    center(1)-100,  center(2)+200 );
                Screen(w,'Flip');
                buttons = 0;
                while ~any(buttons)
                    [x_mouse,y_mouse,buttons] = GetMouse;
                end
                WaitSecs(0.5);
            end
            
            % conditions --------------------------------------------------
            nstims      = conditions(i_trial,2);
            vs          = conditions(i_trial,3);
            constraint  = conditions(i_trial,4);
            stimdur     = conditions(i_trial,5);
            ms          = variables.means(conditions(i_trial,6));
            mc = variables.mc;
            vc = variables.vc;
            
            % stimulus ----------------------------------------------------
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
            x        = sin(theta*pi/180);
            y        = cos(theta*pi/180);
            % decide stimulus category
            Ccat = sign(rand-0.5); % color
            Scat = sign(rand-0.5); % shape
            % actual values for color and shape
            %                      mean          var     nitems  constraint    critmean    critsd    range              nsamples
            C = mypseudorandrange( 0.5+Ccat*mc,  vc,     nstims, constraint,   0.001,      0.001,    stimulus.crange,   1);
            S = mypseudorandrange( 1.5+Scat*ms,  vs,     nstims, constraint,   0.001,      0.001,    stimulus.srange,   1);

            % fixation ----------------------------------------------------
            Screen(w,'FillRect',stimulus.bgcol);
            Screen(w,'FillOval',stimulus.whitecol,tinyrect);
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip'); % show fixation
            when = StimulusOnsetTime + experiment.fixdur;

            % exposition -------------------------------------------------
            % graphical construction
            Harea = stimulus.harea;          % area is kept constant for all shapes
            Ht = 0:(2*pi/180):(2*pi);   % shapes are defined by parametric curves of Ht (angle)
            Hr = 45*pi/180;             % the rotation of the hyperellipse
            for n = 1:nstims
                Ha = sqrt(abs(Harea./4.*gamma(1.+2./S(n))./(gamma(1.+1./S(n))).^2)); % define the parameter a from the expected area + curvature S(n)
                Hx = abs(cos(Ht)).^(2./S(n)).*Ha.*sign(cos(Ht));                     % x coordinates
                Hy = abs(sin(Ht)).^(2./S(n)).*Ha.*sign(sin(Ht));                     % y coordinates
                Hc = [C(n) 0 1-C(n)];                                                % color
                Hpoly = [Hx*cos(Hr)-Hy*sin(Hr) + center(1)+stimulus.radi*x(n) ; Hx*sin(Hr)+Hy*cos(Hr) + center(2)+stimulus.radi*y(n)]';
                Screen('FillPoly',w,RGBcor(Hc,stimulus.lumiIM),Hpoly,1);
            end
            % fixation
            Screen(w,'FillOval',stimulus.whitecol,tinyrect);
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when); % show stimulus
            fix_terror = StimulusOnsetTime - when;  % fixation time error
            when = StimulusOnsetTime + stimdur;
            
            % masking -----------------------------------------------------
            % noise
            rectsize = [rect(3)-rect(1),rect(4)-rect(2)];
            noiseimg1 = (experiment.maskvar*randn(rectsize(1),rectsize(2)) + 128);
            noiseimg2 = zeros(rectsize(1),rectsize(2));
            noiseimg3 = (experiment.maskvar*randn(rectsize(1),rectsize(2)) + 128);
            noiseimg(:,:,1) = noiseimg1;
            noiseimg(:,:,2) = noiseimg2;
            noiseimg(:,:,3) = noiseimg3;
            mtex = Screen('MakeTexture', w, noiseimg);
            Screen('DrawTexture', w, mtex, [], rect, [], 0);
            Screen('Close', mtex);
            % alpha circle
            circleimg = stimulus.bgcol(1)*ones(size(circle_d));                 % color (background) layer. has to be greyscale!!
            circleimg(:,:,2) = circle_a;                                        % alpha (circular) layer
            ctex = Screen('MakeTexture', w, circleimg);
            Screen('DrawTexture', w, ctex, [], rect, [], 0);
            Screen('Close', ctex);
            % fixation
            Screen(w,'FillOval',stimulus.whitecol,tinyrect);
            % flip
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when); % show masking
            stim_terror = StimulusOnsetTime - when;
            when = StimulusOnsetTime + experiment.maskdur;

            % pause -------------------------------------------------------
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when); % show masking
            mask_terror = StimulusOnsetTime - when;
            when = GetSecs + experiment.pausedur;
            
            % response ---------------------------------------------------
            Screen(w,'DrawText',[num2str(mod(i_trial-1,l_blocks)+1),'/',num2str(l_blocks)],center(1),rect(4)-200);
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when); % show masking
            pause_terror = StimulusOnsetTime - when;
            respcode = [];
            RT = [];
            startresp = GetSecs;
            while isempty(respcode)
                [kdown,ksecs,codes] = KbCheck();
                [x_mouse,y_mouse,buttons] = GetMouse();
                % check escape key
                if kdown && (any(strcmpi(KbName(codes),'escape')) || any(strcmpi(KbName(codes),'esc')))
                    getout = 1;
                    break
                end
                % check mouse
                if any(buttons([1,3]))
                    respcode = find(buttons,1);
                    RT = GetSecs-startresp;
                end           
            end
            if getout
                break
            end

            % accuracy and feedback ---------------------------------------
            % accuracy
            if isempty(respcode);
                respcode = 0;
                cor = 0;
                RT = 0;
            else
                switch(experiment.task)
                    case 'color'
                        % correct
                        cor = (Ccat==1 && respcode==experiment.keyright) || (Ccat==-1 && respcode==experiment.keyleft);   % red=left,  blue=right
                    case 'shape'
                        % correct
                        cor = (Scat==1 && respcode==experiment.keyright) || (Scat==-1 && respcode==experiment.keyleft);   % vert=left, horz=right
                end
            end
            % auditory feedback
            if experiment.audio
                if cor==1
                    note2(1200,0.1);
                else
                    note2(400,0.1);
                end
            end
            % visual feedback
            if experiment.text
                if cor==1
                    Screen(w,'DrawText','YES',center(1),center(2));
                else
                    Screen(w,'DrawText','NO',center(1),center(2));
                end
            end
            Screen(w,'Flip');
            
            % report ------------------------------------------------------
            fprintf('trial %d -- key %d -- RT %.2f -- cor %d \n',i_trial,respcode,RT, cor);
            fprintf('        [ i=%d, s=%d, v=%.2f, c=%d, t=%.2f, m=%.2f] \n', conditions(i_trial,1:5),ms);
            fprintf('\n');
            
            % log data ----------------------------------------------------
            data.sub(i_trial)          = participant.subject;              % subject
            data.i_trial               = i_trial;                          % trial number
            data.RT(i_trial)           = RT;                               % reaction times
            data.cor(i_trial)          = cor;                              % correct
            data.cors{conditions(i_trial,6)}(end+1) = cor;                 % corrects (for each condition)
            data.err(i_trial)          = 1-cor;                            % error
            data.Mc(i_trial)           = mc;                               % mean/color
            data.Ms(i_trial)           = ms;                               % mean/shape
            data.Vc(i_trial)           = vc;                               % var /color
            data.Vs(i_trial)           = vs;                               % var /shape
            data.tasknum(i_trial)      = strcmp(experiment.task,'shape');  % 0='color', 1='shape'
            data.setsize(i_trial)      = nstims;                           % nstims (constant)
            data.C(1:nstims,i_trial)   = C;                                % color value
            data.S(1:nstims,i_trial)   = S;                                % shape value
            data.Ccat(i_trial)         = Ccat;                             % -1=categ1 (blue),   1=categ2 (red)
            data.Scat(i_trial)         = Scat;                             % -1=categ1 (square), 1=categ2 (circle)
            data.SR(i_trial)           = experiment.SR;                    % mouse buttons coding
            data.theta(:,i_trial)  = [theta,zeros(1,max_setsizes-nstims)]; % angle for stimuli
            data.fix_terror(i_trial)   = fix_terror;                       % time flip error in fixation screen
            data.stim_terror(i_trial)  = stim_terror;                      % time flip error in stimulus screen
            data.mask_terror(i_trial)  = mask_terror;                      % time flip error in mask screen
            data.pause_terror(i_trial) = pause_terror;                     % time flip error in pause screen
            data.stimcat(i_trial)      = data.Scat(i_trial).*(data.tasknum(i_trial)==1) + data.Ccat(i_trial).*(data.tasknum(i_trial)==0); % category (for the relevant feature)
            data.respcat(i_trial)      = (respcode-2).*(respcode~=0).*sign(experiment.SR+0.5); % resp = -1 (categ1) or 1 (categ2) or 0 (none)
            data.respcode(i_trial)     = respcode;
            data.condition(i_trial)    = conditions(i_trial,1);            % condition index
            data.psych(i_trial)        = conditions(i_trial,7);            % psychophysic function index
            
            %% save
            save([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'conditions','data','experiment','participant','stimulus','variables');
            
        end
    catch e
        save('tmp_error.mat');
        ShowCursor;
        Screen('CloseAll');
        FlushEvents;
        rethrow(e);
    end
    
    %% eoe
    experiment.endexp(end+1) = GetSecs;
    if length(experiment.startexp)==length(experiment.endexp)
        experiment.ttmin = sum(experiment.endexp-experiment.startexp)/60;
        experiment.tth = sum(experiment.endexp-experiment.startexp)/3600;
    end
    
    %% save
    save([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'conditions','data','experiment','participant','stimulus','variables');
    clearvars -except getout conditions data experiment participant stimulus variables;
    
    %% finalize
    ShowCursor;
    Screen('CloseAll');
    FlushEvents;
    
    if getout
        return;
    end

    %% postprocessing psychophysic functions ==============================
    nb_psych = max(data.psych);
    nb_means = length(variables.means);
    
    data.mcor = zeros(nb_psych,nb_means);
    for i_psych = 1:nb_psych
        for i_mean = 1:nb_means
            data.mcor(i_psych,i_mean) = mean(data.cor(data.psych==i_psych & data.Ms==variables.means(i_mean)));
        end
    end
    

    %% plotting ===========================================================
    close all;
    colors = 'bgrcmyk';
    nb_setsizes = length(variables.setsizes);
    nb_timings  = length(variables.timings);
    nb_variances = length(variables.variances);
    nb_constraints = length(variables.constraints);
    
    %% plot psychophysic functions ----------------------------------------
    figure;
    plot(variables.means,data.mcor);
    
    %% fitting curves -----------------------------------------------------
    fits_p = [];
    fits_y = [];
    minx = min(variables.means);
    maxx = max(variables.means);
    x = minx:(maxx-minx)/20:maxx;
    for i_mcor = 1:length(data.mcor)
        shifted_weibull =  @(p,x) 1-p(4) + p(4).*(p(2)/p(1)).*power((x-p(3))/p(1),p(2)-1).*exp(-power((x-p(3))/p(1),p(2)));
        p = nlinfit(variables.means,data.mcor(i_mcor,:),shifted_weibull,[1,1,0,0]);
        fits_p = [fits_p ; p];
        fits_y = [fits_y ; shifted_weibull(p,x)];
    end
    figure;
    plot(x,fits_y);
    
