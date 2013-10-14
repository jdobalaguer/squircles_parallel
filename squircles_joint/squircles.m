% general squircles

%{
    notes :
        experiment.task = S [shape]
        experiment.task = C [color]
        experiment.task = B [both]
%}
    
    %% initialization =====================================================
    new_subject = 0;
    
    % if a configuration file, run it!
    if exist('config.m','file')
        config();
    end
    
    %% experiment
    if ~exist('experiment','var')
        new_subject = 1;
        experiment = struct();
        experiment.task = input('task: ','s');                             % task
        if length(experiment.task)~=1 || ~any(experiment.task=='SCB')
            clear all;
            error('squircles: error: task doesn''t exist');
        end
        experiment.training = str2double(input('training: ','s'));         % training
        if experiment.training
            experiment.nb_blocks = 12;                                     % number of blocks
        else
            experiment.nb_blocks = 80;                                     % number of blocks
        end
        experiment.audio = 0;                                              % feedback audio
        experiment.text = 1;                                               % feedback text
        if experiment.training
            experiment.SR = round(rand);                                   % mouse button
        else
            experiment.SR = str2num(input('SR: ','s'));                    % mouse button
        end
        experiment.spectrum = str2double(input('spectrum: ','s'));         % spectrum
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
    
    if experiment.audio
        audio_port = tools_audioport();
        audio_port.init();
    end
    
    %% variables
    if ~exist('variables','var')
        new_subject = 1;
        variables = struct();
        if experiment.training
            % training conditions
            variables.constraints   = 1;
            variables.setsizes      = [1,3,6,12];
            variables.radis         = str2double(input('radi(130,200): ','s'));
            variables.timings       = .05;%[.4,.2,.1,.05];
            variables.variances     = .15;%.1;
            variables.means         = -.2:.05:.2;
            variables.counters      = 1;
            if min(variables.means) > .3 - max(variables.variances)
                error('squircles_paris: error: some means too big');
            end
        else
            % experiment conditions
            variables.constraints   = 1;
            variables.setsizes      = [1,3,6,12];
            variables.radis         = str2double(input('radi(130,200): ','s'));
            variables.timings       = [.05,.1,.2,.4];
            variables.variances     = .1;
            variables.means         = -.2:.05:.2;
            variables.counters      = 1:15;
            if min(variables.means) > .3 - max(variables.variances)
                error('squircles_paris: error: some means too big');
            end
        end
    end
    max_setsizes = max(variables.setsizes);

    %% check nb_blocks
    nb_means      = length(variables.means);
    nb_blocks     = experiment.nb_blocks;
    nbmin_blocks  = length(variables.setsizes)*length(variables.radis)*length(variables.timings);
    nb_blocks_ss  = nb_blocks/nbmin_blocks;
    nb_conditions = length(variables.constraints) * length(variables.timings) * length(variables.variances) * length(variables.means) * length(variables.counters) * length(variables.setsizes) * length(variables.radis);
    l_blocks      = nb_conditions / nb_blocks;
    if (l_blocks-round(l_blocks))
        fprintf(['squircles_paris: error: nb_conditions = ',num2str(nb_conditions),' \n',...
                 '                        nb_blocks     = ',num2str(nb_blocks),'\n',...
                 '                        has to be multiple of ',num2str(nbmin_blocks),'\n',...
                 '                        l_blocks      = ',num2str(l_blocks),'\n',...
                 '                        nb_blocks_ss  = ',num2str(nb_blocks_ss),'\n']);
        variables
        clear all;
        error('squircles_paris: error: nb_blocks incongruent');
    end
    if (nb_blocks_ss-round(nb_blocks_ss))
        fprintf(['squircles_paris: error: nb_conditions = ',num2str(nb_conditions),' \n',...
                 '                        nb_blocks     = ',num2str(nb_blocks),'\n',...
                 '                        has to be multiple of ',num2str(nbmin_blocks),'\n',...
                 '                        l_blocks      = ',num2str(l_blocks),'\n',...
                 '                        nb_blocks_ss  = ',num2str(nb_blocks_ss),'\n']);
        variables
        clear all;
        error('squircles_paris: error: nb_blocks incongruent');
    end
    
    %% conditions matrix
    % conditions matrix
    % [i_condition,setsize,i_radi,variance,constraint,timing,i_mean,i_psych,counter]
    if ~exist('conditions','var')
        new_subject = 1;
        % create condition cell (by setsizes)
        i_psych     = 0;
        i_condition = 0;
        conditions_cell = {};
        conditions_matrix = [];
        for setsize = variables.setsizes
            for i_radi = 1:length(variables.radis)
                for timing = variables.timings
                    for variance = variables.variances
                        for constraint = variables.constraints
                            i_psych = i_psych + 1;
                            for i_mean = 1:length(variables.means)
                                for counter = variables.counters
                                    % condition
                                    i_condition = i_condition + 1;
                                    % condition matrix
                                    condition = [i_condition,setsize,i_radi,variance,constraint,timing,i_mean,i_psych,counter];
                                    % trials matrix
                                    conditions_matrix = [conditions_matrix ; condition];
                                end
                            end
                        end
                    end
                    conditions_cell{end+1} = conditions_matrix(randperm(size(conditions_matrix,1)),:);
                    conditions_matrix = [];
                end
            end
        end
        % create matrix condition
        conditions = [];
        if experiment.training
            v_blocks = 1:nb_blocks;
        else
            v_blocks = randperm(nb_blocks);
        end
        for i_block = v_blocks
            ib_block = (mod(i_block-1,nb_blocks_ss)  )*l_blocks + 1;
            ie_block = (mod(i_block-1,nb_blocks_ss)+1)*l_blocks;
            conditions = [conditions ; conditions_cell{ceil(i_block/nb_blocks_ss)}(ib_block:ie_block,:)];
        end
        clear v_blocks
        clear condition
        clear conditions_cell
        clear conditions_matrix
    end
    
    %% participant
    if ~exist('participant','var')
        new_subject = 1;
        participant = struct();
        participant.name = input('Nom: ','s');
        participant.age = input('Age: ','s');
        participant.sex = input('Sexe: ','s');
        participant.hand = input('Droitier: ','s');
        participant.subject = Randi(1000);
        if experiment.training
            participant.filename = {['datatra_',experiment.task],sprintf('squirclespsych_%s',participant.name),1};
        else
            participant.filename = {['dataexp_',experiment.task],sprintf('squirclespsych_%s',participant.name),1};
        end
    end
    
    %% stimulus
    if ~exist('stimulus','var')
        new_subject = 1;
        stimulus = struct();
        % question stimulus
        stimulus.tiny      = 4;                                            % fixation size
        stimulus.rndtheta  = 0;                                            % random locations
        stimulus.rndnstim  = 10;                                           % max nstims to use rand locations
        stimulus.mindtheta = 22.5;                                         % minimum angle distance for rand loc stimuli
        stimulus.harea     = 2500;                                         % area of stimulus (is constant)
        stimulus.whitecol  = tools_RGBcor(1,1);                            % white color
        stimulus.lumiBG    = .2;                                           % background luminance
        stimulus.bgcol     = tools_RGBcor(1,stimulus.lumiBG);                    % background color
        stimulus.lumiIM    = .09;                                          % luminance of each item
        stimulus.vrange    = [-.95 +.95];                                  % value range
        stimulus.crange    = [0.05 0.95];                                  % c=0 blue    -> c=1 red
        stimulus.srange    = [1.05 1.95];                                  % s=1 square  -> s=2 circle
        % response stimulus
        stimulus.n_spectrum = length(variables.means);
        stimulus.x_spectrum = [];
        stimulus.y_spectrum = [];
        stimulus.v_spectrum = [];
        stimulus.s_spectrum = [];
        stimulus.c_spectrum = [];
        % screen
        %stimulus.rect      = [0 0 1280 1024]; %lab screen
        %stimulus.rect      = [0 0 1024 768];  %c125 screen
    end
    
    %% save
    while exist([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'file')
        if new_subject
            clear all;
            error('squircles_paris: error: name already in use');
        end
        participant.filename{3} = participant.filename{3} + 1;
    end
    save([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat']);
    
    %% task ===============================================================
    try
        %% open window
        %Screen('Preference', 'SkipSyncTests', 2);
        if ~isfield(stimulus,'rect') || isempty(stimulus.rect)
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
        switch(experiment.task)
            case 'S'
                switch(experiment.SR)
                    case 0
                        Screen(w,'DrawText','SQUARE or CIRCLE',      300,300,            stimulus.whitecol);
                        stimulus.v_spectrum = variables.means;
                        stimulus.s_spectrum = .5*stimulus.v_spectrum + 1.5;
                        stimulus.c_spectrum = 0.5*ones(1,stimulus.n_spectrum);
                    case 1
                        Screen(w,'DrawText','CIRCLE or SQUARE',      300,300,            stimulus.whitecol);
                         stimulus.v_spectrum = fliplr(variables.means);
                        stimulus.s_spectrum = .5*stimulus.v_spectrum + 1.5;
                        stimulus.c_spectrum = 0.5*ones(1,stimulus.n_spectrum);
               end
            case 'C'
                switch(experiment.SR)
                    case 0
                        Screen(w,'DrawText','BLUE or RED',      300,300,            stimulus.whitecol);
                        stimulus.v_spectrum = variables.means;
                        stimulus.s_spectrum = 1.5*ones(1,stimulus.n_spectrum);
                        stimulus.c_spectrum = .5*stimulus.v_spectrum + 0.5;
                    case 1
                        Screen(w,'DrawText','RED or BLUE',      300,300,            stimulus.whitecol);
                        stimulus.v_spectrum = fliplr(variables.means);
                        stimulus.s_spectrum = 1.5*ones(1,stimulus.n_spectrum);
                        stimulus.c_spectrum = .5*stimulus.v_spectrum + 0.5;
                end
            case 'B'
                % TODO ##############################################
                error('squircles: error: task B still to be done');
        end
        % store [x,y] coordinates in stimulus
        stimulus.x_spectrum = 60*linspace(-.5*stimulus.n_spectrum,+.5*stimulus.n_spectrum,stimulus.n_spectrum) +  rect(3)/2;
        stimulus.y_spectrum = ones(1,stimulus.n_spectrum) * center(2);
        % plot squircles
        tools_plotsquircles(w,0,stimulus,[],[],[],[]);
        % plot text
        Screen(w,'DrawText','Hello and thanks for your collaboration!',                             200,100,   stimulus.whitecol);
        Screen(w,'DrawText','In this task you''ll see some coloured squircles.',                    200,150,   stimulus.whitecol);
        Screen(w,'DrawText','You will have to decide as fast as possible the shape',                200,200,   stimulus.whitecol);
        Screen(w,'DrawText','and to say if it looks more like a ...',                               200,250,   stimulus.whitecol);
        Screen(w,'DrawText','You will hear a high bip when the answer is correct',                  200,500,   stimulus.whitecol);
        Screen(w,'DrawText','and a low bip when it is not.',                                        200,550,   stimulus.whitecol);
        Screen(w,'DrawText','You might want to use the feedback to improve your responses.',        200,600,   stimulus.whitecol);
        Screen(w,'DrawText','Good luck!',                                                           200,700,   stimulus.whitecol);
        Screen(w,'Flip');
imwrite(Screen('GetImage', w, rect), ['image_',num2str(randi(10000)),'.png']);
        buttons = 0;
        while ~any(buttons)
            [x_mouse,y_mouse,buttons]=GetMouse;
        end
        WaitSecs(0.5);

        %% create mask texture
        circle_a = cell(1,length(variables.radis));
        stim_r = sqrt(.3183 * stimulus.harea);
        [circle_x,circle_y] = meshgrid(rect(1):rect(3),rect(2):rect(4));
        circle_d = sqrt((circle_x-center(1)).^2 + (circle_y-center(2)).^2);% distance matrix
        for i_radi = 1:length(variables.radis)
            r_min = variables.radis(i_radi) - stim_r;                       % annulus min radius
            r_max = variables.radis(i_radi) + stim_r;                       % annulus max radius
            circle_a{i_radi} = 255*(circle_d<r_min | circle_d>r_max);      % alpha (circular) layer
        end
        
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
            
            data.i_trial      =  0;
            data.i_block      =  0;

            % experiment
            data.exp_sub                = [];
            data.exp_SR                 = [];
            data.exp_trial              = [];
            data.exp_block              = [];
            % screen
            data.screen_theta           = [];
            data.screen_fix_terror      = [];
            data.screen_stim_terror     = [];
            data.screen_mask_terror     = [];
            data.screen_pause_terror    = [];
            % response
            data.resp_RT                = [];
            data.resp_code              = [];
            data.resp_cat               = [];
            data.resp_cor               = [];
            data.resp_err               = [];
            % variables
            data.vb_condition           = [];
            data.vb_setsize             = [];
            data.vb_radi                = [];
            data.vb_var                 = [];
            data.vb_constraint          = [];
            data.vb_timing              = [];
            data.vb_mean                = [];
            data.vb_psych               = [];
            data.vb_counter             = [];
            % sampling
            data.sample_V               = [];
            data.sample_C               = [];
            data.sample_S               = [];
            data.sample_mv              = [];
            data.sample_mc              = [];
            data.sample_ms              = [];
            data.sample_Vcat            = [];
            data.sample_Ccat            = [];
            data.sample_Scat            = [];
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
            condition   = conditions(i_trial,1);
            nstims      = conditions(i_trial,2);
            i_radi      = conditions(i_trial,3);
            radi        = variables.radis(i_radi);
            var_v       = conditions(i_trial,4);
            constraint  = conditions(i_trial,5);
            stimdur     = conditions(i_trial,6);
            mean_v      = variables.means(conditions(i_trial,7));
            psych       = conditions(i_trial,8);
            counter     = conditions(i_trial,9);
            
            % stimulus ----------------------------------------------------
            % actual values for color and shape
            % tools_mypseudorandrange  (mean,  var,  nitems,constraint,crit_mean,crit_sd,range,          nsamples)
            V = tools_mypseudorandrange(mean_v,var_v,nstims,constraint,0.001,    0.001,  stimulus.vrange,1);
            switch experiment.task
                case 'S'
                    % 1-Dimensional SHAPE task
                    S = .5*V + 1.5;                                                % ignore stimulus.srange
                    C = stimulus.crange(1) + rand(1,nstims)*stimulus.crange(2); % uniform distribution
                case 'C'
                    % 1-Dimensional COLOR task
                    C = .5*V + 0.5;                                                % ignore stimulus.crange
                    S = stimulus.srange(1) + rand(1,nstims)*stimulus.srange(2); % uniform distribution
                case 'B'
                    % 2-Dimensional task
                    minC = max(abs(V),stimulus.crange(1));
                    C = [];
                    S = [];
                    i_CS = 1:nstims;
                    while ~isempty(i_CS)
                        % uniform distributions for C and S
                        C(i_CS) = Ccats(i_CS) .* ((stimulus.crange(2)-minC(i_CS)) .* rand(1,length(i_CS)) + minC(i_CS));
                        S(i_CS) = V(i_CS)./C(i_CS);
                        C(i_CS) = .5*C(i_CS) + 0.5;
                        S(i_CS) = .5*S(i_CS) + 1.5;
                        i_CS = find(S>stimulus.srange(2) | S<stimulus.srange(1) | C>stimulus.crange(2) | C<stimulus.crange(1));
                    end
                    clear minC;
                    clear i_CS;
            end
            mv = mean(V);
            ms = mean(S);
            mc = mean(C);
            Vcat = sign(mv);
            Ccat = sign(mc-0.5);
            Scat = sign(ms-1.5);
            % not really good, but... fixing: when ss=1 and m=0 Vcat is 0!
            if Vcat==0 Vcat = sign(randi(2)-1.5); end
            if Scat==0 Scat = sign(randi(2)-1.5); end
            if Ccat==0 Ccat = sign(randi(2)-1.5); end
            
            % fixation ----------------------------------------------------
            Screen(w,'FillRect',stimulus.bgcol);
            Screen(w,'FillOval',stimulus.whitecol,tinyrect);
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip'); % show fixation
            when = StimulusOnsetTime + experiment.fixdur;
    
            % exposition -------------------------------------------------
            % plot the set of squircles
            [theta] = tools_plotsquircles(w,1,stimulus,nstims,radi,S,C);

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
            circleimg = stimulus.bgcol(1)*ones(size(circle_d));            % color (background) layer. has to be greyscale!!
            circleimg(:,:,2) = circle_a{i_radi};                           % alpha (circular) layer
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
            [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when); % show pause
            mask_terror = StimulusOnsetTime - when;
            when = GetSecs + experiment.pausedur;
            
            % response ---------------------------------------------------
            respcode = [];
            respcat = [];
            RT = [];
            % quantitative response
            if experiment.spectrum
                % pause
                [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when);
                pause_terror = StimulusOnsetTime - when;
                % randomize mouse position...
                i_spectrum = randi(stimulus.n_spectrum);
                SetMouse(round(stimulus.x_spectrum(i_spectrum)),round(stimulus.y_spectrum(i_spectrum)));
                % getting the answer...
                xrect = rect([1,3]);
                mxrect = mean(xrect);
                dxrect = .5*diff(xrect);
                yrect = rect([1,3]);
                myrect = mean(yrect);
                dyrect = .5*diff(yrect);
                startresp = GetSecs;
                while isempty(respcode)
                    [kdown,ksecs,codes] = KbCheck();
                    % map mouse values for standard value [-1,+1] response
                    [x_mouse,y_mouse,buttons] = GetMouse();
                    % plot text
                    Screen(w,'DrawText',[num2str(mod(i_trial-1,l_blocks)+1),'/',num2str(l_blocks)],center(1),rect(4)-200);
                    % detect the answer
                    [~,i_spectrum] = min(abs(stimulus.x_spectrum - x_mouse));
                    % plot the answer
                    xans = stimulus.x_spectrum(i_spectrum);
                    yans = stimulus.y_spectrum(i_spectrum);
                    Screen('FillRect',w,100,[xans-30,yans-30,xans+30,yans+30]);
                    % plot squircle
                    tools_plotsquircles(w,0,stimulus,[],[],[],[]);
                    % plot flip
                    [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip');
                    % check escape key
                    if kdown && (any(strcmpi(KbName(codes),'escape')) || any(strcmpi(KbName(codes),'esc')))
                        getout = 1;
                        break
                    end
                    % check mouse
                    if any(buttons([1,3]))
                        respcode = find(buttons,1);
                        respcat = stimulus.v_spectrum(i_spectrum);
                        RT = GetSecs-startresp;
                    end           
                end
                if getout
                    break
                end
            % categorial response
            else
                Screen(w,'DrawText',[num2str(mod(i_trial-1,l_blocks)+1),'/',num2str(l_blocks)],center(1),rect(4)-200);
                [VBLTimestamp,StimulusOnsetTime] = Screen(w,'Flip',when); % show response
                pause_terror = StimulusOnsetTime - when;
                
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
                        respcat = -(respcode-2).*(respcode~=0).*sign(experiment.SR-0.5);
                        RT = GetSecs-startresp;
                    end           
                end
                if getout
                    break
                end
            end

            % check response
            if isempty(respcode);
                respcode = 0;
                respcat = 0;
                RT = 0;
            end
            
            % accuracy and feedback ---------------------------------------
            % quantitative feedback
            if experiment.spectrum
                % accuracy
                err = respcat-mean_v;
                cor = (err==0);
                % auditory feedback
                if experiment.audio
                    audio_port.beep(cor);
                end
                % plot participant's answer
                Screen('FillRect',w,100,[xans-30,yans-30,xans+30,yans+30]);
                % visual feedback
                [~,i_spectrum] = min(abs(stimulus.v_spectrum - mean_v));
                % plot the actual answer
                xans = stimulus.x_spectrum(i_spectrum);
                yans = stimulus.y_spectrum(i_spectrum);
                Screen('FillRect',w,200,[xans-30,yans-30,xans+30,yans+30]);
                % plot squircle
                tools_plotsquircles(w,0,stimulus,[],[],[],[]);
                % plot flip
                Screen(w,'Flip');
imwrite(Screen('GetImage', w, rect), ['image_',num2str(randi(10000)),'.png']);
            % categorial response
            else
                % accuracy
                experiment.SR
                cor = (Vcat==respcat);
                err = (cor==0);
                % auditory feedback
                if experiment.audio
                    audio_port.beep(cor);
                end
                % visual feedback
                if experiment.text
                    % quantitative feedback
                    if cor==1
                        Screen(w,'DrawText','YES',center(1),center(2));
                    else
                        Screen(w,'DrawText','NO',center(1),center(2));
                    end
                    Screen(w,'Flip');
                end
            end
            % wait
            if experiment.text
                while any(buttons)
                    [~,~,buttons] = GetMouse();
                end
                while ~any(buttons)
                    [~,~,buttons] = GetMouse();
                end
            end

            % report ------------------------------------------------------
            fprintf('trial %d -- key %d -- RT %.2f -- err %.2f \n',i_trial,respcode,RT, err);
            fprintf('        [ i=%d, ss=%d, i_r=%d, v=%.2f, c=%d, t=%.2f, m=%.2f] \n', conditions(i_trial,1:6),mv);
            fprintf('\n');
            
            % log data ----------------------------------------------------
            data.i_trial = i_trial;                                                 % trial number
            
            % experiment
            data.exp_sub(i_trial)               = participant.subject;              % subject
            data.exp_SR(i_trial)                = experiment.SR;                    % mouse buttons coding
            data.exp_trial(i_trial)             = i_trial;                          % trial number
            data.exp_block(i_trial)             = i_block;                          % block number
            % stimulus
            data.screen_theta(:,i_trial)        = [theta,zeros(1,max_setsizes-nstims)]; % angle for stimuli
            data.screen_fix_terror(i_trial)     = fix_terror;                       % time flip error in fixation screen
            data.screen_stim_terror(i_trial)    = stim_terror;                      % time flip error in stimulus screen
            data.screen_mask_terror(i_trial)    = mask_terror;                      % time flip error in mask screen
            data.screen_pause_terror(i_trial)   = pause_terror;                     % time flip error in pause screen
            % response
            data.resp_RT(i_trial)               = RT;                               % reaction times
            data.resp_code(i_trial)             = respcode;
            data.resp_cat(i_trial)              = respcat;                          % resp = -1 (categ1) or 1 (categ2) or 0 (none)
            data.resp_cor(i_trial)              = cor;                              % correct
            data.resp_err(i_trial)              = err;                              % error
            % variables
            data.vb_condition(i_trial)          = condition;                        % condition index
            data.vb_setsize(i_trial)            = nstims;                           % nstims (constant)
            data.vb_radi(i_trial)               = radi;                             % circle radius
            data.vb_var(i_trial)                = var_v;                            % value var  (a priori)
            data.vb_constraint(i_trial)         = constraint;                       % constraint
            data.vb_timing(i_trial)             = stimdur;                          % stimdur
            data.vb_mean(i_trial)               = mean_v;                           % value mean (a priori)
            data.vb_psych(i_trial)              = psych;                            % psychophysic function index
            data.vb_counter(i_trial)            = counter;                          % condition counter
            % sampling
            data.sample_V(1:nstims,i_trial)     = V;                                % value value
            data.sample_C(1:nstims,i_trial)     = C;                                % color value
            data.sample_S(1:nstims,i_trial)     = S;                                % shape value
            data.sample_V(nstims+1:end,i_trial) = nan;
            data.sample_C(nstims+1:end,i_trial) = nan;
            data.sample_S(nstims+1:end,i_trial) = nan;
            data.sample_mv(i_trial)             = mv;                               % mean value value
            data.sample_mc(i_trial)             = mc;                               % mean color value
            data.sample_ms(i_trial)             = ms;                               % mean shape value
            data.sample_Vcat(i_trial)           = Vcat;                             % -1=categ1 (blue),   1=categ2 (red)
            data.sample_Ccat(i_trial)           = Ccat;                             % -1=categ1 (blue),   1=categ2 (red)
            data.sample_Scat(i_trial)           = Scat;                             % -1=categ1 (square), 1=categ2 (circle)
            
            %% save
            save([participant.filename{1},filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'conditions','data','experiment','participant','stimulus','variables');
            
        end
    catch e
        save('tmp_error.mat');
        ShowCursor;
        Screen('CloseAll');
        FlushEvents;
        % send a mail
        tools_alertmail();
        % stop sound
        if experiment.audio
            audio_port.stop();
        end
        % rethrow error
        rethrow(e);
    end
    
    %send a mail
    tools_alertmail(participant.name);
    % stop sound
    if experiment.audio
        audio_port.stop();
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
