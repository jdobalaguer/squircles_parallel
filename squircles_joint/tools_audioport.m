classdef tools_audioport < handle
    properties
        port
        freq_low
        freq_high
        sampling
        duration
    end 
    methods
        % constructor
        function obj = tools_audioport()
            obj.port = [];
            obj.freq_low  = 2500;
            obj.freq_high = 7500;
            obj.sampling  = [];
            obj.duration  = 0.1;
        end
        % init
        function init(obj)
            PsychPortAudio('Verbosity',0);
            InitializePsychSound;
            try
                obj.port = PsychPortAudio('Open', [], [], 1, [], 1);
            catch e
                psychlasterror('reset');
                obj.port = PsychPortAudio('Open', [], [], 1, [], 1);
            end
            status = PsychPortAudio('GetStatus',obj.port);
            obj.sampling = status.SampleRate;
        end
        % beep
        function beep(obj,high)
            % if not initialised, return
            if isempty(obj.port)
                return
            end
            
            % freq
            if high freq = obj.freq_high;
            else    freq = obj.freq_low;
            end
            % sampling
            sampling = obj.sampling;
            % duration
            duration = obj.duration;
            % create pitch array
            i_ymax = round(duration*sampling);
            y = sin(linspace(0,duration*freq,i_ymax));
            % start the audioport
            PsychPortAudio('FillBuffer',obj.port,y);
            PsychPortAudio('Start',obj.port);
        end
        % close
        function stop(obj)
            PsychPortAudio('Close', obj.port);
            obj.port = [];
        end
    end
end