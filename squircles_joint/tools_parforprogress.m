function percent = tools_parforprogress(N)
%tools_parforprogress Progress monitor (progress bar) that works with parfor.
%   tools_parforprogress works by creating a file called tools_parforprogress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   tools_parforprogress(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   tools_parforprogress updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   tools_parforprogress(0) deletes tools_parforprogress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = tools_parforprogress which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      tools_parforprogress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         tools_parforprogress;
%      end
%      tools_parforprogress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

error(nargchk(0, 1, nargin, 'struct'));

if nargin < 1
    N = -1;
end

percent = 0;
w = 50; % Width of progress bar

if N > 0
    f = fopen('tools_parforprogress.txt', 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fclose(f);
    
    if nargout == 0
        disp(['  0%[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    delete('tools_parforprogress.txt');
    percent = 100;
    
    if nargout == 0
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
else
    if ~exist('tools_parforprogress.txt', 'file')
        error('tools_parforprogress.txt not found. Run tools_parforprogress(N) before tools_parforprogress to initialize tools_parforprogress.txt.');
    end
    
    f = fopen('tools_parforprogress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('tools_parforprogress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
end
