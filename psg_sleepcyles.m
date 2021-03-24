% PSG_SLEEPCYCLES calculates the sleep cycles and further classifies each
% sleep cycle in descending sleep, deep sleep, ascending sleep and REM
% sleep episodes.
%
% Usage:
%     >> C = psg_sleepcycles(stages, definition, varargin);
%
% Inputs:
%     'stages'     - [vector] or {vector} contains the sleep stage for each
%                    30-second epoch
%     'definition' - [vector] or {vector} contains the integers or 
%                    charactersused to define each sleep epoch in the order 
%                    'Wake', 'NREM 1', 'NREM 2', 'NREM 3', 'REM'.
% Optional:
%     'RecStart' - [datenum] of the start of the recording
% Output:
%     'C.sleepcycle'   - [vector] contains an integer value for each sleep 
%                        epoch that denotes the sleep cycle
%     'C.sleepepisode' - [vector] contains an integer value for each sleep 
%                        epoch that denotes the sleep period.
%                        -1 = Descending sleep episode
%                        -2 = Deep sleep episode
%                         1 = Ascending (transition to REM) sleep episode
%                         2 = REM sleep episode
%     'C.sleepstage'   - [vector] contains the translated integer value for 
%                        each sleep stage
% Examples:
%     >> stages = {'W', 'W', ..., 'N1', 'N2', 'N3', ..., 'R' ...};
%     >> definition = {'W', 'N1', 'N2', 'N3', 'R'};
%     >> C = psg_sleepcycles(stages, definition);
%
%     >> stages = [0, 0, ..., -1, -2, -3, ..., 1 ...};
%     >> definition = [0, -1, -2, -3, 1];
%     >> C = psg_sleepcycles(stages, definition);
%
% Authors: 
%     Rick Wassing, rickwassing@gmail.com
%     Yishul Wei,   ysw.1990@gmail.com
% History: 
%     05/03/2020 - Created

function C = psg_sleepcyles(stages, definition, varargin)

% check input
if ~any(size(stages)) == 1
    error('Input ''stage_num'' must be a vector')
end
if ~any(size(definition)) == 1
    error('Input ''definition'' must be a vector')
end
if length(definition) ~= 5
    error('Input ''definition'' must contain 5 elements. An integer or character that denotes each sleep stage (wake, N1, N2, N3, and REM).')
end
if ...
        (iscell(stages) && ~iscell(definition)) || ...
        (~iscell(stages) && iscell(definition)) || ...
        (isnumeric(stages) && ~isnumeric(definition)) || ...
        (~isnumeric(stages) && isnumeric(definition))
    error('Inputs ''stage_num'' and ''definition'' must be the same class.')
end

% force row vectors
if size(stages, 2)     == 1; stages     = stages';     end
if size(definition, 2) == 1; definition = definition'; end

% Parse the varargin
p = inputParser;
addParameter(p, 'RecStart', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', 'numel', 1}) ...
);
parse(p,varargin{:});
RecStart   = p.Results.RecStart;

% Save the time vector to the output
C.times = (0:30/(60*60*24):(length(stages)-1)*(30/(60*60*24))) + RecStart;

% translate the stage_num [vector] or {vector} to the cycle_num [vector]
stages_num = nan(1, length(stages));
if iscell(stages)
    stages_num(strcmp(stages, definition{1})) = 0; % Wake
    stages_num(strcmp(stages, definition{2})) = 1; % N1
    stages_num(strcmp(stages, definition{3})) = 2; % N2
    stages_num(strcmp(stages, definition{4})) = 3; % N3
    stages_num(strcmp(stages, definition{5})) = 5; % REM
elseif isnumeric(stages)
    stages_num(stages == definition(1)) = 0; % Wake
    stages_num(stages == definition(2)) = 1; % Wake
    stages_num(stages == definition(3)) = 2; % Wake
    stages_num(stages == definition(4)) = 3; % Wake
    stages_num(stages == definition(5)) = 5; % Wake
end
C.sleepstage = stages_num;

% Crop to the end of sleep, i.e. remove the final morning wakefulness period
% fin_awake = find(stages_num > 0,1,'last');
% stages_num = stages_num(1:fin_awake);

% remove periods of wake and N1 sleep directly after REM periods until N2 starts
stages_num(1) = nan;
for ep=2:numel(stages_num);
    if stages_num(ep)<2 && ~(stages_num(ep-1)<5);
        stages_num(ep)=nan;
    end
end

% when do the NREM cycles start and end
[~,cycle_beg,cycle_end] = get_bouts(~(stages_num < 5));

% if the cycle is less than 15 minutes, remove those periods.
for b=2:numel(cycle_beg);
    if cycle_beg(b)-cycle_end(b-1) <= 30; % 30 epochs = 15 minutes, NREM cycles must be at least 15 minutes
        stages_num((cycle_end(b-1)+1):(cycle_beg(b)-1)) = nan;
    end;
end
% and then find the start and ends again
[~,cycle_beg,cycle_end] = get_bouts(stages_num<5);

REM_dur = cycle_beg(2:end)-cycle_end(1:(end-1));
idx = REM_dur <= 10; % 10 epochs = 5 minutes, REM cycles must be at least 5 minutes
idx(1)=false; % except the first REM period which can be any duration.
cycle_beg([false idx])=[];
cycle_end([idx false])=[];

% The cycles should be at least 30 minutes
cysel = (cycle_end-cycle_beg+1)>29;
cycle_beg  = cycle_beg(cysel);
cycle_nend = cycle_end(cysel);
cycle_rend = [(cycle_beg(2:end)-1) numel(stages_num)];

C.sleepcycle   = nan(size(stages));
C.sleepepisode = nan(size(stages));

for c = 1:length(cycle_beg)
    C.sleepcycle(cycle_beg(c):cycle_rend(c)) = c;
    C.sleepepisode(cycle_nend(c)+1:cycle_rend(c)) = 1; % rem period
    
    % find the first two-consecutive N3 stages
    idx = find((stages_num(1:end-1) == 3 & C.sleepcycle(1:end-1) == c) & (stages_num(2:end) == 3 & C.sleepcycle(2:end) == c),1,'first');
    idx = [idx,find((stages_num(1:end-1) == 3 & C.sleepcycle(1:end-1) == c) & (stages_num(2:end) == 3 & C.sleepcycle(2:end) == c),1,'last')+1];
    
    if ~isempty(idx)
        C.sleepepisode(cycle_beg(c):idx(1)-1) = -1;
        C.sleepepisode(idx(1):idx(2)) = -2;
        C.sleepepisode(idx(2)+1:cycle_nend(c)) = 1;
    else
        % if SWS was not reached in this sleep cycle, find the first and last NREM stage
        idx = [find((stages_num == 1 | stages_num == 2 | stages_num == 3) & C.sleepcycle == c,1,'first'), ...
               find((stages_num == 1 | stages_num == 2 | stages_num == 3) & C.sleepcycle == c,1,'last')];
        C.sleepepisode(cycle_beg(c):idx(1)-1) = -1;
        C.sleepepisode(idx(1):idx(2)) = -2;
        C.sleepepisode(idx(2)+1:cycle_nend(c)) = 1;
    end
    C.sleepepisode(cycle_nend(c)+1:cycle_rend(c)) = 2; % rem period
end

end

% ------------
% Subfunctions
% ------------

function [b,ibeg,iend] = get_bouts(OnOffSeq)
    % Get bout lengths (in order) from a sequence of logicals
    beg_end = diff([0 OnOffSeq 0]);
    ibeg = sort(find(beg_end>0));
    iend = sort(find(beg_end<0));
    b = iend - ibeg;
    iend = iend - 1;
end
