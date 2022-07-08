% expressionGen.m Last modifications: 26/05/18
function r_str = expressionGen_modified(edgeVar,ig,dsjnt)
% input: edgeVar: variable representaion, ig: rows contain independent
% terms, dsjnt: rows contain disjoint terms
% output: r_str: desired expression
% Here there can be three cases:
% 1) ig=[],dsjnt ~=[]; 2) ig ~=[],dsjnt=[]; 3) ig ~=[],dsjnt ~=[]. Both
% ig=[] and dsjnt=[] cannot happen.
% The way expressions will be generated in all these cases are as follow:
% if case 1 happens, then expressionGen.m will generate " (_*_*_ + _*_*_) "
% if case 2 happens, then expressionGen.m will generate " _*_*_*_ "
% if case 3 happens, then expressionGen.m will generate " (_*_*_ + _*_*_) *_*_*_*_ "
% e.g., case 1: dsjnt = [0 0 1 0 0;0 2 -1 2 0]; r_str = ((1-x3)+x3*(1-x2*x4))
% case 2: ig = [0 0 1 0 0;1 1 0 0 1]; r_str = (1-x3)*(1-x1*x2*x5)
% case 3: dsjnt = [0 0 1 0 0;0 2 -1 2 0]; ig = [0 0 1 0 0;1 1 0 0 1]; then,
% r_str = ((1-x3)+x3*(1-x2*x4))*(1-x3)*(1-x1*x2*x5)
% ------------------------------------------------------------------------

r_str = [];
terms_dsjnt = size(dsjnt,1);
terms_ig = size(ig,1);

for i = 1:terms_dsjnt
    % -ve numbers in dsjnt are non-bar part.
    % e.g. dsjnt = [-2 -2 0 1 0 -2]; so, r_nonbar = x1*x2*x6*
    pos_nonbar = find(dsjnt(i,:) < 0 ); % the negative elements are non-bar group
    pos_nonbar_len = length(pos_nonbar);
    r_nonbar = [];
    for j = 1:pos_nonbar_len
        var_name = edgeVar{pos_nonbar(j)};
        r_nonbar = strcat(r_nonbar, var_name, '*');
        % nonbar part will never occur alone. so the whole exp can be
        % computed in a loop. no need to put the last expression outside the
        % loop  a*b*...c*
    end
    
    % +ve numbers in dsjnt are bar part. In the next for loop following is
    % happening:
    % e.g. dsjnt = [3 0 3 -2 4 0 4]; so, r_bar = (1-x1*x3)*(1-x5*x7)*
    pos_bar = find(dsjnt(i,:) > 0 );
    pos_bar_len = length(pos_bar);
    r_bar = [];
    % the following loop will always run bcz dsjnt will always have sm bar
    % elements. So line77 r_bar(end)=[] will not cause problem bcz there
    % will always be * sign to remove.
    for j = 1:pos_bar_len
        
        if dsjnt(i,pos_bar(j)) == 0
            % In case, the element of dsjnt is 0, it will skip the rest and
            % bypass the current iteration
            continue
        end
        
        % e.g., Lets say, dsjnt = (AC)'D(EG)' = [3 0 3 -2 4 0 4]; So,
        % pos_bar_1group  will store position of all same +ve number to
        % compute the expression in the following for loop. For +ve no. 3,
        % it will be (1-x1*x3) and for 4 it will be (1-x5*x7)
        pos_bar_1group = find(dsjnt(i,:) == dsjnt(i, pos_bar(j)));
        pos_bar_1group_len = length(pos_bar_1group);
        
        % As the position of +ve no. is found, the no. from dsjnt is
        % deleted (=0), so that it will not be revisited, therefore there
        % is condition in line 48
        dsjnt(i, pos_bar_1group) = 0;
        
        r_grp = [];
        
        if pos_bar_1group_len == 1
            
            %             r_grp = strcat(r_grp, edgeVar{pos_bar_1group(pos_bar_1group_len)});
            
            r_bar = strcat(r_bar, 'y', int2str(pos_bar_1group(pos_bar_1group_len)), '*');
            
        else
            
            for k = 1:pos_bar_1group_len - 1
                var_name = edgeVar{pos_bar_1group(k)};
                r_grp = strcat(r_grp, var_name,'*');
            end
            r_grp = strcat(r_grp, edgeVar{pos_bar_1group(pos_bar_1group_len)});
            
            r_bar = strcat(r_bar, '(1-', r_grp, ')*');
            
        end
        
    end
    % r_bar looks like  (1-e*f*..g)*(1-h)* . so last * sign has to be
    % removed
    r_bar(end) = [];
    r_str = strcat(r_str, r_nonbar, r_bar, '+');
end

if terms_dsjnt ~= 0
    % In case 1 is active only
    % e.g., ig=[], dsjnt = [0 0 1 0 0;0 2 -1 2 0]; r_str = (1-x3)+x3*(1-x2*x4)+
    % so last + sign has to be removed
    % And finally, r_str = ((1-x3)+x3*(1-x2*x4))
    r_str(end) = [];
    r_str = strcat('(', r_str, ')');
end
% At the end of the previous loop, expression created will be like:
% Lets say, dsjnt = [0 0 1 0 0;0 2 -1 2 0]; r_str = ((1-x3)+x3*(1-x2*x4))

% -----------------------------evaluating ig-----------------------------
r_bar_ig_grp = [];
for i = 1:terms_ig
    % In ig, non-zero nos. are only '1'
    pos_bar_ig = find(ig(i,:) == 1 ); % '1's are bar group
    pos_bar_ig_len = length(pos_bar_ig);
    
    % The following expression will create expression which will look like:
    % Lets say, ig = [1 1 0 0 1]; r_grp = (1-x1*x2*x5)*
    r_grp = [];
    
    if pos_bar_ig_len == 1
        %         r_grp = strcat('(1-', r_grp, edgeVar{pos_bar_ig(pos_bar_ig_len)}, ')*');
        
        r_bar_ig_grp = strcat(r_bar_ig_grp, 'y', int2str(pos_bar_ig(pos_bar_ig_len)), '*');
        
    else
        for k = 1:pos_bar_ig_len - 1
            var_name = edgeVar{pos_bar_ig(k)};
            r_grp = strcat(r_grp, var_name,'*');
        end
        r_grp = strcat('(1-', r_grp, edgeVar{pos_bar_ig(pos_bar_ig_len)}, ')*');
        
        r_bar_ig_grp = strcat(r_bar_ig_grp, r_grp);
    end
    
    
end
% At the end of the previous loop, expression created will be like:
% Lets say, ig = [0 0 1 0 0;1 1 0 0 1]; r_str = (1-x3)*(1-x1*x2*x5)*

% Now trimming
if terms_dsjnt ~= 0 && terms_ig ~= 0
    % In case 3 is active, then remove * sign from end of r_bar_ig_grp first
    % e.g., dsjnt = [0 0 1 0 0;0 2 -1 2 0]; ig = [0 0 1 0 0;1 1 0 0 1]; then,
    % r_str = ((1-x3)+x3*(1-x2*x4)); r_bar_ig_grp = (1-x3)*(1-x1*x2*x5)*
    % Finally, r_str = ((1-x3)+x3*(1-x2*x4))*(1-x3)*(1-x1*x2*x5)
    r_bar_ig_grp(end) = [];
    r_str = strcat(r_str, '*' , r_bar_ig_grp);
    
elseif terms_ig ~= 0
    % In case 2 is active only, then remove * sign from end of r_bar_ig_grp first
    % e.g., dsjnt = []; ig = [0 0 1 0 0;1 1 0 0 1]; r_str = [],
    % r_bar_ig_grp = (1-x3)*(1-x1*x2*x5)*
    % And finally, r_str = (1-x3)*(1-x1*x2*x5)
    r_bar_ig_grp(end) = [];
    r_str = strcat(r_str, r_bar_ig_grp);
end
