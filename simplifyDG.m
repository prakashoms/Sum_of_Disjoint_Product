% simplifyDG.m Last modifications: 14/07/18
function dsjnt = simplifyDG(dg)
% simplifyDG.m takes input dg from sumOfDisjoint.m, which is a matrix. Each
% row contains binary representaion of an intersecting (dependent) complemented
% group. e.g. if dependent groups are (12)'(13)'(345)'. Then
% dg=[1 1 0 0 0;1 0 1 0 0;0 0 1 1 1]
% Output: dsjnt is a matrix where each row conatins disjoint terms
% e.g. (12)'(13)'(345)'=(1)'(345)' + 1(2)'(3)'
% dsjnt will be [1 0 3 3 3;-1 1 2 0 0]
% Here, all complemented (bar) part is represented by +ve integers and
% non-complemented (non-bar) part is represented by -ve integers
% ------------------------------------------------------------------------

[term_nm, edge_nm] = size(dg);
markr = 1:1:term_nm+1;

for i = 1:term_nm
    % All rows of dg has binary representaion. So, rows are multiplied with
    % markr to distinguish among each other. Like 1st, 2nd, 3rd rows are
    % multiplied with 1,2,3 respectively
    dg(i,:) = dg(i,:) * markr(i);
end

% The following loops find those two rows in dg which have maximum
% intersection (give their row pos, dg_i and dg_j). So, rule R3 will be
% applied first to these two rows. This way it minimizes the no. of
% disjoint terms in final exp; R3:  (AC)'(AB)' = (A)' + A(C)'(B)'
max_cmn = 0;
max_idx = [];
for i = 1:term_nm
    for j = i+1:term_nm
        mul_dg = dg(i,:) .* dg(j,:);
        pos_vec = find(mul_dg > 0);
        N = length(pos_vec);
        if N > max_cmn
            dg_i = i;
            dg_j = j;
            max_idx = pos_vec;
            max_cmn = N;
        end
    end
end

if max_cmn == 0 % if none of the rows in dg have any intersections
    dsjnt = sum(dg,1);
    display('simplifyDG:none of the rows in dg have any intersections')
    pause
    return % dg will always have sth common. This condition may not get activated ever
end

% R3:  (AC)'(AB)' = (A)' + A(C)'(B)'
lft_term = zeros(1,edge_nm); % lft_term contains elements in (A)'
lft_term(1,max_idx) = markr(dg_i); % replacing the common pos (for A') with their marker

rt_term = dg(dg_i,:) + dg(dg_j,:); % rt_term contains elements in A(C)'(B)'
rt_term(1,max_idx) = - markr(dg_i); % replacing the common pos (for A) with their -ve marker

dg([dg_i, dg_j],:) = []; % deleting rows from dg which are being operated
dsjnt = [lft_term ; rt_term]; % dsjnt contains the disjoint terms, generated after operation R3

rw_dsjnt = 2; % row size of dsjnt now

% The following loop picks one row from dg and operates with all rows of
% dsjnt. In each iteration two groups are send to simplifyTwoGroups.m, which
% returns disjoint terms. This continues untill all rows of dg are exausted

% term_nm signifies rows of dg. After above operations it is reduced by 2
for i = 1:term_nm-2 % rows of dg (terms in multiplied form)
    dsjnt_gen = [];
    for j = 1:rw_dsjnt % rows of dsjnt (terms in added form)
        
        simplify_exp = simplifyTwoGroups(dsjnt(j,:), dg(i,:));
        dsjnt_gen = [dsjnt_gen;simplify_exp];

    end
    dsjnt = dsjnt_gen;
    rw_dsjnt = size(dsjnt,1);
end

