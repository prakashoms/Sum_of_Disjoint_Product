% function r = sumOfDisjointPdt(A, edgeVar)

% Code to compute r = P(A1 U A2 U A3 U ...)
% ------------------------------------------------------------------------
% For example
load('ammonia'); % loads collection of sets (*.mat file), e.g., A = {A1, A2, A3, ...}, where Ai = {...xi...}
edgeVar = {'x1','x2','x3','x4','x5','x6','x7','x8'}; % variable names to be used to generate expression

% ------------------------------------------------------------------------


num_var = length(edgeVar); % number of variables
m = length(A); % m = no. of rows in A, (A1, A2, ... Am)

r = []; % r is initialised. It contains the final expression
% The following loop generates the expression for first disjoint product(DP):
% The kind of expession being generated here is, P(A1) = a11*a12*a13*...a1n
A_1st_len = length(A{1,1});
for i = 1:A_1st_len-1
    var_name = edgeVar{A{1,1}(1,i)}; % Storing required elements from edgeVar
    r = strcat(r, var_name,'*');
end
r = strcat(r,edgeVar{A{1,1}(1,A_1st_len)});

r_star = [];
% The following loop generates the expression for prob P(A2)
% P(A1'A2) = P(A1'|A2)*P(A2); A1' is A1 bar(complement)
% P(A2): (r_star) prior prob = a21*a22*a23...a2n
A_star_len = length(A{2,1});
for i = 1:A_star_len-1
    var_name = edgeVar{A{2,1}(1,i)};
    r_star = strcat(r_star, var_name,'*');
end
r_star = strcat(r_star, edgeVar{A{2,1}(1,A_star_len)});

% In general, P(A1A2A3...Am) = P(A1'A2'A3'...Am-1'|Am)*P(Am). Joint
% probability of A1',A2'...Am-1' need to be evaluated given that Am is true.
% Set A1...Am-1 are represented by bar_grp and Am is rerepresented by
% A_star
% bar_grp is initialised . Here it contains binary representaion of set A1
bar_grp = zeros(1,num_var);
bar_grp(1,A{1,1}) = 1;
% A_star is initialised . Here it contains binary representaion of set A2
A_star = zeros(1,num_var);
A_star(1,A{2,1}) = 1;

% A_common is a vector containing colomn position of non-zero elements(=1)
A_common = find(A_star);
D_grp = bar_grp; % swapping of values. Will work with D_grp now

% For 2nd disjoint product (DP), which is
% P(A1'A2) = P(A1'|A2)*P(A2); A2 (r_star)  is already found above
% Here only one rule needs to be activated to solve P(A1'|A2)
% which is, Rule R0: (AB)'(AC) = (B)'(AC); (AB)' is AB whole complement

% Here R0 is activated. All non-zero(=1) positions in A_star is made '0' in
% bar_grp (=D_grp)
D_grp(:,A_common) = 0;

D_common = find(D_grp); % return colomn position of non-zero elements(=1) in D_grp
D_commom_len = length(D_common);
r_grp = [];
% The following loop generates the expression for prob P(A1'|A2)
for i = 1:D_commom_len - 1
    var_name = edgeVar{D_common(i)};
    r_grp = strcat(r_grp, var_name,'*');
end
r_grp = strcat(r_grp, edgeVar{D_common(D_commom_len)}); % print bar group in 2nd term

% r contains expression of first 2 DPs: P(A1)+P(A1'A2)
r = strcat(r, '+(1-', r_grp, ')*', r_star);

% In the following while loop, given expression will be generated:
% P(A1'A2'A3)+P(A1'A2'A3'A4)+...+P(A1'A2'A3'A4'...Am)

% fprintf(rId,'obj_exp = %s...\n',r); %%% xX
% % unrel_value(1) = steam_metering_unrel_exp1(xi);

k = 2; % k =2 bcoz first 2 DPs are evaluated. And it will run untill k<m
count_disjoint = 2; % bcz P(A1) and P(A1'|A2) are already calculated
r_term_length=0;
file_num=1;

while k ~= m
    %     tic
    r_term = [];
    
    % The following loop generates the expression for prior prob P(Ak+1)
    A_star_len = length(A{k+1,1});
    r_star = [];
    for i = 1:A_star_len-1
        var_name = edgeVar{A{k+1,1}(1,i)};
        r_star = strcat(r_star, var_name,'*');
    end
    r_star = strcat(r_star, edgeVar{A{k+1,1}(1,A_star_len)});
    
    % Here bar_grp will keep storing binary representaion of set Ak,
    % bar_grp will keep increasing as it is stacking set Ak and it will be
    % used in every iteration to find conditional prob, P(A1'A2'...Ak'|Ak+1)
    bar_grp(k,A{k,1}) = 1;
    % A_star contains binary representaion of set Ak+1
    A_star = zeros(1,num_var);
    A_star(1,A{k+1,1}) = 1;
    
    A_common = find(A_star); % find colomn positions of non-zero elements(=1)
    D_grp = bar_grp; % swapping of values. Will work with D_grp and store bar_grp for next iteration
    D_grp(:,A_common) = 0; % Here R0 is activated. All non-zero(=1) positions in A_star is made '0' in bar_grp (=D_grp)
    
    % Arrange all complemented (D_grp=bar_grp) group according to their
    % increasing order of cardinality. This is done for progressive
    % checking of superset
    [~, D_sort_idx] = sort(sum(D_grp,2)); % store index acoording to increasing order of cardinality
    D_grp = D_grp(D_sort_idx,:); % arranged complemented groups acoording to increasing order of cardinality
    D_grp_len = k; % size(D_grp,1); because k decides the row size of bar_grp; And size of bar_grp and D_grp is same
    
    % In the following while loop, each complemented group is checked
    % progressively with the rest of groups, for superset or equal set. If
    % such set is found, it will be deleted from D_grp; Rule is:
    % If (ABC) is superset of (AB):- (AB)'(ABC)'=(AB)'. OR,
    % if both row have same elements:- (AB)'(AB)'=(AB)';
    % e.g. Lets D_grp is a set of [11;12;13;14]. Then 11 will be checked
    % with 12, 13, 14; similarly 12 will be checked with 13 and 14.
    
    % Min. value of D_grp_len will be 2, as it is decided by k. w is
    % initialised as 1 to enter into the loop
    
    % % The following loop removes the superset only
    w = 1;
    while w < D_grp_len
        j = w + 1;
        % for progressively checking; used while loop bcz, as a superset or
        % equal set is found, it is deleted from D_grp. So its length may
        % change in each iteration
        while j <= D_grp_len
            D_grp_mulwj = sum((D_grp(w,:).* D_grp(j,:)),2);
            % after multiplying both row, only common pos will have 1;
            % And after adding in row, D_grp_mulwj will give no. of common elements
            D_grp_addw = sum(D_grp(w,:),2); % simply adding in row to get the no of elements in D_grp(w,:)
            
            if D_grp_mulwj == D_grp_addw % jth row is superset or equal set of wth row
                % If a set is a superset or an equal set, it will activate the
                % current if-condition
                D_grp(j,:) = []; % delete jth row; remove superset or equal set
                
                D_grp_len = D_grp_len - 1; % after deleting, length of D_grp_len is reduced by 1
                j = j - 1; % the pointer also needs to be reduced by 1
                
                % e.g. wrt above example, When 11 is being checked with 12.
                % And 12 is found superset, then 12 will be deleted from
                % D_grp. now, D_grp_len = 3. and j will be reduced to 1. So
                % that in next iteration, it will again start checking with
                % j = 2.
            end
            j = j + 1;
        end
        w = w + 1;
    end
    
    % D_grp_flag stores bianry no. '0' classifies that row in D_grp as
    % non-intersecting (independent) complemented group. '1' classifies
    % that row in D_grp as intersecting (dependent) complemented group. And
    % 0 classifies that row in D_grp as non-intersecting (independent)
    % complemented group.
    % So, further rules(R1,R2,R3) are needed only to simplify intersecting
    % comlemented groups. And non-itersecting groups will be kept aside
    D_grp_flag = zeros(1,D_grp_len); % D_grp_len is changed in the above loop
    
    % % The following while loop tag (=1) the dependent rows
    w = 1;
    while w < D_grp_len
        j = w + 1;
        % for progressively checking; used while. Here for loop can also be
        % used as there is no deletion
        while j <= D_grp_len
            D_grp_mulwj = sum((D_grp(w,:).* D_grp(j,:)),2);
            % after multiplying both row, only common pos will have 1;
            % And after adding in row, D_grp_mulwj will give no. of common elements
            D_grp_addw = sum(D_grp(w,:),2); % simply adding in row to get the no of elements in D_grp(w,:)
            
            if (D_grp_mulwj ~= 0) && (D_grp_mulwj < D_grp_addw)
                % (D_grp_mulwj ~= 0) => Sth is common in btw two sets and
                % (D_grp_mulwj < D_grp_addw), then only this conditon will
                % be activated
                % e.g. [12] and [13] will get in, means they are intersecting(dependent).
                % But [12] and [34] should not get in, as they are
                % non-intersecting (independent)
                D_grp_flag(1,w) = 1; % '1' categorises both wth and jth pos set as dependent groups
                D_grp_flag(1,j) = 1;
            end
            j = j + 1;
        end
        w = w + 1;
    end
    
    % The following loop separates intersecting (1=dependent) and
    % non-intersecting (0=independent) complemented groups
    dg = [];
    dg_len = 1;
    ig = [];
    ig_len = 1;
    for i = 1:D_grp_len
        if D_grp_flag(i) == 1
            dg(dg_len,:) = D_grp(i,:); % dependent
            dg_len = dg_len + 1;
        else
            ig(ig_len,:) = D_grp(i,:); % independent
            ig_len = ig_len + 1;
        end
    end
    % dg contains rows of dependent multiplicable terms and ig contains rows of independent
    % multiplicable terms which will be multiplied with disjoint terms(dsjnt) got by operations on dg ( +(.+.+dsjnt+.+.)*.*ig*.*rstar )
    dsjnt = [];
    
    if size(dg,1) ~= 0 % if dg has no elemnts then no need to simplify
        % simplifyDG.m takes input dg, which is a matrix. Each row contains
        % binary representaion of an intersecting (dependent) complemented
        % group. e.g. if dependent groups are (12)'(13)'(345)'. Then
        % dg=[1 1 0 0 0;1 0 1 0 0;0 0 1 1 1]
        % Output: dsjnt is a matrix where each row conatins disjoint terms
        % e.g. (12)'(13)'(345)'=(1)'(345)' + 1(2)'(3)'
        % dsjnt will be [2 0 4 4 4;-2 2 3 0 0]
        dsjnt = simplifyDG(dg);
    end
    count_disjoint = count_disjoint+size(dsjnt,1);
    if size(dsjnt,1) == 0
        count_disjoint = count_disjoint + 1; % for independent terms
    end
    % expressionGen generates expression. Here there can be three cases:
    % 1) ig=[],dsjnt ~=[]; 2) ig ~=[],dsjnt=[]; 3) ig ~=[],dsjnt ~=[]. Both
    % ig=[] and dsjnt=[] cannot happen.
    % For different cases, the expression generated will look like:
    % e.g., case 1: dsjnt = [0 0 1 0 0;0 2 -1 2 0]; r_str = ((1-x3)+x3*(1-x2*x4))
    % case 2: ig = [0 0 1 0 0;1 1 0 0 1]; r_str = (1-x3)*(1-x1*x2*x5)
    % case 3: dsjnt = [0 0 1 0 0;0 2 -1 2 0]; ig = [0 0 1 0 0;1 1 0 0 1]; then,
    % r_str = ((1-x3)+x3*(1-x2*x4))*(1-x3)*(1-x1*x2*x5)
    r_term = expressionGen_shorten(edgeVar,ig,dsjnt); % expressionGen_short: y = 1-x
    
    r_term_length= r_term_length + length(r_term);
    
    r = strcat(r, '+' , r_term , '*' , r_star);

    k = k+1
end