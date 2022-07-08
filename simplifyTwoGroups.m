% simplifyTwoGroups.m Last modifications: 06/05/18
function dg = simplifyTwoGroups(pdt, dg)
% simplifyTwoGroups.m takes pdt and dg as input from simplifyDG.m and it
% returns dg which BECOMES, a row of disjoint terms after operations.
% NOTE: pdt and dg are rows from dsjnt and dg matrix of simplifyDG.m,
% respetively. So, pdt is a row of disjoint term and dg is a row
% complemented group.
% e.g. Lets say two groups are: (ABC(DEF)'(GH)'(I)') ((BCDEHJ)')
% Then, pdt = [-2 -2 -3 1 1 1 3 3 4 0]; dg = [0 5 5 5 5 0 0 5 0 5];

% simplifyTwoGroups.m only takes care of maintaining the tree. It picks two
% required goups and sends to applyaRule.m, which applies the desired rule
% and returns the simplified form.

% e.g.,      (ABC(DEF)'(GH)'(I)')  ((BCDEHJ)')
%                                      | (Rule, R2: applyaRule.m)
%             (ABC(DEF)'(GH)'(I)')  ((DEHJ)')
%                                      | (R3)
%                   ---------------------------------------
%                   |                                     |
% (ABC(GH)'(I)') [ (DE)'               +            (DE(F)'(HJ)') ]
%                   |                                     | (R3)
%       *1* ABC(GH)'(I)'(DE)'              ------------------------------
%                                          |                            |
%                       (ABCDE(F)'(I)') [ (H)'            +       (H(G)'(J)') ]
%                                          |                            |
%                               *2* ABCDE(F)'(H)'(I)'     +  *3* ABCDEH(F)'(G)'(I)'(J)'

% *1,2,3* groups (disjoint product) will form the final rows of dg at the end
% ------------------------------------------------------------------------

pdt_str = pdt;
row_dg = size(dg,1); % track the no. of dg grps (rows)
strt_of_branch = 1;
end_of_branch = row_dg;

% end_of_branch index should always be less than equal to the no. of
% row_dg(total row in current dg)
while end_of_branch <= row_dg
    k = 0;
    cumm_leaf_of_branch = 0;
    dg_group = [];
    for i = strt_of_branch : end_of_branch
        % strt_of_branch and end_of_branch are two extreme points of baby
        % branches of same parent branch.
        % e.g., In the above tree, in 3rd level, index of (DE)' is the
        % start_of_branch and index of (DE(F)'(HJ)') is end_of_branch
        
        % applyaRule.m is called to operate on pdt and dg
        % e.g., in 3rd level, pdt is (ABC(GH)'(I)') and dg is [(DE)' ; (DE(F)'(HJ)')]
        [pdt1,dg1] = applyaRule_mod(pdt, dg(i,:));
        
        dg_group = [dg_group; dg1];
        
        current_leaf_of_a_branch = size(dg1,1);
        % current_leaf_of_a_branch stores the current size of bifurcation
        % in dg after applying a rule.
        
        if sum(abs(pdt1)) ~= 0 % abs bcoz some may be -ve
            % this condition checks whether pdt1 is not empty. If it is not
            % empty, it will break from the loop and will enter into the
            % leaves (along depth) of current bifurcation and operate on
            % those before traversing along the width
            k = 1;
            % 0 or 1 based on whether this condition gets activated
            break
        end
        
        % cumm_leaf_of_branch gets calculated only if the above condition
        % doesn't get active. It stores the total bifurcation of all baby
        % branches of same parent branch. It will happen when pdt1 is empty
        % and there is nothing to get operated. e.g. in the last level
        cumm_leaf_of_branch = cumm_leaf_of_branch + current_leaf_of_a_branch;
        
        %         leaf_of_a_branch = size(dg1,1);
        %         cumm_leaf_of_branch = cumm_leaf_of_branch + leaf_of_a_branch;
    end
    
    dg = [dg(1:strt_of_branch-1,:) ; dg_group ; dg(i+1:end,:)];
    % dg1 will be inserted in btw before strt_of_branch and i (till i runs)
    
    row_dg = size(dg,1); % track the no. of dg grps (rows)
    
    pdt_str = [pdt_str ; pdt]; % stack the pdt
    
    if k == 0
        % implies it moved along the width (the above condition in the loop
        % doesn't get activated)
        % Here strt_of_branch and end_of_branch will be modified for the
        % next operation in the tree
        if i == end_of_branch
            pdt_str(end,:) = [];
            % if i hits the end_of_branch, which implies, operations are
            % done with all baby branches of a parent branch then it will
            % delete the last pdt stacked for that branch. Bcz that will no
            % longer be used
            % e.g. In 2nd level, star and end of a branch are both 1.
            % pdt_str will store (ABC(DEF)'(GH)'(I)'). Since i hits end_of
            % branch, (ABC(DEF)'(GH)'(I)') will be deleted from pdt_str
        end
        
        pdt = pdt_str(end,:);
        % pdt to be used in next iteration is modified
        strt_of_branch = strt_of_branch + cumm_leaf_of_branch;
        % index of strt_of_branch is modified as we are moving along the
        % width
        end_of_branch = strt_of_branch;
        % end_of_branch will be same as strt_of_branch as we are going to
        % operate on next one group in dg along the width
        
    else
        % k=1
        % If i hits the end of branch latest pdt stacked in pdt_str will be
        % deleted as it will no longer be used.
        % e.g. In 3rd level, as i hits end_of_branch (2:- (DE(F)'(HJ)')),
        % it gets bifurcated into (H)' and (H(G)'(J)') with pdt1 =
        % (ABCDE(F)'(I)'). So last pdt (= ABC(GH)'(I)') stored in pdt_str
        % will be deleted as it will no longer be used. And new pdt1 will
        % be stored into pdt for next iteration
        if i == end_of_branch
            pdt_str(end,:) = [];
        end
        
        pdt = pdt1;
        
        % while breaking out of the loop, i = strt_of _branch, then the new
        % str_of_branch will start from there only
        % if i~= strt_of_branch, implies i=end_of_branch, bcz here
        % there is at max two baby branches (R3) in a parent branch.
        % Then strt_of branch has to accomodate index of previous baby
        % branches created
        if i == strt_of_branch
            strt_of_branch = i;
        else
            strt_of_branch = (i + cumm_leaf_of_branch) - 1;
        end
        
        end_of_branch = (strt_of_branch + current_leaf_of_a_branch) - 1;
    end
end
% In case this rule gets activated: For, AB*(AB bar) = nul. applyaRule.m
% will return pdt1 = zeros(1,edge_nm); dg1 = zeros(1,edge_nm); So '0' vector
% will get stored in dg. This has to be removed from dg. If this will be
% done within loop, then it will create problem because of the way the
% index is managed. So, it is done outside the loop.

indx_0rows = find(sum(abs(dg),2)==0);

dg(indx_0rows,:) = [];

