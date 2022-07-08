% applyaRule.m Last modifications: 16/06/18
function [pdt1, dg1] = applyaRule_mod(pdt, dg)
% applyaRule.m takes pdt and dg as input from simplifyTwoGroups.m and it
% returns pdt1 and dg1 just after applying one of these rules (R1,R2,R3)
% Priority at which rules are being applied:
% Rule 1: (A)'(AB)' = (A)'
% Rule 2: (A) (AB)' = (A)(B)'
% Rule 3: (AC)'(AB)' = (A)' + (A)(C)'(B)'
% ------------------------------------------------------------------------

edge_nm = size(dg,2); % no. of variables

mul_opr = pdt .* dg;
% if nothing is common in between pdt and dg (e.g. (AB)'(CD)'= (AB)'(CD)')
% Then the following codition will be active only
if sum(abs(mul_opr)) == 0
    pdt1 = zeros(1,edge_nm);
    dg1 = pdt + dg;
    return
end

% e.g. in pdt = ABC(DEF)'(GH)'(I)' and dg = (BCDEHJ)'; common non-bar part
% is (BC) and their pos is stored in pos_non_bar_part
pos_non_bar_part = find(mul_opr < 0);
non_bar_part = zeros(1,edge_nm);
non_bar_part(1,pos_non_bar_part) = pdt(1,pos_non_bar_part);

prirty_opt = 7; % randon initialization larger than 3 (bcz of prirty={1,2,3})
prirty = 5; % randon initialization smaller than prirty_opt but larger than 3
large_s1_len = 0;

% In the end of the following loop, priority at which a rule needs to be
% evaluated, is found.
for i = 1:edge_nm
    part_pdt1 = zeros(1,edge_nm);
    part_dg1 = zeros(1,edge_nm);
    pos_pdt = [];
    pos_dg = [];
    s2 = [];
    s2_len = 0;
    
    if mul_opr(i) < 0
        % mul_opr(i) cannot be = 0 for all i, bcz that case is already
        % considered in the above condition
        % e.g. Lets say, pdt = ABC(DEF)'(GH)'(I)' and dg = (BCDEHJ)'. Their
        % representation: pdt = [-2 -2 -3 1 1 1 3 3 4 0]; dg = [0 5 5 5 5 0 0 5 0 5];
        % After applying R2: (ABC(DEF)'(GH)'(I)') ((BCDEHJ)') = (ABC(DEF)'(GH)'(I)') (DEHJ)'
        % The common elements between non-bar and bar will be removed from
        % bar group
        
        % here first obj is to check whether cases like (AB)(AB)'=nul
        % exists
        pos_dg = find(dg == dg(i));
        part_dg1(1,pos_dg) = dg(1,pos_dg); % part of dg which is common with non-bar part is extracted
        mul_non_bar_dg1 = non_bar_part .* part_dg1;
        pos_mul_non_bar_dg1 = find(mul_non_bar_dg1 ~= 0); % pos of non-bar part which is common with extracted dg part
        
        % if length of extracted dg and non-bar part is equal, then
        % following case is active: (AB)(AB)'=nul;
        % e.g. Lets say, pdt = ABC, dg = (BC)'(DEF)'; Then extracted dg
        % part is (BC)' and extracted non-bar part is BC
        if length(pos_mul_non_bar_dg1) == length(pos_dg)
            pdt1 = zeros(1,edge_nm);
            dg1 = zeros(1,edge_nm);
            return
        end
        
        prirty = 2; % Rule 2: (A) (AB)' = (A)(B)'
        % positions in bar part common with non-bar part is stored, which
        % will be used to remove A from (AB)'
        s2 = pos_non_bar_part;
        
    elseif mul_opr(i) > 0
        % In both R1 and R3 multiplication at common pos will result in +ve
        % Here obj is to find whether R1 ( (A)'(AB)' = (A)') or R3 (
        % (AC)'(AB)' = (A)' + (A)(C)'(B)') will be active. In case of R1,
        % elements of the smaller group is subset of the other group. But
        % that is not the case in R3.
        
        pos_pdt = find(pdt == pdt(i));
        pos_dg = find(dg == dg(i));
        
        part_pdt1(1,pos_pdt) = pdt(1,pos_pdt); % bar-part of pdt which is common with bar-part of dg is extracted
        part_dg1(1,pos_dg) = dg(1,pos_dg);
        part_mult_opr = part_pdt1 .* part_dg1;
        s2 = find(part_mult_opr ~= 0);
        
        if (isempty( setdiff( pos_pdt, s2 ) ) == 1) ||...
                (isempty( setdiff( pos_dg, s2 ) ) == 1)
            % cases considered here are: 1. (A)'(A)'=A' , 2. (A)'(AB)'=(A)'
            % and 3. (AB)'(A)'=(A)'. And all the output will be part of pdt
            % Because in case of R1, (A)' and (AB)' can be part of either
            % pdt or dg
            
            prirty = 1; % Rule 1: (A)'(AB)' = (A)'
            uncommon_in_pdt = setdiff( pos_pdt, s2);
            % uncommon_in_pdt will give me the index of those elements
            % which are not coomon in both pdt and dg and exists in pdt.
            % e.g., pdt = (ABC)' , dg = (AB)'; therefore, pos_pdt=[1,2,3].
            % dg=[1,2], s2=[1,2]; uncommon_in_pdt=[3]. This will be used to
            % delete those elements from pdt.
            % But if pdt = (AB)' , dg = (ABC)'; then uncommon_in_pdt=[].
            % here there is no need to delete bcz it doesn't exist in pdt
            
        else
            prirty = 3; % Rule 3: (AC)'(AB)' = (A)' + (A)(C)'(B)'
        end
        
    end
    
    if prirty < prirty_opt % don't make less than equal to. Otherwise simplifyTwoGroups.m will go into an infinite loop
        prirty_opt = prirty;
        part_pdt1_opt = part_pdt1;
        part_dg1_opt = part_dg1;
        pos_pdt_opt = pos_pdt;
        pos_dg_opt = pos_dg;
        pos_part_cmn = s2; % position of common elements in part_pdt1 and part_dg1
        markr_pdt = pdt(i);
    end
    
    if prirty == 1
        % The moment the position where R1 needs to be applied is found,
        % it will break from the loop. There is no need to check further.
        break
    end
end

if prirty_opt == 1  % Rule 1: (A)'(AB)' = (A)'. It also takes care of (A)'(A)' = (A)'

    dg(1,pos_dg_opt) = 0;
    dg1 = dg;
    
    pdt(1,uncommon_in_pdt) = 0;
    pdt1 = pdt;

elseif prirty_opt == 2 % Rule 2: (A) (AB)' = (A)(B)'. It also takes care of (AB)(AB)'= nul
    % At any stage, pdt should have no elements in common. So whenever
    % splitting (R3) happens the branches created, will have common
    % elements (as in, A and A'). So it can be proved that A will always be
    % part of pdt and (AB)' be a part of dg. And the following works only
    % when that happen
    dg(1,pos_part_cmn) = 0;
    
    dg1 = dg;
    pdt1 = pdt;
    
elseif prirty_opt == 3 % Rule 3: (AC)'(AB)' = (A)' + (A)(C)'(B)'
    pdt(1,pos_pdt_opt) = 0;
    dg(1,pos_dg_opt) = 0;
    
    pdt1 = pdt + dg;
    
    lft_term = zeros(1,edge_nm); % lft_term contains elements in (A)'
    lft_term(1,pos_part_cmn) = markr_pdt; % replacing the common pos (for A') with their marker
    
    rt_term = part_pdt1_opt + part_dg1_opt; % rt_term contains elements in A(C)'(B)'
    rt_term(1,pos_part_cmn) = - markr_pdt; % replacing the common pos (for A) with their -ve marker
    
    dg1 = [lft_term ; rt_term];
    
end
