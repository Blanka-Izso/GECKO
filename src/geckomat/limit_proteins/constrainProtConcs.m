function model = constrainProtConcs(model)
% constrainProtConcs
%   Constrain enzyme usages by their concentration as provided in
%   model.ec.concs. For enzymes with non-NaN entries in model.ec.concs,
%   their protein pool draw reaction is replaced with a protein exchange
%   reaction, bound by the 
%
% Input:
%   model       an ec-model with enzyne levels in model.ec.concs
%
% Output:
%   model       an ec-model constraint with available protein levels
%
% Note: to populate model.ec.concs you should run getProteomics.s

%Enzyme with NaN entry in model.ec.concs => draw from prot_pool
%Enzyme with numeric entry in model.ec.concs => exchange reaction with
%enzyme level as UB

%Get indices of usage reactions 
usageRxns = strcat('usage_prot_',model.ec.enzymes);
[~, usageRxnsIdx] = ismember(usageRxns, model.rxns);

if any(usageRxnsIdx == 0)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end
%Get index of protein pool exchange rxn
protPoolIdx = ismember(model.rxns,'prot_pool_exchange');
if ~any(protPoolIdx)
    error('Cannot find protein pool exchange reaction.')
end

%Protein that should be constraint by UB
protCons = ~isnan(model.ec.concs);

%Set all reactions to draw from prot_pool
model.S(protPoolIdx, usageRxnsIdx(~protCons)) = -1;
model.ub(usageRxnsIdx(protCons)) = Inf;

%If non-NaN in model.ec.concs, then constrain by UB
if any(protCons)
    model.S(protPoolIdx, usageRxnsIdx(protCons)) = 0;
    model.ub(usageRxnsIdx(protCons)) = model.ec.concs(protCons);
end
end