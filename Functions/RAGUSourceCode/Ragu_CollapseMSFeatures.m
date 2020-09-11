function out = Ragu_CollapseMSFeatures(in,BetweenDesign,WithinDesign)

[ngroups,nconds,~,~] = size(in);
out = nan(size(in));

if isempty(BetweenDesign)
    BetweenDesign = ones(ngroups,1);
end

if isempty(WithinDesign)
    WithinDesign = ones(nconds,1);
end

if ngroups ~= numel(BetweenDesign)
    error('Group sizes mismatch in Ragu_CollapseMSFeatures');
end

if nconds ~= numel(WithinDesign)
    error('Condition sizes mismatch in Ragu_CollapseMSFeatures');
end

nDifferentConditions = unique(WithinDesign);
nDifferentGroups     = unique(BetweenDesign);

for g = nDifferentGroups
    for c = nDifferentConditions
    gindex = find(BetweenDesign == g);
    cindex = find(WithinDesign  == c);
    m = mean(mean(in(gindex,cindex,:,:),1),2);
    out(gindex,cindex,:,:) = repmat(m,[numel(gindex),numel(cindex),1,1]);
   end
end


