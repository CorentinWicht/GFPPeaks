function [OutFeature,ndiff] = RaguMSFeatures(MSClass,n,MSFit,gfp,NoCenter)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

OutFeature = nan(size(MSClass,1),size(MSClass,2),n,6);
on   = nan(size(MSClass,1),size(MSClass,2),n);
off  = nan(size(MSClass,1),size(MSClass,2),n);
dur  = nan(size(MSClass,1),size(MSClass,2),n);
auc  = nan(size(MSClass,1),size(MSClass,2),n);
cog  = nan(size(MSClass,1),size(MSClass,2),n);
mgfp = nan(size(MSClass,1),size(MSClass,2),n);


MSClassR = reshape(MSClass,size(MSClass,1)*size(MSClass,2),size(MSClass,3));

%ndiff = sum(var(MSClassR)> 0);

h = hist(MSClassR,1:n) / size(MSClassR,1);

mscomplexity = h.*log(h);
mscomplexity(isnan(mscomplexity)) = 0;
mscomplexity = -sum(mscomplexity);
ndiff = sum(mscomplexity);

x = 1:size(MSClass,3);

for g = 1:size(MSClass,1)
    for c = 1:size(MSClass,2)
        for m = 1:n
            o = find((MSClass(g,c,:) == m),1,'first');
            f = find((MSClass(g,c,:) == m),1,'last');
            if ~isempty(gfp)
                mgfp(g,c,m) = squeeze(mean(gfp(g,c,1,MSClass(g,c,:) == m)));
            end
            dur(g,c,m) = sum(MSClass(g,c,:) == m);
            if ~isempty(o)
                on(g,c,m)  = o;
            end
            if ~isempty(f)
                off(g,c,m) = f;
            end
            MSGFP = squeeze(MSFit(g,c,:));
            MSGFP(MSClass(g,c,:) ~= m) = 0;
            auc(g,c,m) = sum(MSGFP);
            cog(g,c,m) = sum(MSGFP .* x') / auc(g,c,m);
%            auc(g,c,m) = sum(MSFit(g,c,MSClass(g,c,:) == m));
        end
    end
end

if nargin < 5
    NoCenter = true;
end

if NoCenter == true
    for m = 1:n
        on(:,:,m)  = on(:,:,m)  - squeeze(mean(mean(on(:,:,m),1),2));
        off(:,:,m) = off(:,:,m) - squeeze(mean(mean(off(:,:,m),1),2));
        dur(:,:,m) = dur(:,:,m) - squeeze(mean(mean(dur(:,:,m),1),2));
        auc(:,:,m) = auc(:,:,m) - squeeze(mean(mean(auc(:,:,m),1),2));
        cog(:,:,m) = cog(:,:,m) - squeeze(mean(mean(cog(:,:,m),1),2));
        mgfp(:,:,m) = mgfp(:,:,m) - squeeze(mean(mean(mgfp(:,:,m),1),2));
    end
end


OutFeature(:,:,:,1) = on;
OutFeature(:,:,:,2) = off;
OutFeature(:,:,:,3) = dur;
OutFeature(:,:,:,4) = auc;
OutFeature(:,:,:,5) = cog;
OutFeature(:,:,:,6) = mgfp;