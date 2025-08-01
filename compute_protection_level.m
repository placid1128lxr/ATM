function [vpl, alloc] = compute_protection_level(sigma, threshold_plus_bias, pfault, phmi, pl_tol, alloc_max)

if nargin<6
    alloc_max = ones(length(sigma),1);
end

%%%% Exclude sigmas that are infinite and evaluate their integrity contribution
index_Inf = find(sigma == Inf);
index_fin = setdiff(1:length(sigma),index_Inf);
p_not_monitorable = sum(pfault(index_Inf));

if p_not_monitorable>=phmi    
    vpl = Inf;
    alloc = zeros(length(sigma),1);
else
    sigma = sigma(index_fin);
    threshold_plus_bias = threshold_plus_bias(index_fin);
    pfault = pfault(index_fin);
    phmi = phmi - p_not_monitorable;
    maxCount = 10;

    Klow=-norminv(min(1,phmi./(pfault.*alloc_max)));
    vpl_low=max(threshold_plus_bias + Klow.*sigma);

    Khigh = max(0,-norminv(phmi./(length(sigma)*pfault)));
    vpl_high=max(threshold_plus_bias + Khigh.*sigma);

    log10phmi=log10(phmi);

    count=0;
    while ((vpl_high-vpl_low>pl_tol)&&(count<maxCount))
        count = count+1;
        vpl_half = (vpl_low+vpl_high)/2;
        cdfhalf = log10(sum(pfault.*min(modnormcdf((threshold_plus_bias-vpl_half)./sigma),alloc_max)));
        if cdfhalf>log10phmi
           vpl_low = vpl_half;
        else
           vpl_high = vpl_half;
        end
    end

   vpl = vpl_high;
   alloc = min(modnormcdf((threshold_plus_bias-vpl)./sigma),alloc_max);
end
%End

