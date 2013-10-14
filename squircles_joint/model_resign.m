
function sdata = model_resign(sdata)
    % resign means
    sdata.vb_mean   = sdata.vb_mean   .* sdata.sample_Vcat;
    sdata.sample_mv = sdata.sample_mv .* sdata.sample_Vcat;
end