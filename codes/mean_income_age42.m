function [ out ] = mean_income_age42( par,mu )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;
j_pos     = find(par.age == 42,1,'first');

S           = par.grids{1,j_pos};
inc_val     = zeros(length(S)*length(FE_pos)*length(INNO_pos),length(EDUC));
inc_prob    = zeros(length(S)*length(FE_pos)*length(INNO_pos),length(EDUC));

for educ = 1:length(EDUC)
    S          = par.grids{educ,j_pos};
    INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,1,length(INNO_pos)),length(S),length(FE_pos)));
    FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},1,length(FE_pos)),length(S),1,length(INNO_pos)));
    AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
    inc_val(:,educ)   = exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));

    prob_aux          = mu{j_pos}(:,:,educ,:,:);
    prob_aux            = squeeze(sum(prob_aux,5));
    inc_prob(:,educ)    = prob_aux(:);
end
inc_prob                = inc_prob./sum(inc_prob(:));

out              = inc_prob(:)'*inc_val(:);
    

end