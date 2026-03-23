function [p, stats] = partially_overlapping_ttest(x_A, x_B, subj_A, subj_B)
%PARTIALLY_OVERLAPPING_TTEST T_new2 test for partially overlapping samples
%
%   [p, stats] = partially_overlapping_ttest(x_A, x_B, subj_A, subj_B)
%
%   Implements the T_new2 statistic from Derrick, Russ, Toher & White
%   (2017) for comparing two groups where some subjects contribute
%   observations to both groups (partially overlapping samples). T_new2 is
%   a Welch-like test that uses all available data, accounting for the
%   covariance between paired observations.
%
%   When no subjects overlap (nc = 0), reduces to Welch's t-test.
%   When all subjects are paired (na = nb = 0), reduces to paired t-test.
%
%   Expects one observation per subject per group (apply subject-level
%   aggregation first if subjects have multiple observations per group).
%
%   Inputs:
%     x_A    - [nA x 1] observations from group A (one per subject)
%     x_B    - [nB x 1] observations from group B (one per subject)
%     subj_A - {nA x 1} unique subject IDs for group A
%     subj_B - {nB x 1} unique subject IDs for group B
%
%   Outputs:
%     p     - two-sided p-value (t-distribution)
%     stats - struct with fields:
%       .T       - T_new2 test statistic
%       .df      - approximate degrees of freedom (v_new2)
%       .nc      - number of paired subjects (in both groups)
%       .na      - subjects exclusive to group A
%       .nb      - subjects exclusive to group B
%       .r       - Pearson correlation between paired observations
%       .method  - 'T_new2' | 'welch' | 'insufficient'
%
%   Reference:
%     Derrick, B., Russ, B., Toher, D., & White, P. (2017). Test
%     statistics for the comparison of means for two samples that include
%     both paired and independent observations. Journal of Modern Applied
%     Statistical Methods, 16(1), 137-157. doi: 10.22237/jmasm/1493597280

    x_A = x_A(:); x_B = x_B(:);
    subj_A = subj_A(:); subj_B = subj_B(:);

    % Remove NaN observations
    ok = ~isnan(x_A); x_A = x_A(ok); subj_A = subj_A(ok);
    ok = ~isnan(x_B); x_B = x_B(ok); subj_B = subj_B(ok);

    n1 = length(x_A);  % total observations in Sample 1
    n2 = length(x_B);  % total observations in Sample 2

    if n1 <= 1 || n2 <= 1
        p = NaN;
        stats = struct('T', NaN, 'df', NaN, 'nc', 0, ...
            'na', n1, 'nb', n2, 'r', NaN, 'method', 'insufficient');
        return;
    end

    % Identify paired subjects (appear in both groups)
    [is_paired_A, idx_in_B] = ismember(subj_A, subj_B);

    nc = sum(is_paired_A);   % number of pairs
    na = n1 - nc;            % observations exclusive to Sample 1
    nb = n2 - nc;            % observations exclusive to Sample 2

    % Sample means — all observations in each sample (Eq. Table 1)
    xbar1 = mean(x_A);
    xbar2 = mean(x_B);

    % Sample variances — Bessel-corrected (Eq. Table 1)
    S1_sq = var(x_A);   % uses (n1-1) denominator
    S2_sq = var(x_B);
    S1 = sqrt(S1_sq);
    S2 = sqrt(S2_sq);

    % Pearson correlation between paired observations
    if nc >= 2
        paired_xA = x_A(is_paired_A);
        paired_xB = x_B(idx_in_B(is_paired_A));
        r = corr(paired_xA, paired_xB);
    else
        r = 0;  % cannot estimate correlation with < 2 pairs
    end

    % Handle zero-variance edge case
    if S1_sq + S2_sq == 0
        p = NaN;
        stats = struct('T', NaN, 'df', NaN, 'nc', nc, ...
            'na', na, 'nb', nb, 'r', r, 'method', 'zero_variance');
        return;
    end

    % --- T_new2 statistic (Derrick et al. 2017, p.144) ---
    %
    %                    X_bar1 - X_bar2
    %  T_new2 = ------------------------------------
    %           sqrt( S1^2/n1 + S2^2/n2 - 2r*S1*S2*nc/(n1*n2) )

    denom_sq = S1_sq/n1 + S2_sq/n2 - 2*r*S1*S2*nc/(n1*n2);

    if denom_sq <= 0
        % Can occur with very high positive correlation; fall back to
        % Welch's test (ignoring pairing) as a conservative default
        denom_sq = S1_sq/n1 + S2_sq/n2;
        method = 'welch';
    else
        method = 'T_new2';
    end

    T = (xbar1 - xbar2) / sqrt(denom_sq);

    % --- Degrees of freedom: v_new2 (Derrick et al. 2017, p.144) ---
    %
    % gamma = Welch-Satterthwaite df (same as Welch's test)
    % v_new2 interpolates between v1 = nc-1 and gamma:
    %
    %   v_new2 = (nc-1) + (gamma - (nc-1)) / (na+nb+2*nc) * (na+nb)
    %
    % Defaults: na=nb=0 -> v1 (paired);  nc=0 -> gamma (Welch)

    gamma = (S1_sq/n1 + S2_sq/n2)^2 / ...
            ((S1_sq/n1)^2/(n1-1) + (S2_sq/n2)^2/(n2-1));

    if nc > 0 && (na + nb) > 0
        df = (nc - 1) + ((gamma - (nc - 1)) / (na + nb + 2*nc)) * (na + nb);
    elseif nc > 0
        df = nc - 1;        % all paired
    else
        df = gamma;          % all independent (Welch)
    end

    df = max(df, 1);  % ensure df >= 1

    % Two-sided p-value from t-distribution
    p = 2 * tcdf(-abs(T), df);

    stats.T = T;
    stats.df = df;
    stats.nc = nc;
    stats.na = na;
    stats.nb = nb;
    stats.r = r;
    stats.method = method;
end
