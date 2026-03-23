# partially_overlapping_ttest

MATLAB implementation of the T_new2 test statistic from Derrick, Russ, Toher & White (2017) for comparing the means of two samples that include both paired and independent observations.

## The problem

Standard two-sample tests assume observations are either fully independent (e.g., Welch's t-test) or fully paired (e.g., paired t-test). In practice, data often fall between these extremes: some subjects appear in both groups (paired observations) while others appear in only one group (independent observations). This is known as **partially overlapping samples**.

Common examples:
- Longitudinal studies where population membership changes over time but some subjects are retained
- Two groups with a shared subset of participants (e.g., students taking one or both of two elective courses)
- Pre/post designs with dropout or new enrollment

Discarding the unpaired observations to run a paired test, or ignoring the pairing to run an independent test, both waste data and can distort Type I error rates and power.

## The solution

T_new2 uses **all available data** — both paired and independent observations — in a single test statistic. It is a Welch-like t-statistic whose standard error incorporates the covariance between paired observations. It reduces to:

- **Welch's t-test** when no subjects overlap (n_c = 0)
- **Paired t-test** when all subjects are paired (n_a = n_b = 0)

T_new2 has been validated by Monte Carlo simulation across 24,000+ parameter combinations and is Type I error robust under both equal and unequal variances (Derrick et al., 2017, Figures 1-2, Tables 3-6).

## Usage

```matlab
[p, stats] = partially_overlapping_ttest(x_A, x_B, subj_A, subj_B)
```

### Inputs

| Parameter | Description |
|-----------|-------------|
| `x_A` | Numeric vector of observations from group A (one per subject) |
| `x_B` | Numeric vector of observations from group B (one per subject) |
| `subj_A` | Cell array of subject identifiers for group A |
| `subj_B` | Cell array of subject identifiers for group B |

Subject identifiers are used to determine which observations are paired (same subject in both groups) and which are independent (subject in one group only). Each subject should appear at most once per group.

### Outputs

| Field | Description |
|-------|-------------|
| `p` | Two-sided p-value from the t-distribution |
| `stats.T` | T_new2 test statistic |
| `stats.df` | Approximate degrees of freedom (v_new2) |
| `stats.nc` | Number of paired subjects |
| `stats.na` | Number of subjects exclusive to group A |
| `stats.nb` | Number of subjects exclusive to group B |
| `stats.r` | Pearson correlation between paired observations |
| `stats.method` | `'T_new2'`, `'welch'` (fallback), or `'insufficient'` |

### Example

```matlab
% 5 subjects in group A, 4 in group B, subjects 3 and 4 are in both
x_A = [2.1, 3.4, 5.0, 4.2, 1.8]';
x_B = [6.3, 5.8, 7.1, 6.0]';
subj_A = {'s1','s2','s3','s4','s5'}';
subj_B = {'s3','s4','s6','s7'}';

[p, stats] = partially_overlapping_ttest(x_A, x_B, subj_A, subj_B);
fprintf('T=%.3f, df=%.1f, p=%.4f (nc=%d, na=%d, nb=%d)\n', ...
    stats.T, stats.df, p, stats.nc, stats.na, stats.nb);
```

## Verification

`test_partially_overlapping_ttest.m` reproduces the worked example from Table 2-3 of Derrick et al. (2017, p. 145). Expected output:

```
T_new2:  -2.276  (expected: -2.276)
v_new2:  10.365  (expected: 10.365)
p-value:  0.045  (expected: 0.045)
```

The R package [`Partiallyoverlapping`](https://cran.r-project.org/package=Partiallyoverlapping) provides an independent implementation for cross-validation.

## When to use subject-level aggregation

If subjects contribute multiple observations within a group (e.g., multiple recording sessions), reduce to one value per subject per group before calling this function. Taking the within-subject mean is a standard approach for eliminating intracluster correlation (Galbraith, Daniel & Vissel, 2010, *J Neuroscience*).

## Reference

Derrick, B., Russ, B., Toher, D., & White, P. (2017). Test statistics for the comparison of means for two samples that include both paired and independent observations. *Journal of Modern Applied Statistical Methods*, 16(1), 137-157. doi: [10.22237/jmasm/1493597280](https://doi.org/10.22237/jmasm/1493597280)

## License

MIT
