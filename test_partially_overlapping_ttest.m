%% test_partially_overlapping_ttest.m
%
% Verification against the worked example in Derrick et al. (2017), Table 2-3.
%
% 16 students took one or both of two optional exam modules:
%   Mathematical Statistics (Sample 1): students 1-10
%   Operational Research (Sample 2): students 3-16
%   Students 3-10 took both (8 pairs)
%   Students 1-2 took only Math Stats (na = 2)
%   Students 11-16 took only Operational Research (nb = 6)
%
% Expected results (Table 3, p.145):
%   T_new2  = -2.276
%   v_new2  = 10.365
%   p-value = 0.045
%
% Also verifiable intermediate values (p.145):
%   xbar1 = 63.300,  xbar2 = 75.786
%   s1^2  = 263.789, s2^2  = 179.874
%   na = 2, nb = 6, nc = 8, n1 = 10, n2 = 14
%   r = 0.366
%   gamma = 17.095

%% Data from Table 2 (p.145)
%
% Student:                 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
% Mathematical Statistics: 73  82  74  59  49   -  42  71   -  39   -   -   -   -  59  85
% Operational Research:    72   -  89  78  64  83  42  76  79  89  67  82  85  92  63   -
%
% Sample 1 (Math Stats):      students 1,2,3,4,5,7,8,10,15,16  (n1=10)
% Sample 2 (Oper Research):   students 1,3,4,5,6,7,8,9,10,11,12,13,14,15  (n2=14)
% Paired (both modules):      students 1,3,4,5,7,8,10,15  (nc=8)
% Math Stats only:            students 2,16  (na=2)
% Oper Research only:         students 6,9,11,12,13,14  (nb=6)

x_A    = [73, 82, 74, 59, 49, 42, 71, 39, 59, 85]';
subj_A = {'1','2','3','4','5','7','8','10','15','16'}';

x_B    = [72, 89, 78, 64, 83, 42, 76, 79, 89, 67, 82, 85, 92, 63]';
subj_B = {'1','3','4','5','6','7','8','9','10','11','12','13','14','15'}';

%% Run our implementation

[p, stats] = partially_overlapping_ttest(x_A, x_B, subj_A, subj_B);

%% Display results and compare to expected values

fprintf('=== Verification of partially_overlapping_ttest.m ===\n');
fprintf('Data: Derrick et al. (2017), Table 2 worked example\n\n');

fprintf('--- Intermediate values (paper p.145) ---\n');
fprintf('  n1 (Math Stats):    %d      (expected: 10)\n', stats.na + stats.nc);
fprintf('  n2 (Oper Research): %d      (expected: 14)\n', stats.nb + stats.nc);
fprintf('  na (exclusive to 1): %d      (expected: 2)\n', stats.na);
fprintf('  nb (exclusive to 2): %d      (expected: 6)\n', stats.nb);
fprintf('  nc (pairs):          %d      (expected: 8)\n', stats.nc);
fprintf('  xbar1:             %7.3f  (expected: 63.300)\n', mean(x_A));
fprintf('  xbar2:             %7.3f  (expected: 75.786)\n', mean(x_B));
fprintf('  s1^2:              %7.3f  (expected: 263.789)\n', var(x_A));
fprintf('  s2^2:              %7.3f  (expected: 179.874)\n', var(x_B));
fprintf('  r:                 %7.3f  (expected: 0.366)\n', stats.r);

%% Compute gamma independently for verification
S1_sq = var(x_A); S2_sq = var(x_B);
n1 = length(x_A); n2 = length(x_B);
gamma = (S1_sq/n1 + S2_sq/n2)^2 / ((S1_sq/n1)^2/(n1-1) + (S2_sq/n2)^2/(n2-1));
fprintf('  gamma:             %7.3f  (expected: 17.095)\n', gamma);

fprintf('\n--- Test results (paper Table 3, p.145) ---\n');
fprintf('  T_new2:            %7.3f  (expected: -2.276)\n', stats.T);
fprintf('  v_new2 (df):       %7.3f  (expected: 10.365)\n', stats.df);
fprintf('  p-value:           %7.3f  (expected: 0.045)\n', p);
fprintf('  method:            %s\n', stats.method);

%% Pass/fail check
tol = 0.01;  % tolerance for rounding (paper reports 3 decimal places)
checks = [
    abs(stats.T - (-2.276)) < tol, ...
    abs(stats.df - 10.365) < tol, ...
    abs(p - 0.045) < tol, ...
    stats.nc == 8, ...
    stats.na == 2, ...
    stats.nb == 6
];
check_names = {'T_new2', 'v_new2', 'p-value', 'nc', 'na', 'nb'};

fprintf('\n--- Verification ---\n');
all_pass = true;
for i = 1:length(checks)
    if checks(i)
        fprintf('  %s: PASS\n', check_names{i});
    else
        fprintf('  %s: FAIL\n', check_names{i});
        all_pass = false;
    end
end

if all_pass
    fprintf('\nAll checks passed. Implementation matches Derrick et al. (2017).\n');
else
    fprintf('\nSome checks failed. Review implementation.\n');
end
