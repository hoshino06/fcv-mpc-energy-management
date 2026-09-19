function [alpha_sol] = horizon_param(N,T_end,dt_first)

% --- 設定値 ---
% N = 5;               % ステップ数
% T_end = 0.5;         % ホライゾン全体の時間 [s]
% dt_first = 0.05;     % 最初のステップ幅 [s]

% --- 方程式の定義 ( f(alpha) = 0 となる alpha を探す ) ---
% 0.05 * sum(exp((0:N-1)/alpha)) - 0.5 = 0
f = @(alpha) dt_first * sum(exp((0:N-1)./alpha)) - T_end;

% --- 数値計算を実行 ---
% initial_guess = 3.0; % 適当な初期値
initial_guess = [0.01, 100];
alpha_sol = fzero(f, initial_guess);

fprintf('導出された alpha: %.4f\n', alpha_sol);

% --- 確認用のステップ幅計算 ---
sum_exp = sum(exp((1:N)./alpha_sol));
tsp_str = (exp((1:N)./alpha_sol) / sum_exp) * T_end;

disp('各ステップの幅 [s]:');
disp(tsp_str);
fprintf('合計時間: %.4f s (目標 %.4f)\n', sum(tsp_str), T_end);
fprintf('最初のステップ: %.4f s (目標 %.4f)\n', tsp_str(1), dt_first);

end