clear; close all; clc;

% 音ファイルを読み込み
[src, fs] = audioread("pfcl.wav");  % 今回はsrc→66077サンプル，fs→44.1 [kHz]

% STFTをする
win_len = 4096;     % 窓長
shift_len = 1024;   % シフト長
STFT = DGTtool("windowLength", win_len, "windowShift", shift_len, "windowName", "hann");     % DGTtool使ってみた
X_comp = STFT(src);         % 複素スペクトログラム
X_obs = abs(X_comp);     % Xが観測信号の振幅スペクトログラム I×J
X = X_obs;              
phase = angle(X_comp);  % 位相を保存（istft前に使用）
I = size(X, 1);         % Xの行数
J = size(X, 2);         % Xの列数

% 乱数シード
seed = 123456;
rng(seed);

% WとHの定義
K = 2;      % 分解ランク（基底の数）．今回は2（2音源だから）
W = rand(I, K);
H = rand(K, J);

% 反復更新
max_itr = 100;           % 反復回数
epsilon = 10e-10;            % フロアリング用
loss = zeros(max_itr, 1);  % 損失保存用変数

for i = 1:max_itr
    W_new = W .* ((X_obs * H.') ./ (W * (H * (H.'))));     % Wの更新
    W = max(W_new, epsilon);            % フロアリング（小さい数になったりゼロ割が起きないように）
    H_new = H .* ((W.' * X_obs) ./ ((W.' * W) * H));       % Hの更新
    H = max(H_new, epsilon);            % フロアリング
    loss(i) = sum(abs(X_obs - W * H).^2, "all");       % 損失計算（目的関数）
end

plot(loss);       % 損失確認用

% Xを推定
X_est = W * H;          % 最終的なX
X_est_comp = X_est .* exp(phase * 1j);
est = STFT.pinv(X_est_comp);             % iSTFTで波形に戻す
audiowrite("est.wav", est, fs);     % 音を保存

% 個別の音を推定
X_est_A = W(:, 1) * H(1, :);
X_est_comp_A = X_est_A .* exp(phase * 1j);
est_A = STFT.pinv(X_est_comp_A);
audiowrite("est_scrA.wav", est_A, fs);

X_est_B = W(:, 2) * H(2, :);
X_est_comp_B = X_est_B .* exp(phase * 1j);
est_B = STFT.pinv(X_est_comp_B);
audiowrite("est_scrB.wav", est_B, fs);
