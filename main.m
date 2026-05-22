% 音源ファイルの読み込み
[srcSig, fs] = audioread("pfcl.wav");

%stftに使う各種パラメータの設定
windowLength = 4096; % 窓長：100msに最も近い2のべき乗の点数(1s = fs個の点)
shiftLength = windowLength / 2; % シフト長：窓長の半分
windowName = "Hann"; %窓の種類：ハン窓

% 非負値行列の作成
F = DGTtool("windowShift", shiftLength, "windowLength", windowLength, "windowName", windowName);
X = F(srcSig); % stftの実行
phaseX = angle(X); % 複素スペクトログラムの位相を取得
X = abs(X); % 絶対値の取得

% 反復更新に用いる初期設定
K = 2; % 基底：K本の基底で混合信号を分解
t = 1000; % 反復回数
loss = zeros(1, t); % コスト関数
seed = 1;
rng(seed);
W = rand(size(X, 1), K); % 基底行列の初期化
H = rand(K, size(X, 2)); % 重み行列の初期化

% MMアルゴリズムに基づく反復更新測
for i = 1 : t
    W = W .* ((X * H.') ./ (W * (H * H.')));
    W = max(W, eps); % ゼロ割対策
    H = H .* ((W.' * X) ./ ((W.' * W) * H));
    H = max(H, eps); % ゼロ割対策
    loss(i) = sum((X - (W * H)).^2, 'all'); % コスト関数の計算
end

% 誤差の推移の表示
% title("誤差の推移");
% plot(loss);

% 複素スぺクトログラムの計算
nmfX = W * H; % 振幅スぺクトログラム
nmfX = nmfX .* exp(1i .* phaseX); % 複素スペクトログラム(for文で'i'を使っている場合虚数の'i'は'1i'を使用

% 複素スペクトログラムから元の音源を作成
nmfSig = F.pinv(nmfX);
% F.plot(nmfSig, fs); 
% sound(nmfSig, fs);

% 音源の出力
outputDir = "./output/"; % 出力先を指定
if ~exist(outputDir, 'dir') % 指定のフォルダがない場合作成
    mkdir(outputDir);
end

audiowrite(outputDir+"nmfSignal.wav", nmfSig, fs);


%追加課題
% 分離した音源の作成
sigLen = length(srcSig); % 信号源の長さを取得
estX = zeros(size(W, 1), size(H, 2), K);
estSig = zeros(sigLen, K);
for i = 1 : K % 規定数分繰り返し
    estX(:, :, i) = W(:, i) * H(i, :); % 振幅スぺクトログラム
    estX(:, :, i) = estX(:, :, i) .* exp(1i .* phaseX); % 複素スペクトログラム
    tmpSig = F.pinv(estX(:, :, i)); % 一時的に分離音源を保存（ゼロパディング付き）
    estSig(:, i) = tmpSig(1:sigLen); % ゼロパディングによるサイズの違いを考慮しながら分離した音源を求める
    audiowrite(sprintf(outputDir+"estimatedSignal%d.wav", i), estSig(:, i), fs); % 分離した音源
end