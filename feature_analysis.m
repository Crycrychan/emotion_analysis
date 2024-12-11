clearvars;
RMSSD = [];
ECG_LF_PSD = [];
ECG_HF_PSD = [];
ECG_LH_RATIO = [];
EEG_THETA = [];
EEG_ALPHA = [];
EEG_BETA = [];
relaxed_ecg = [];
relaxed_eeg = [];
relaxed_feature = table(RMSSD, ECG_LF_PSD, ECG_HF_PSD, ECG_LH_RATIO, EEG_THETA, EEG_ALPHA, EEG_BETA);
happy_ecg = [];
happy_eeg = [];
happy_feature = table(RMSSD, ECG_LF_PSD, ECG_HF_PSD, ECG_LH_RATIO, EEG_THETA, EEG_ALPHA, EEG_BETA);
annoyed_ecg = [];
annoyed_eeg = [];
annoyed_feature = table(RMSSD, ECG_LF_PSD, ECG_HF_PSD, ECG_LH_RATIO, EEG_THETA, EEG_ALPHA, EEG_BETA);
read = readmatrix("./data/Andy_relaxed_1.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_relaxed_2.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_relaxed_3.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_relaxed_4.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_relaxed_5.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_relaxed_1.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_relaxed_2.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_relaxed_3.csv");
relaxed_ecg = [relaxed_ecg; transpose(read(:, 3))]; relaxed_eeg = [relaxed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_happy_1.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_happy_2.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_happy_3.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_happy_4.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_happy_5.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_happy_1.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_happy_2.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_happy_3.csv");
happy_ecg = [happy_ecg; transpose(read(:, 3))]; happy_eeg = [happy_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_annoyed_1.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_annoyed_2.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_annoyed_3.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_annoyed_4.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/Andy_annoyed_5.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_annoyed_1.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_annoyed_2.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];
read = readmatrix("./data/YEH_annoyed_3.csv");
annoyed_ecg = [annoyed_ecg; transpose(read(:, 3))]; annoyed_eeg = [annoyed_eeg; transpose(read(:, 2))];

fs = 2000;
L = 60001;
entries = 8;

for i = 1:entries
    relaxed_ecg(i, :) = relaxed_ecg(i, :) - mean(relaxed_ecg(i, :));
    [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(relaxed_ecg(i, :), fs, 0);
    rmssd = rms(diff(qrs_i_raw(2:end-1)) .* 0.0005);
    [pxx, f] = pwelch(relaxed_ecg(i, :), [], [], [], fs);
    ecg_lf_psd = sum(pxx .* (f > 0) .* (f < 5)) * (f(2) - f(1));
    ecg_hf_psd = sum(pxx .* (f > 5) .* (f < 30)) * (f(2) - f(1));
    [pxx, f] = pwelch(relaxed_eeg(i, :), [], [], [], fs);
    eeg_theta = sum(pxx .* (f > 4) .* (f < 8)) * (f(2) - f(1));
    eeg_alpha = sum(pxx .* (f > 8) .* (f < 13)) * (f(2) - f(1));
    eeg_beta = sum(pxx .* (f > 13) .* (f < 20)) * (f(2) - f(1));
    relaxed_feature(end+1, :) = {rmssd, ecg_lf_psd, ecg_hf_psd, ecg_hf_psd / ecg_lf_psd, eeg_theta, eeg_alpha, eeg_beta};

    happy_ecg(i, :) = happy_ecg(i, :) - mean(happy_ecg(i, :));
    [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(happy_ecg(i, :), fs, 0);
    rmssd = rms(diff(qrs_i_raw(2:end-1)) .* 0.0005);
    [pxx, f] = pwelch(happy_ecg(i, :), [], [], [], fs);
    ecg_lf_psd = sum(pxx .* (f > 0) .* (f < 5)) * (f(2) - f(1));
    ecg_hf_psd = sum(pxx .* (f > 5) .* (f < 30)) * (f(2) - f(1));
    [pxx, f] = pwelch(happy_eeg(i, :), [], [], [], fs);
    eeg_theta = sum(pxx .* (f > 4) .* (f < 8)) * (f(2) - f(1));
    eeg_alpha = sum(pxx .* (f > 8) .* (f < 13)) * (f(2) - f(1));
    eeg_beta = sum(pxx .* (f > 13) .* (f < 20)) * (f(2) - f(1));
    happy_feature(end+1, :) = {rmssd, ecg_lf_psd, ecg_hf_psd, ecg_hf_psd / ecg_lf_psd, eeg_theta, eeg_alpha, eeg_beta};

    annoyed_ecg(i, :) = annoyed_ecg(i, :) - mean(annoyed_ecg(i, :));
    [qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(annoyed_ecg(i, :), fs, 0);
    rmssd = rms(diff(qrs_i_raw(2:end-1)) .* 0.0005);
    [pxx, f] = pwelch(annoyed_ecg(i, :), [], [], [], fs);
    ecg_lf_psd = sum(pxx .* (f > 0) .* (f < 5)) * (f(2) - f(1));
    ecg_hf_psd = sum(pxx .* (f > 5) .* (f < 30)) * (f(2) - f(1));
    [pxx, f] = pwelch(annoyed_eeg(i, :), [], [], [], fs);
    eeg_theta = sum(pxx .* (f > 4) .* (f < 8)) * (f(2) - f(1));
    eeg_alpha = sum(pxx .* (f > 8) .* (f < 13)) * (f(2) - f(1));
    eeg_beta = sum(pxx .* (f > 13) .* (f < 20)) * (f(2) - f(1));
    annoyed_feature(end+1, :) = {rmssd, ecg_lf_psd, ecg_hf_psd, ecg_hf_psd / ecg_lf_psd, eeg_theta, eeg_alpha, eeg_beta};
end

%% 
corr_coef = [];
for i = 1:7
    corr_coef = [corr_coef 0];
    x = [ones(entries, 1); -1.*ones(entries, 1)];
    y = [relaxed_feature{:, i}; happy_feature{:, i}];
    temp = corrcoef(x, y);
    if(abs(temp(2)) > abs(corr_coef(i)))
        corr_coef(i) = temp(2);
    end
    y = [relaxed_feature{:, i}; annoyed_feature{:, i}];
    temp = corrcoef(x, y);
    if(abs(temp(2)) > abs(corr_coef(i)))
        corr_coef(i) = temp(2);
    end
    y = [annoyed_feature{:, i}; happy_feature{:, i}];
    temp = corrcoef(x, y);
    if(abs(temp(2)) > abs(corr_coef(i)))
        corr_coef(i) = temp(2);
    end
end
temp = corr_coef;
corr_coef = table(RMSSD, ECG_LF_PSD, ECG_HF_PSD, ECG_LH_RATIO, EEG_THETA, EEG_ALPHA, EEG_BETA);
corr_coef = [corr_coef; array2table(temp, 'VariableNames', corr_coef.Properties.VariableNames)];
