clear; clc;
%% SETUP
addpath(genpath('C:\Users\yoonm\OneDrive\Documents\MATLAB\RIR-Generator-master')); % Path to RIR-Generator toolbox
% Parameters
fs = 16000;
room_dim = [6 6 6];
rt60 = 0.3;
mic_spacing = 0.16;
mic_center = [3 3 1.75];
src_distance = 1.5;
snr_list = [5]; % SNR dB
c = 343;
n = 4096; % number of points for rir
mtype = 'bidirectional'; % mic type
diffuse = 1; 
%% File paths
%addpath('C:\Users\yoonm\OneDrive\Documents\ETSI Noise Handsfree'); % path to ETSI noise files
%addpath('C:\Users\yoonm\OneDrive\Documents\LibriSpeech\test-clean\121\123852') % path to LibriSpeech
clean_file = 'C:\Users\yoonm\OneDrive\Documents\MATLAB\Harvard_clean_16kHz.wav'; %
noise_file = 'C:\Users\yoonm\OneDrive\Documents\MATLAB\Pub_handsfree_16kHz.wav';

% Load signals
[clean, fs1] = audioread(clean_file);
[noise, fs2] = audioread(noise_file);  % 8-ch
if fs1 ~= fs2 || fs1 ~= fs || fs2 ~= fs
    disp("audio files must be in 16kHz SR.");
    if fs2 ~= fs
        noise = resample(noise,fs,fs2);
    end
    if fs1 ~= fs
        clean = resample(clean,fs,fs1);
    end
end

% Dipole mic positions (16 cm apart along x-axis)
mic1 = mic_center - [mic_spacing/2 0 0];
mic2 = mic_center + [mic_spacing/2 0 0];
mic_array = [mic1; mic2];

% Clean speech source (0° = front)
clean_src = mic_center + [0 src_distance 0];

% Angles for noise source (30° to 330°)
angles_deg = 90:90:270;
angles_rad = deg2rad(angles_deg);

% Extract channel 5 of ETSI noise
noise_ch5 = noise(:,5);
% noise_ch5 = noise_ch5(1:length(clean)); % Match length
%%
if diffuse == 1
    % Define fixed noise directions (diffuse)
    diffuse_angles_deg = [90, 180, 270];
    diffuse_angles_rad = deg2rad(diffuse_angles_deg);

    % Compute RIRs and convolve for each angle, then sum them
    noise_mic = zeros(length(noise_ch5) + n-1,2);  % Preallocate with margin

    for k = 1:length(diffuse_angles_rad)
        theta = diffuse_angles_rad(k);
        noise_src = mic_center + src_distance * [cos(theta), sin(theta), 0];
    
        rir_k = rir_generator(c, fs, mic_array, noise_src,room_dim, rt60, n, mtype);
    
        tmp = zeros(length(noise_ch5) + size(rir_k,2) - 1,2);
        tmp(:,1) = conv(noise_ch5, rir_k(1,:));
        tmp(:,2) = conv(noise_ch5, rir_k(2,:));
    
        % Sum contributions
        noise_mic(1:size(tmp,1),:) = noise_mic(1:size(tmp,1),:) + tmp;
    end

    % Convolve clean speech from the front (as before)
    rir_clean = rir_generator(c, fs, mic_array, clean_src,room_dim, rt60, n, mtype);
    % Convolve
    clean_mic = [conv(clean(:,1), rir_clean(1,:)) conv(clean(:,2), rir_clean(2,:))];
    

    % Truncate both to same length
    min_len = min(size(clean_mic,1), size(noise_mic,1));
    clean_mic = clean_mic(1:min_len, :);
    noise_mic = noise_mic(1:min_len, :);
    for snr_db = snr_list
        for ch = 1:2
            clean_ch = clean_mic(:, ch);
            noise_ch = noise_mic(:, ch);

            % Compute power
            clean_power = mean(clean_ch.^2);
            noise_power = mean(noise_ch.^2);

            % Scale noise for target SNR
            scaling_factor = sqrt(clean_power / (noise_power * 10^(snr_db/10)));
            noise_scaled = noise_ch * scaling_factor;

            % Mix
            mix(:, ch) = clean_ch + noise_scaled;
        end
        mix = mix / max(abs(mix), [], 'all') * 0.99;
        % Save to WAV
        out_name = sprintf('testvec_diffuse_snr%ddB.wav', snr_db);
        audiowrite(out_name, mix, fs);
    end
else
    for angle_idx = 1:length(angles_deg)
        theta = angles_rad(angle_idx);
        noise_src = mic_center + src_distance * [cos(theta) sin(theta) 0];

        % Generate RIRs
        rir_clean = rir_generator(c, fs, mic_array, clean_src, room_dim, rt60, n, mtype);
        rir_noise = rir_generator(c, fs, mic_array, noise_src, room_dim, rt60, n, mtype);

        % Convolve
        clean_mic = [conv(clean(:,1), rir_clean(1,:)) conv(clean(:,2), rir_clean(2,:))];
        noise_mic = [conv(noise_ch5, rir_noise(1,:)) conv(noise_ch5, rir_noise(2,:))];

        min_len = min(size(clean_mic,1), size(noise_mic,1));
        clean_mic = clean_mic(1:min_len,:);
        noise_mic = noise_mic(1:min_len,:);

        for snr_db = snr_list
            for ch = 1:2
                clean_ch = clean_mic(:, ch);
                noise_ch = noise_mic(:, ch);

                % Compute power
                clean_power = mean(clean_ch.^2);
                noise_power = mean(noise_ch.^2);

                % Scale noise for target SNR
                scaling_factor = sqrt(clean_power / (noise_power * 10^(snr_db/10)));
                noise_scaled = noise_ch * scaling_factor;

                % Mix
                mix(:, ch) = clean_ch + noise_scaled;
            end
            mix = mix / max(abs(mix), [], 'all') * 0.99; % Should I do this independently?
            % Save to WAV
            out_name = sprintf('testvec_angle%03d_snr%ddB.wav', angles_deg(angle_idx), snr_db);
            audiowrite(out_name, mix, fs);
        end
    end
end
disp('Test vector generation complete.');
