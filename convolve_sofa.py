#!/usr/bin/env python
# -*- coding: utf-8 -*-

import soundfile
import sofa
import numpy as np
from scipy import signal

def apply_IR(wavefile, sofafile, output, measurement, R_left=0, R_right=1, E=0):
    # HRTFデータの読み込み
    HRTF = sofa.Database.open(sofafile)

    IR_left = HRTF.Data.IR.get_values(indices={"M": measurement, "R": R_left, "E": E})
    IR_right = HRTF.Data.IR.get_values(indices={"M": measurement, "R": R_right, "E": E})

    # FFTに使用するため、要素数を2のべき乗に。(後ろを0埋め)
    #N = int(IR_right.size / 2)
    # overlap methodを使用するため、前半をさらに0埋め
    #N = 2**math.ceil(math.log2(data.size))
    #np.r_[np.zeros(N), data, np.zeros(N-data.size)]

    # WAVEファイルの読み込み
    wave_data, framerate = soundfile.read(wavefile)

    # Nごとに入力を区切る
    #frame_num = math.ceil((wave_data.size + N) / N)
    #frame_size = int(N * frame_num)
    # 最初のNフレームにIRを適用するため、N個の0を最初に追加。
    # 2^Nの余剰分だけ、後ろに0を追加
    #wave_data = np.r_[np.zeros(N), wave_data, np.zeros(frame_size - wave_data.size)]

    # convolution
    #IR_right_fft = np.fft.fft(IR_right, n=None)
    #IR_left_fft = np.fft.fft(IR_right, n=None)
    #out_right = np.array([])
    #out_left = np.array([])
    #
    #for i in range(frame_num):
    #    wave_fft_part = np.fft.fft(wave_data[i*N:(i+2)*N])
    #    convolution_right = wave_fft_part * IR_right_fft
    #    convolution_left = wave_fft_part * IR_left_fft
    #    out_right = np.r_[out_right, np.fft.ifft(convolution_right)[0:N].real]
    #    out_left = np.r_[out_left, np.fft.ifft(convolution_left)[0:N].real]

    out_right = signal.fftconvolve(wave_data, IR_right)
    out_left = signal.fftconvolve(wave_data, IR_left)

    # save wave file
    out_stereo = np.zeros((out_right.size, 2))
    out_stereo[:,0] = out_left
    out_stereo[:,1] = out_right
    soundfile.write(output+"_right.wav", out_right, framerate)
    soundfile.write(output+"_left.wav", out_left, framerate)
    soundfile.write(output+".wav", out_stereo, framerate)

if __name__ == '__main__':
    wavefile = "./sample.wav"
    #sofafile = "./HRTF-Database/SOFA/MRT02.sofa"
    sofafile = "./hrtf_nh2.sofa"

    for measurement in range(0, 1500, 100):
        apply_IR(
                wavefile=wavefile,
                sofafile=sofafile,
                output="./sample{0}".format(measurement),
                measurement=measurement,)
