#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import soundfile
import sofa
from scipy import signal

def _apply_IR(wavefile, HRTF, output, measurement, R_left=0, R_right=1, E=0):
    """
    wavefileにsofafileのインパルス応答を畳み込む.
    scipyを使わず、numpyで自分で畳み込み積分を実装。基本的に使わない。

    wavefile: 入力の音源ファイル
    sofafile: 入力の.sofaファイル
    output: 出力の音源ファイル
    measurement: sofafileに含まれるインパルス応答のうち、使用する観測番号
    R_left, R_right, E: 詳細はsofa.org参照
    """
    IR_left = HRTF.Data.IR.get_values(indices={"M": measurement, "R": R_left, "E": E})
    IR_right = HRTF.Data.IR.get_values(indices={"M": measurement, "R": R_right, "E": E})

    # FFTに使用するため、要素数を2のべき乗に。(後ろを0埋め)
    N = int(IR_right.size / 2)
    # overlap methodを使用するため、前半をさらに0埋め
    N = 2**math.ceil(math.log2(data.size))
    np.r_[np.zeros(N), data, np.zeros(N-data.size)]

    # WAVEファイルの読み込み
    wave, framerate = soundfile.read(wavefile)

    # Nごとに入力を区切る
    frame_num = math.ceil((wave.size + N) / N)
    frame_size = int(N * frame_num)
    # 最初のNフレームにIRを適用するため、N個の0を最初に追加。
    # 2^Nの余剰分だけ、後ろに0を追加
    wave = np.r_[np.zeros(N), wave, np.zeros(frame_size - wave.size)]

    # convolution
    IR_right_fft = np.fft.fft(IR_right, n=None)
    IR_left_fft = np.fft.fft(IR_right, n=None)
    wave_conv_right = np.array([])
    wave_conv_left = np.array([])

    for i in range(frame_num):
        wave_fft_part = np.fft.fft(wave[i*N:(i+2)*N])
        convolution_right = wave_fft_part * IR_right_fft
        convolution_left = wave_fft_part * IR_left_fft
        wave_conv_right = np.r_[wave_conv_right, np.fft.ifft(convolution_right)[0:N].real]
        wave_conv_left = np.r_[wave_conv_left, np.fft.ifft(convolution_left)[0:N].real]

    #wave_conv_right = signal.fftconvolve(wave, IR_right)
    #wave_conv_left = signal.fftconvolve(wave, IR_left)

    # save wave file
    out_stereo = np.zeros((wave_conv_right.size, 2))
    out_stereo[:,0] = wave_conv_left
    out_stereo[:,1] = wave_conv_right
    #soundfile.write(output+"_right.wav", wave_conv_right, framerate)
    #soundfile.write(output+"_left.wav", wave_conv_left, framerate)
    soundfile.write(output+".wav", out_stereo, framerate)


def apply_IR(wavefile, HRTF, output, measurement, R_left=0, R_right=1, E=0):
    """
    wavefileにsofafileのインパルス応答を畳み込む

    wavefile: 入力の音源ファイル
    sofafile: 入力の.sofaファイル
    output: 出力の音源ファイル
    measurement: sofafileに含まれるインパルス応答のうち、使用する観測番号
    R_left, R_right: 左耳、右耳のreceiverのindex。通常は、左耳が0, 右耳が1
    E: emitterのindex
    """
    IR_left = HRTF.Data.IR.get_values(indices={"M": measurement, "R": R_left, "E": E})
    IR_right = HRTF.Data.IR.get_values(indices={"M": measurement, "R": R_right, "E": E})

    # WAVEファイルの読み込み
    wave, framerate = soundfile.read(wavefile)

    # waveファイルがステレオの場合
    if wave.ndim > 1:
        # 左音源を使用
        # TODO: 左右音源を考慮したconvolutionの実装(何かしらの仮定が必要?)
        wave = wave[:, 0] # 左音源

    wave_conv_right = signal.fftconvolve(wave, IR_right)
    wave_conv_left = signal.fftconvolve(wave, IR_left)

    # save wave file
    out_stereo = np.zeros((wave_conv_right.size, 2))
    out_stereo[:,0] = wave_conv_left
    out_stereo[:,1] = wave_conv_right
    soundfile.write(output+".wav", out_stereo, framerate)

def plot_source_positions(HRTF, output):
    """
    system = spherical: [az (deg), el (deg), r (m)] 
    system = cartesian: [x, y, z] (m)

    0 < azimuth < 2pi
    az = 0 -> along x-axis

    -pi/2 < el < pi/2
    el = 0 -> z = 0 horizon

    x = r*cos(el)*cos(az)
    y = r*cos(el)*sin(az)
    z = r*sin(el)
    """
    source_positions = HRTF.Source.Position.get_values(system="cartesian")

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(source_positions[:,0], source_positions[:,1], source_positions[:,2], "o", ms=2)
    fig.savefig(output)

def find_measurement(HRTF, az, el):
    """
    azimuth, elevationで指定された方向に最も近い観測のmeasurement indexを返す
    (az, el方向の単位ベクトルとの内積が最も大きいindexを返す)

    azimuth, elevetation: degree
    """
    az *= np.pi / 180.
    el *= np.pi / 180.
    unit_vector = np.array([
            np.cos(el) * np.cos(az),
            np.cos(el) * np.sin(az),
            np.sin(el),
            ])

    measurements = HRTF.Source.Position.get_values(system="cartesian")
    return np.argmax(np.dot(measurements, unit_vector)) # np.dotで内積を計算し、その最大値のindexを返す


if __name__ == '__main__':
    wavefile = "./sample.wav"
    sofafile = "./hrtf_nh2.sofa"

    HRTF = sofa.Database.open(sofafile)
    #plot_source_positions(HRTF, "test.png")
    measurements = []
    measurements.append(find_measurement(HRTF=HRTF, az=0., el=0.))
    measurements.append(find_measurement(HRTF=HRTF, az=90., el=0.))
    measurements.append(find_measurement(HRTF=HRTF, az=180., el=0.))
    measurements.append(find_measurement(HRTF=HRTF, az=270., el=0.))

    print(measurements)
    for m in measurements:
        apply_IR(
                wavefile=wavefile,
                HRTF=HRTF,
                output="./sample{0}".format(m),
                measurement=m,)
