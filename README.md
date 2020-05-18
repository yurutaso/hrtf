## 音声ファイルにSOFA形式のHRTFを適用

- convolutionはnp.fftを使ってもよいが、データの整形も必要で面倒。scipy.signal.fftconvolveを使うと楽。
- モノラル音源に左右のHRTF(R_right, R_left)を適用。合成したステレオのwavファイルを作成する。
- sofaファイルについては、[sofaconventions.org](http://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics))を参照。
- HRTFのデータは、[ARI](https://www.kfs.oeaw.ac.at/index.php?view=article&id=608&lang=en)などを使用
- convolutionの実装については、[Qiitaの記事](https://qiita.com/stringamp/items/078e4f962e2073119b01)が参考になる。(この記事の実装内容は、scipyを使えばほとんど不要ではあるが。)

