# FDPSを利用した弾性体・岩石破壊・粉体摩擦・流体などを扱えるSPHコード (杉浦コード)

以下に記載のFDPSフレームワークを利用してMPI・OpenMP並列化を施したコードである. このコードを使用して論文を書いた時は, FDPSに関する以下の2つの論文を必ず引用すること.

- Iwasawa, M., Tanikawa, A., Hosono, N., et al. 2015, in Proceedings of the 5th International Workshop on Domain-Specific Languages and High-Level Frameworks for High Performance Computing, WOLFHPC ’15 (New York, NY: ACM), 1:1

- Iwasawa, M., Tanikawa, A., Hosono, N., et al. 2016, PASJ, 68, 541

またこのコードは基本的に以下の論文に基づいているため, この論文も必ず引用すること.

- Sugiura, K., Kobayashi, H., Inutsuka, S., 2018, AA, 620, A168

## 国立天文台アテルイIIでのコンパイル時の注意

FDPSバージョン7ではC++17の機能を明示的に使っているが, アテルイIIのデフォルトコンパイラのcrayではC++17を指定する方法が分からなかった. そのためアテルイIIにsshログインした後にコンパイルする時は以下のコマンドで環境をgnuに変えてから行うこと.

$ module switch PrgEnv-cray PrgEnv-gnu

# FDPS

FDPS is a general-purpose, high-performance library for particle simulations.

The current version is 7.1. The previous versions are [here](https://github.com/FDPS/FDPS/releases).

We maintain this from subversion-over-github interface.

If you have some questions, please do not hesitate to contact us. Our
e-mail address is fdps-support-@-mail.jmlab.jp (please replace -@- by @).

Tutorial of FDPS is here
[doc/doc_tutorial_cpp_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_cpp_en.pdf?raw=true)
, Specification is here
[doc/doc_specs_cpp_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_cpp_en.pdf?raw=true).
Tutorial and specification of Fortran/C interface to FDPS are
[doc/doc_tutorial_ftn_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_ftn_en.pdf?raw=true)
and
[doc/doc_specs_ftn_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_ftn_en.pdf?raw=true),
respectively.
You can also find a two-page handout here
[doc/doc_SC15_handout.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_SC15_handout.pdf?raw=true).


FDPSのチュートリアルは
[doc/doc_tutorial_cpp_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_cpp_ja.pdf?raw=true)
、仕様書は
[doc/doc_specs_cpp_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_cpp_ja.pdf?raw=true)
にあります。
また、FDPSのFortran/C言語インターフェースのチュートリアルは
[doc/doc_tutorial_ftn_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_ftn_ja.pdf?raw=true)
、Fortran/C言語インターフェースの仕様書は
[doc/doc_specs_ftn_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_ftn_ja.pdf?raw=true)
にあります。
2ページでFDPSがわかるハンドアウトはこちらです
[doc/doc_SC15_handout.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_SC15_handout.pdf?raw=true)。

ご質問などお問い合わせはfdps-support-@-mail.jmlab.jpにお願いいたします (-@-を@に置換して下さい)。
