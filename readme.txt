readme.txt 2021/7/20 杉浦圭祐
杉浦SPHコードの大まかな使い方などの説明です。
main.cpp param.h 及び my-src 以下の全てのソースコードは杉浦が作成したものですので、
改変などを自由にしていただいて大丈夫です。
src 以下のソースコードはFDPSのライブラリ用のコードなので、書き換えない方がいいと思います。
用いているモデルの論文などはdocuments以下に置いてありますので、そちらをご覧ください。


1. コンパイル
Makefileのあるディレクトリで make をすることでコンパイルできます。
コンパイルには main.cpp, param.h と my-src および src のディレクトリ以下のソースコードが必要なので、
これらが一緒にコンパイルされるように PS_PATH を書き換えてください。
コンパイラは CC = の右辺を書き換えることで適宜書き換えてください。
OpenMP もしくは MPI を用いた並列計算を行いたいときは、それぞれ
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL (コンパイラによってはここに -fopenmp オプションが必要)
もしくは
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
の行のコメントアウトを外してからコンパイルをお願いします。


2. 実行
実行時には、第一変数に実行ファイル名を、第二変数に初期条件ファイル名を、第三変数に計算結果の連番ファイルが格納されるディレクトリ名を指定します。
デフォルトの Makefile でコンパイルすると sph.out という実行ファイルが出来上がるので、
例えば init.bin という初期条件ファイルから計算を始めたい場合は (初期条件ファイルの作り方については make-initial-conditionsディレクトリ内の書類などをご覧ください)、
./sph.out init.bin data-directory
として実行してください。
MPI並列計算をする場合は、例えば4MPI並列の計算を実行する場合には、
mpirun -np 4 ./sph.out init.bin data-directory
などとして実行してください (mpirunのコマンドはインストールしたものを使用してください)。
天文台XC50で実行するやり方は、 「XC50への計算の投げ方.txt」をご覧ください。


3. 主なソースコードの簡単な説明
3.1. main.cpp
時間発展とデータの書き出しをする部分のファイルです。

14行目の
char filename[]="FDPS-basalt-sphere-collision-0000.txt";
で、全SPH粒子の位置などの情報を記した連番のASCIIファイルの名前を指定しています。(このファイルの書式の指定方法は後述)

20行目の
double t_end=2000.0;
で計算終了する数値計算内の時刻を指定します。

23行目の
double t_interval=10.0;
で連番ファイルを出力する時間間隔を指定します。

35行目からの部分で連番ファイルなどを格納するディレクトリを作成していますが、
すでにそのディレクトリが存在する場合は作成できなくて計算を終了させてしまうので、ご注意ください。

55行目および73行目からの
sprintf(filename,"FDPS-basalt-sphere-collision-%04d.txt",l);
sprintf(dataname,"%s/%s",datadir,filename);
header.time = getTime();
header.Nbody = sph_system.getNumberOfParticleGlobal();
sph_system.writeParticleAscii(dataname,header);
で連番ファイルの書き出しを行なっています。
書き出すデータの種類や書式を変えるには、my-src/class.h の13行目からの
void writeAscii(FILE* fp) const{
    fprintf(fp, "%e\n", time);
    fprintf(fp, "%d\n", Nbody);
} #ヘッダ部分
および178行目からの
void writeAscii(FILE* fp) const{
    double km=1.0e3;
    fprintf(fp,"%.2f %.2f %.2f %.2f %lld %d %.4e %.4e %.4e %.4e \n",this->pos.x/km,this->pos.y/km,this->pos.z/km,this->smth/km,this->id,this->property_tag,this->dens,this->eng,this->pres,this->temp);
} #SPH粒子の情報の部分
を書き換えてください。
各SPH粒子に割り当てられている質量などの情報とその変数名は、my-src/class.h の83行目からの場所で確認することができます。

63行目のfor文がメインループです。
数値計算内の時刻に達したらfor文を抜ける他、
天文台XC50のように現実での経過時間で強制的に計算が終了させられてしまう前に計算を止めて中間ファイルを出力するために、
(sim_now-sim_start)<60*60*11.75
で現実での時間がどれくらい経過したら終了するか指定します。

65行目の
one_time_development_roop(sph_system,dinfo,dens_tree,gradients_tree,elast_tree,grav_tree);
の関数で1タイムステップ時間発展します。

89行目からの
header.time = getTime();
header.Nbody = sph_system.getNumberOfParticleGlobal();
sph_system.writeParticleBinary(midfile,header);
で、midfile.bin という名前の中間ファイルを出力します。
この中間ファイルは計算を継続するための全ての情報を記したファイルで、
このファイルを初期条件ファイルとして計算を投げ直すと、計算を継続することができます。

最後に、93行目からの箇所で、
数値計算内の時刻がt_endに達していた場合は、最終状態の全ての情報を記録した.binファイルを、連番ファイルと同じディレクトリに書き出します。
この.binファイルは最終状態の解析などに使用します。

3.2. param.h
パラメータや使用するモデル、手法を指定するファイルです。
弾性体か流体か、ひび割れモデルを使用するか否か、などの手法の変更は基本的にはこのファイルを書き換えるだけでできます。
それぞれの手法などに必要なパラメータもここで指定してください。
2種類以上の物質を使用したい場合は、例えば2種類使用する場合は36行目で、
const PS::S32 N_MATERIAL = 2; //the number of type of materials
とした上で、例えば44行目のせん断弾性率のように{}を使用して配列形式で書いてあるパラメータについては、
const PS::F64 MU_SHEAR[N_MATERIAL] = {2.0e8, 2.0e10};//shear modulus
というように使用したい物質の数だけ指定してください。
sph_system[i].property_tag = 0 となっているi番目のSPH粒子は上記配列の0番目のパラメータを、
sph_system[i].property_tag = 1 となっているi番目のSPH粒子は上記配列の1番目のパラメータを使用することになります。
粉体の計算として最も大切なパラメータである粉体の摩擦係数は, 50行目の MU_D_FRIC です.
49行目の MU_I_FRIC は一枚岩のパラメータなので少し違います.

3.3. my-src/class.h
4.1でも説明したように、それぞれのSPH粒子に付随する変数の定義と書き出し、書き込みに必要な関数の定義を行っているファイルです。
SPH粒子に付随する各変数の意味については、83行目からの部分に書いてあります。

3.4. my-src/force.h
全ての相互作用関数 (運動方程式、連続の式、エネルギー方程式など) が定義されているファイルです。
