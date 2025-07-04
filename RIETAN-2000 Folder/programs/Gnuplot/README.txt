Microsoft Windows 95/98/NT/2000 版 gnuplot-3.7.1+1.2.0 インストールマニュアル

               Copyright (C) 1993 - 2001   Masahito Yamaga


1. 含まれているファイル

        README.txt          このファイル
        Copyright.txt       オリジナルの gnuplot の著作権に関する文書
        Copyright_plus.txt  gnuplot-3.7.1+1.2 の著作権に関する文書
        plus.pdf            gnuplot-3.7.1+1.2 による機能拡張部分のマニュアル
        wgnuplot.exe        実行バイナリ
        wgnuplot.hlp        ヘルプファイル
        wgnuplot.mnu        メニューファイル (英語)
        wgnuplot.jmn        メニューファイル (日本語)
        VFlib.dll           VFlib (pbm, png, gif などの出力で
                            日本語表示を可能にするためのライブラリ)
        vfontcap_95         VFlib の 設定ファイル (Win95/98 用)
        vfontcap_NT         VFlib の 設定ファイル (WinNT/2000 用)
        vfontcap_95.pc98    VFlib の 設定ファイル (NEC PC98 + Win95/98 用)
        vfontcap_NT.pc98    VFlib の 設定ファイル (NEC PC98 + WinNT 用)
        zlib.dll            png デバイスに必要なライブラリ


2. コンパイル環境

Miscrosoft Visual C++ 6.0 (Windows 2000 サービスパック 1)

デフォルトのデバイスに加えて gd (gif) と png デバイスを追加してあります.

     gd (gif) = gd1.4 + gd1.4.kanji.patch + VFlib 2.24.1 + freetype-1.2
     png      = libpng-1.0.3 + zlib-1.1.3


3. インストール方法

 (1) アーカイブされている全てのファイルを同じディレクトリに置いてください.

 (2) Windows 95/98 か NT/2000 に応じてそれぞれ以下のように実行してください.

     (a) Windows 95/98 の場合
           vfontcap_95 を C:\WINDOWS\vfontcap として置いてください.

     (b) Windows NT/2000 の場合
           vfontcap_NT を C:\WINNT\vfontcap として置いてください.

     (c) NEC PC98 + Windows 95/98 の場合
           vfontcap_95.pc98 を A:\WINDOWS\vfontcap として置いてください.

     (d) NEC PC98 + Windows NT の場合
           vfontcap_NT.pc98 を A:\WINNT\vfontcap として置いてください.

 (3) メニューを日本語にしたい場合は wgnuplot.jmn を wgnuplot.mnu に名前を
     変更してください.

 (4) 画面上に日本語を表示したい場合は以下の手順を踏んでください.

       (a) gnuplot を起動して現れたウィンドウで右ボタンをクリック
       (b) メニューから Choose Font を選択
       (c) そこで日本語を扱えるフォントを選び OK をクリック
       (d) 再び右ボタンをクリックし, メニューから Update wgnuplot.ini を選択

     グラフを表示するウィンドウで日本語表示する場合も同様にします.
     最後に Update wgnuplot.ini を実行することを忘れないようにしてください.


4. 注意事項

著作権に関する文書 Copyright.txt と Copyright_plus.txt には必ず目を
通しておいてください.


5. Web サイト

最新情報は <http://www.yama-ga.com/gnuplot/> から取得できます.


                  平成 13年 1月 12日 山賀正人/Masahito Yamaga (ma@yama-ga.com)
