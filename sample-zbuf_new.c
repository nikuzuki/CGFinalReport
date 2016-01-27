//
// sample-zbuf.c
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "imgio.h"
#include "imgioX11.h"
#include "VecMat.h"

#define RESX 512 /* 画像サイズ(横) */
#define RESY 512 /* 画像サイズ(縦) */

#define FL   256 /* 透視投影変換の焦点距離 */

double zbuf[RESX][RESY]; // Zバッファ用配列

// 座標の変換
int ConvIJtoXY(int i, int j, double *x, double *y)
{
    *x = (double)(i - RESX/2);
    *y = (double)(RESY/2 - j);
    return 0;
}

// 座標の変換
int ConvXYtoIJ(double x, double y, int *i, int *j)
{
    *i = (int)(x + RESX/2);
    *j = (int)(RESY/2 - y);
    return 0;
}

// Zバッファの初期化
int InitZbuf()
{
    int i, j;

    for (j=0; j<RESY; j++) {
        for (i=0; i<RESX; i++) {
            zbuf[i][j] = -99999.99;
        }
    }
    return 0;
}

// 画像全体(IMG *img_out)を 指定したRGB値で塗りつぶす
int InitImage(IMG *img_out, int R, int G, int B)
{
    int i, j;
    for (j=0; j<img_out->height; j++) {
        for (i=0; i<img_out->width; i++) {
            R(*img_out, i, j) = R;
            G(*img_out, i, j) = G;
            B(*img_out, i, j) = B;
        }
    }
    return 0;
}

// 指定位置の画素に、指定色で点を打つ
int DrawPoint(IMG *img_out, int X, int Y, int r, int g, int b)
{
    int x, y;

    // 画像中央にスクリーンの原点を置き、右にx軸、左にy軸を取るように変換
    x = X + RESX/2; // 塗りたい画素の x座標値の計算
    y = RESY/2 -Y;  // 塗りたい画素の y座標値の計算（上下反転に注意）

    if (x>=0 && x<W(*img_out) && y>=0 && y<H(*img_out)) {
        R(*img_out, x, y) = r;
        G(*img_out, x, y) = g;
        B(*img_out, x, y) = b;
    }
    return 0;
}

int DrawLine(IMG 			*img_out, int x0, int y0, int x1, int y1, int r, int g, int b)
{
    int x, y, tmp;

    if (x0 >= x1) {
        tmp = x1; x1 = x0; x0 = tmp;
        tmp = y1; y1 = y0; y0 = tmp;
    }
    if (abs(y1 - y0)/(x1-x0+0.0001)<1) {
        for (x=x0; x<=x1; x++) {
            y = (int)((float)(y1 - y0)/(float)(x1 - x0)*(x - x0) + y0);
            DrawPoint(img_out, x, y, r, g, b);
        }
    } else {
        if (y0 >= y1) {
            tmp = x1; x1 = x0; x0 = tmp;
            tmp = y1; y1 = y0; y0 = tmp;
        }
        for (y=y0; y<=y1; y++) {
            x = (int)((float)(x1 - x0)/(float)(y1 - y0)*(y - y0) + x0);
            DrawPoint(img_out, x, y, r, g, b);
        }
    }

    return 0;
}



// 三角形の表示
int DispOneTriangle(IMG *img, Mat mat, Vec p0, Vec p1, Vec p2,
                    int R, int G, int B)
{
    int i, j;
    Vec q0, q1, q2, q10, q20, qq;
    Vec v0, v1, v2;
    Vec z01, z12, z20;
    Vec N;
    double t;

    //平行光線の定義
    Vec L;
    L[0] =  0;
    L[1] =  100;
    L[2] =  100;
    L[3] =  1;

    //視点ベクトルの定義
    Vec W;
    W[0] = 1;
    W[1] = 1;
    W[2] = 1;
    W[3] = 1;


    /*-----拡散反射の追加部分---------------------------------------*/
    int Rd, Gd, Bd;
    double kd = 1;
    double Ii = 10.0;
    double cos_alpha;



     /*-----鏡面反射の追加部分--------------------------------------*/

    Vec RR;
     RR[0] = 1;   RR[1] = 1;   RR[2] = 1;   RR[3] = 1;


    Vec Rn;
    Rn[0] = 1;  Rn[1] = 1;  Rn[2] = 1;  Rn[3] = 1;
    double Rn_norm;

    Vec Wn;
    Wn[0] = 1;  Wn[1] = 1;  Wn[2] = 1;  Wn[3] = 1;
    double Wn_norm;


    Vec Ld;
    Ld[0] = 1;  Ld[1] = 1;  Ld[2] = 1;  Ld[3] = 1;
    double Ld_norm;

    Vec Ln;
    Ln[0] = 1;  Ln[1] = 1;  Ln[2] = 1;  Ln[3] = 1;
    double Ln_norm;

    Vec Nn;
    Nn[0] = 1;  Nn[1] = 1;  Nn[2] = 1;  Nn[3] = 1;
    double Nn_norm;

    double p=50, cos_n;
    double ks = 10;

    /*-------------------------------------------------------------*/

    // 同次座標表現のための処理（一応）
    p0[3] = p1[3] = p2[3] = 1;

    // ３頂点を幾何変換
    MultipleVec(mat, p0, q0);
    MultipleVec(mat, p1, q1);
    MultipleVec(mat, p2, q2);

    // 法線ベクトルの計算（正規化なし）
    SubVec(q1, q0, q10);
    SubVec(q2, q0, q20);
    OuterProduct(q10, q20, N);

    /*拡散反射の係数計算*/
    cos_alpha = InnerProduct(N, L) / ( sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2])
					     * sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) );
    // printf("%lf \n",cos_alpha);

    /*鏡面反射の係数計算*/
    Ln_norm = sqrt( L[0]*L[0] + L[1]*L[1] + L[2]*L[2] );
    Ln[0] = L[0] / Ln_norm;
    Ln[1] = L[1] / Ln_norm;
    Ln[2] = L[2] / Ln_norm;

    Nn_norm = sqrt( N[0]*N[0] + N[1]*N[1] + N[2]*N[2] );
    Nn[0] = N[0] / Nn_norm;
    Nn[1] = N[1] / Nn_norm;
    Nn[2] = N[2] / Nn_norm;

    Ld_norm = fabs( InnerProduct(Ln, Nn) );
    Ld[0] = (-1) * Ln[0] / Ld_norm;
    Ld[1] = (-1) * Ln[1] / Ld_norm;
    Ld[2] = (-1) * Ln[2] / Ld_norm;

    //光の正反射方向ベクトル
    RR[0] = Ld[0] + 2 * Nn[0];
    RR[1] = Ld[1] + 2 * Nn[1];
    RR[2] = Ld[2] + 2 * Nn[2];


    Rn_norm = sqrt( RR[0]*RR[0] + RR[1]*RR[1] + RR[2]*RR[2] );
    Rn[0] = RR[0] / Rn_norm;
    Rn[1] = RR[1] / Rn_norm;
    Rn[2] = RR[2] / Rn_norm;

    Wn_norm = sqrt( W[0]*W[0] + W[1]*W[1] + W[2]*W[2] );
    Wn[0] = W[0] / Wn_norm;
    Wn[1] = W[1] / Wn_norm;
    Wn[2] = W[2] / Wn_norm;

    cos_n = InnerProduct(Rn, Wn);
    printf("%lf \n",cos_n);

    // 画面に対応する領域全体を走査して、三角形の内部か判定
    for (j=0; j< RESY; j++) {
        for (i=0; i<RESX; i++) {

            // レイの生成
            ConvIJtoXY(i, j, &qq[0], &qq[1]);
            qq[0] += 0.5; qq[1] += 0.5; qq[2] = -FL;

            // レイと三角形の乗る平面との交点計算
            t = InnerProduct(q0, N)/InnerProduct(N, qq);
            qq[0] *= t; qq[1] *= t; qq[2] *= t;

            // 交点から３頂点へのベクトルの外積のz成分を利用して内外判定
            SubVec(q0, qq, v0);
            SubVec(q1, qq, v1);
            SubVec(q2, qq, v2);
            OuterProduct(v0, v1, z01);
            OuterProduct(v1, v2, z12);
            OuterProduct(v2, v0, z20);

            // 三角形内部なら
            if ( (z01[2] >0 && z12[2] > 0 && z20[2] > 0) ||
                 (z01[2] <0 && z12[2] < 0 && z20[2] < 0)    ){

                // zbuf との比較
                if (t>0 && qq[2] > zbuf[i][j]) {
                    // 投影点を描画
			Rd = 0;//R * Ii * ( kd * cos_alpha + ks * pow(cos_n,p));
			Gd = (int)(G * Ii * ( kd * cos_alpha + ks * pow(cos_n,p)));
			Bd = 0;//B * Ii * ( kd * cos_alpha + ks * pow(cos_n,p));
                    DrawPoint(img, -FL*qq[0]/qq[2], -FL*qq[1]/qq[2], Rd,Gd,Bd);
                    // zbuf を更新
                    zbuf[i][j] = qq[2];
                }

            }
        }
    }
    return 0;
}

/* === main === */
int main(int argc, char *argv[])
{
    IMG img_out;
    FILE *fp;
    Mat mat; // 変換
    double th;
    Vec p[10][3];
    int c[10][3];
    int i;
    Vec p0, p1, p2, p3;

    fprintf(stderr, "START\n");
    if (argc != 2) {
        fprintf(stderr, "USAGE: %s output.ppm\n", argv[0]);
        exit(0);
    }

    imgioX11_InitWindow(); // 途中経過表示用の Window を生成
    cIMG(RESX, RESY, &img_out, COLOR); // IMG型のデータ生成
/*
	i=0;
        p[i][0][0] = -2; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] =  2; p[i][2][1] =  0; p[i][2][2] =  2;
	i=1;
        p[i][0][0] =  2; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] =  2; p[i][2][1] =  0; p[i][2][2] = -2;
	i=2;
        p[i][0][0] =  2; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] = -2; p[i][2][1] =  0; p[i][2][2] = -2;
	i=3;
        p[i][0][0] = -2; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] = -2; p[i][2][1] =  0; p[i][2][2] =  2;

	i=4;
        p[i][0][0] = -2; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] =  2; p[i][2][1] =  0; p[i][2][2] =  2;
	i=5;
        p[i][0][0] =  2; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] =  2; p[i][2][1] =  0; p[i][2][2] = -2;
	i=6;
        p[i][0][0] =  2; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] = -2; p[i][2][1] =  0; p[i][2][2] = -2;
	i=7;
        p[i][0][0] = -2; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] = -2; p[i][2][1] =  0; p[i][2][2] =  2;

*/

  i=0;
        p[i][0][0] = -2+5; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0+5; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] =  2+5; p[i][2][1] =  0; p[i][2][2] =  2;
	i=1;
        p[i][0][0] =  2+5; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0+5; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] =  2+5; p[i][2][1] =  0; p[i][2][2] = -2;
	i=2;
        p[i][0][0] =  2+5; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0+5; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] = -2+5; p[i][2][1] =  0; p[i][2][2] = -2;
	i=3;
        p[i][0][0] = -2+5; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0+5; p[i][1][1] = -5; p[i][1][2] =  0;
        p[i][2][0] = -2+5; p[i][2][1] =  0; p[i][2][2] =  2;

	i=4;
        p[i][0][0] = -2+5; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0+5; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] =  2+5; p[i][2][1] =  0; p[i][2][2] =  2;
	i=5;
        p[i][0][0] =  2+5; p[i][0][1] =  0; p[i][0][2] =  2;
        p[i][1][0] =  0+5; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] =  2+5; p[i][2][1] =  0; p[i][2][2] = -2;
	i=6;
        p[i][0][0] =  2+5; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0+5; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] = -2+5; p[i][2][1] =  0; p[i][2][2] = -2;
	i=7;
        p[i][0][0] = -2+5; p[i][0][1] =  0; p[i][0][2] = -2;
        p[i][1][0] =  0+5; p[i][1][1] =  5; p[i][1][2] =  0;
        p[i][2][0] = -2+5; p[i][2][1] =  0; p[i][2][2] =  2;


  for (i=0; i<8; i++){
    c[i][0] = 0;
    c[i][1] = 256;
    c[i][2] = 0;
  }

   //for (th=0; th<=360; th += 10.0){

  InitImage(&img_out, 30, 30, 255); // 画像内のクリア
  InitZbuf();

  InitMat(mat);                  // 変換行列の初期化(単位行列)
  ScaleMat(50, 50, 50, mat);     // 拡大
      //  RotateYMat( 0, mat);           // Y軸まわりの回転
     //   RotateXMat( 10, mat);           // Z軸まわりの回転
  TranslateMat(0, 0, -512, mat); // 立方体を視野の前方へ

  for (i=0; i<8; i++) {
    DispOneTriangle(&img_out, mat, p[i][0], p[i][1], p[i][2], c[i][0], c[i][1], c[i][2]);
    imgioX11_DisplayImage(img_out); /* 途中経過表示の更新 */
  }

 //  }

  wIMG(img_out, argv[1]);  /* 画像の書き出し */
  fprintf(stderr, "SAVED, and WAIT PUSH Q,q \n");
  imgioX11_RedrawLoop(img_out); /* イベント再描画設定。Q, q で終了 */
  return 0;
}
