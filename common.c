#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mrp.h"

const CPOINT dyx[] = {
	/* 1 */
	{ 0,-1}, {-1, 0},
	/* 2 */
	{ 0,-2}, {-1,-1}, {-2, 0}, {-1, 1},
	/* 3 */
	{ 0,-3}, {-1,-2}, {-2,-1}, {-3, 0}, {-2, 1}, {-1, 2},
	/* 4 */
	{ 0,-4}, {-1,-3}, {-2,-2}, {-3,-1}, {-4, 0}, {-3, 1}, {-2, 2}, {-1, 3},
	/* 5 */
	{ 0,-5}, {-1,-4}, {-2,-3}, {-3,-2}, {-4,-1}, {-5, 0}, {-4, 1}, {-3, 2},
	{-2, 3}, {-1, 4},
	/* 6 */
	{ 0,-6}, {-1,-5}, {-2,-4}, {-3,-3}, {-4,-2}, {-5,-1}, {-6, 0}, {-5, 1},
	{-4, 2}, {-3, 3}, {-2, 4}, {-1, 5},
	/* 7 */
	{ 0,-7}, {-1,-6}, {-2,-5}, {-3,-4}, {-4,-3}, {-5,-2}, {-6,-1}, {-7, 0},
	{-6, 1}, {-5, 2}, {-4, 3}, {-3, 4}, {-2, 5}, {-1, 6},
	/* 8 */
	{ 0,-8}, {-1,-7}, {-2,-6}, {-3,-5}, {-4,-4}, {-5,-3}, {-6,-2}, {-7,-1},
	{-8, 0}, {-7, 1}, {-6, 2}, {-5, 3}, {-4, 4}, {-3, 5}, {-2, 6}, {-1, 7},
	/* 9 */
	{ 0,-9}, {-1,-8}, {-2,-7}, {-3,-6}, {-4,-5}, {-5,-4}, {-6,-3}, {-7,-2},
	{-8,-1}, {-9, 0}, {-8, 1}, {-7, 2}, {-6, 3}, {-5, 4}, {-4, 5}, {-3, 6},
	{-2, 7}, {-1, 8},
	/* 10 */
	{ 0,-10}, {-1,-9}, {-2,-8}, {-3,-7}, {-4,-6}, {-5,-5}, {-6,-4}, {-7,-3},
	{-8,-2}, {-9,-1}, {-10, 0}, {-9, 1}, {-8, 2}, {-7, 3}, {-6, 4}, {-5, 5},
	{-4, 6}, {-3, 7}, {-2, 8}, {-1, 9},

};

double sigma_a[] = {0.15, 0.26, 0.38, 0.57, 0.83, 1.18, 1.65, 2.31,
3.22, 4.47, 6.19, 8.55, 11.80, 16.27, 22.42, 30.89};

double qtree_prob[] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};

double zerocoef_prob[NUM_ZMODEL] = {
	0.003,0.010,0.020,0.033,0.046,0.062,0.079,0.097,0.116,0.135,0.156,0.177,0.200,
	0.222,0.246,0.270,0.294,0.319,0.344,0.370,0.396,0.421,0.447,0.474,0.500,0.526,
	0.552,0.578,0.604,0.630,0.656,0.681,0.706,0.730,0.754,0.778,0.800,0.823,0.844,
	0.865,0.884,0.903,0.921,0.938,0.954,0.967,0.980,0.990,0.997,
};


FILE *fileopen(char *filename, char *mode)
{
	FILE *fp;
	fp = fopen(filename, mode);
	if (fp == NULL) {
		fprintf(stderr, "Can\'t open %s!\n", filename);
		exit(1);
	}
	return (fp);
}

void *alloc_mem(size_t size)
{
	void *ptr;
	if ((ptr = (void *)malloc(size)) == NULL) {
		fprintf(stderr, "Can\'t allocate memory (size = %d)!\n", (int)size);
		exit(1);
	}
	return (ptr);
}

void **alloc_2d_array(int height, int width, int size)
{
	void **mat;
	char *ptr;
	int k;

	mat = (void **)alloc_mem(sizeof(void *) * height + height * width * size);
	ptr = (char *)(mat + height);
	for (k = 0; k < height; k++) {
		mat[k] =  ptr;
		ptr += width * size;
	}
	return (mat);
}

/*
int ***alloc_3d_array( int n0, int n1, int n2 ) {
  int ***b;
  int i,j;
  b= (int***) malloc( sizeof(int**) * n0 *64);
  for ( i= 0; i < n0; i++ ) {
    b[i]= (int**) malloc( sizeof(int*) * n1  *64);//なぜ64を掛けないといけないのか?後日要検証
  }
  b[0][0]= (int*) malloc( sizeof(int) * n0 * n1 * n2 *64);
  for ( i= 0; i < n0; i++ ) {
    for(j=0;j<n1;j++){
      b[i][j]= (int*) ( b[0][0] + ( j + i*n1 )*n2 );
    }
  }
  return b;
}
*/
void ***alloc_3d_array(int height, int width, int depth, int size)
{
  void ***mat;
  char **ptr, *ptrp;
  int k,d;

  mat = alloc_mem(sizeof(void **) * height
      + sizeof(void *) * height * width
      + height * width * depth * size);
  ptr = (char **)(mat + height);
  ptrp = (char *)(ptr + height * width);

  for (k = 0; k < height; k++) {
    mat[k] = (void **)ptr;
    ptr += width;
    for (d = 0; d < width; d++) {
      mat[k][d] = ptrp;
      ptrp += depth * size;
    }
  }

  return (mat);
}
 



IMAGE *alloc_image(int width, int height, int maxval)
{
	IMAGE *img;
	img = (IMAGE *)alloc_mem(sizeof(IMAGE));
	img->width = width;
	img->height = height;
	img->maxval = maxval;
	img->val = (img_t **)alloc_2d_array(img->height, img->width,
		sizeof(img_t));
	return (img);
}

/*
Natural logarithm of the gamma function
cf. "Numerical Recipes in C", 6.1
http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf.html
*/
double lngamma(double xx)//Γ関数
{
	int j;
	double x,y,tmp,ser;
	double cof[6] = {
		76.18009172947146,	-86.50532032941677,
		24.01409824083091,	-1.231739572450155,
		0.1208650973866179e-2,	-0.5395239384953e-5
	};

	y = x = xx;
	tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++)
		ser += (cof[j] / ++y);
	return (log(2.5066282746310005 * ser / x) - tmp);
}

double calc_ggprob(double beta, double shape, double h, double x)
{
	double p;

	if (x < 0.0) x = -x;//絶対値
	if (x < 1E-6) {// 1 / 10^6.百万分の一
		p = exp(-pow(beta * h, shape)) + exp(0.0);//(beta * h)^shape + 1.x = 0に近似.式(2-4)上右側.xは誤差e
	} else {
		p = exp(-pow(beta * (x - h), shape))//確率モデルの面積を求める際,全ての高さを足すが1/8ずらしているので同じ値でない
			+ exp(-pow(beta * (x + h), shape));//よってhによって補完する
	}
	return (p);
}

void set_freqtable(PMODEL *pm, int size, int ssize,
				   double shape, double sigma, double h, double off)
{
	double beta, norm;
	int i;

	pm->size = size; //511
	pm->freq = (uint *)alloc_mem((2 * size + 1) * sizeof(uint));//一次元配列
	pm->cumfreq = &pm->freq[size];

	/* Generalized Gaussian distribution */
	beta = exp(0.5*(lngamma(3.0/shape)-lngamma(1.0/shape))) / sigma;//一般化ガウス関数.ηのみ?.論文式(2-4).shape = cn
	norm = 0.0;//初期化
	for (i = 0; i < size; i++) {
		norm += calc_ggprob(beta, shape, h, i - off);//式(2-4)上式のexpとその中身
	}
	norm = (double)(MAX_TOTFREQ - size * MIN_FREQ) / norm;//MAX_TOTFREQ = (1 << 14) . MIN_FREQ = 1.
	norm += 1E-8;	/* to avoid machine dependent rounding errors */ //量子化するとき誤差が出ないようにするため
	pm->norm = norm;
	pm->cumfreq[0] = 0;
	for (i = 0; i < size; i++) {
		pm->freq[i] = (uint)(norm * calc_ggprob(beta, shape, h, i - off) + MIN_FREQ);//誤差511個分の確率モデルの高さが入る
		pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];//確率モデルの高さの合計
	}
	if (ssize > 0) {
		pm->cost = (float *)alloc_mem((size + ssize) * sizeof(float));
		pm->subcost = &pm->cost[size];
	}
	return;
}

PMODEL ***init_pmodels(int num_group, int num_pmodel, int pm_accuracy,
					   int *pm_idx, double *sigma, int size)
{
	PMODEL ***pmodels, *pmbuf, *pm;
	int gr, i, j, num_subpm, ssize;
	double delta_c, c, s, sw, off;

	num_subpm = 1 << pm_accuracy;//num_subpm is 8.pm_accuracy is 3
	ssize = size;//size is 256
	size = size + ssize - 1;//size = 511
	sw = 1.0 / (double)num_subpm;// 1/8の量子化

	delta_c = 3.2 / (double)num_pmodel;//cnは0.2ずつ動くので. num_pmodel is 16 . delta_c = 0.2 
	off = (double)(ssize - 1);//off = 255
	if (pm_idx != NULL) {
		num_pmodel = 1;
		ssize = 0;
	}
	pmodels = (PMODEL ***)alloc_2d_array(num_group, num_pmodel, sizeof(PMODEL *));
	pmbuf = (PMODEL *)alloc_mem(num_group * num_pmodel * num_subpm * sizeof(PMODEL));
	for (gr = 0; gr < num_group; gr++) {//gr は分散の番号16種類
		s = sigma[gr];
		for (i = 0; i < num_pmodel; i++) {//確率モデルの番号16種類
			pmodels[gr][i] = pmbuf;//3次元配列にした
			for (j = 0; j < num_subpm; j++) {//1/8の場所どれにするか.8種類
				pm = pmbuf++;
				pm->id = i;
				if (num_pmodel > 1) {
					c = delta_c * (double)(i + 1);//常にこれを使用する　idとcnが配列の使用によってずれてしまうためidに1を足す
				} else if (pm_idx != NULL) {
					c = delta_c * (double)(pm_idx[gr] + 1);
				} else {
					c = 2.0;
				}
				if (c < 0.1) c = 0.1;
				set_freqtable(pm, size, ssize, c, s, sw/2.0, off - sw * j);
			}
		}
	}
	return (pmodels);
}

/* probaility model for coefficients and thresholds */
void set_spmodel(PMODEL *pm, int size, int m)
{
	int i, sum;
	double p;

	pm->size = size;
	if (m >= 0) {
		p = 1.0 / (double)(1 << (m % 8));
		sum = 0;
		for (i = 0; i < pm->size; i++) {
			pm->freq[i] = (uint)(exp(-p * i) * (1 << 10));
			if (pm->freq[i] == 0) pm->freq[i]++;
			sum += pm->freq[i];
		}
		if (m & 8) pm->freq[0] = (sum - pm->freq[0]);	/* weight for zero */
	} else {
		for (i = 0; i < pm->size; i++) {
			pm->freq[i] = 1;
		}
	}
	pm->cumfreq[0] = 0;
	for (i = 0; i < pm->size; i++) {
		pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
	}
	return;
}





double *init_ctx_weight_double(void)
{
	double *ctx_weight;
  int k;
	double dy, dx;

	ctx_weight = (double *)alloc_mem(NUM_UPELS * sizeof(double));
	for (k = 0; k < NUM_UPELS; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;
#if !CTX_WEIGHT
		ctx_weight[k] = 64;
#elif MHD_WEIGHT
		ctx_weight[k] = (64.0 / (fabs(dy) + fabs(dx)) + 0.5);
#else
		ctx_weight[k] = (64.0 / sqrt(dy * dy + dx * dx) + 0.5);
#endif
	}
	return (ctx_weight);
}




int *init_ctx_weight(void)
{
	int *ctx_weight, k;
	double dy, dx;

	ctx_weight = (int *)alloc_mem(NUM_UPELS * sizeof(int));
	for (k = 0; k < NUM_UPELS; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;
#if !CTX_WEIGHT
		ctx_weight[k] = 64;
#elif MHD_WEIGHT
		ctx_weight[k] = (int)(64.0 / (fabs(dy) + fabs(dx)) + 0.5);
#else
		ctx_weight[k] = (int)(64.0 / sqrt(dy * dy + dx * dx) + 0.5);
#endif
	}
	return (ctx_weight);
}

void mtf_classlabel(char **class, int *mtfbuf, int y, int x,
					int bsize, int width, int num_class)//move to front
{
	int i, j, k, ref[3];

	if (y == 0) {
		if (x == 0) {
			ref[0] = ref[1] = ref[2] = 0;//(0,0)の時は周りを0とする
		} else {
			ref[0] = ref[1] = ref[2] = class[y][x-1];//(0,x)の時はひとつ前の値(0,x-1)にする
		}
	} else {
		ref[0] = class[y-1][x];//上
		ref[1] = (x == 0)? class[y-1][x] : class[y][x-1];//xが0なら上，そうでないなら一個前
		ref[2] = (x + bsize >= width)?
		class[y-1][x] : class[y-1][x+bsize];//最終ブロックなら上，そうでないなら右上
		if (ref[1] == ref[2]) {
			ref[2] = ref[0];
			ref[0] = ref[1];
		}
	}
	/* move to front */
	for (k = 2; k >= 0; k--) {
		if ((j = mtfbuf[ref[k]]) == 0) continue;
		for (i = 0; i < num_class; i++) {
			if (mtfbuf[i] < j) {
				mtfbuf[i]++;
			}
		}
		mtfbuf[ref[k]] = 0;
	}
	return;
}

double cpu_time(void)
{
#ifndef HAVE_CLOCK
#  include <sys/times.h>
	struct tms t;
#endif
#ifndef CLK_TCK
#  define CLK_TCK 60
#endif
	static clock_t prev = 0;
	clock_t cur, dif;

#ifdef HAVE_CLOCK
	cur = clock();
#else
	times(&t);
	cur = t.tms_utime + t.tms_stime;
#endif
	if (cur > prev) {
		dif = cur - prev;
	} else {
		dif = (unsigned)cur - prev;
	}
	prev = cur;

#ifdef HAVE_CLOCK
	return ((double)dif / CLOCKS_PER_SEC);
#else
	return ((double)dif / CLK_TCK);
#endif
}
