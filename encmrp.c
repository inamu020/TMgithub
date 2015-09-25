/***** Encoder *****/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "mrp.h"

extern CPOINT dyx[];
extern double sigma_a[];
extern double qtree_prob[];
extern double zerocoef_prob[];

float ****calc_entropy_of_conditional_probability(PMODEL ***pmodels, int num_group,
												  int num_pmodel, int pm_accuracy, 
												  int maxval)
{
	int i, j, k, l, gr, total;
	uint *freq_p, *cumfreq_p;
	PMODEL *pm_p, *pm;
	double *logfreq, logtotal, log2 = log(2.0);
	float **ptr1, *ptr2, ****c_prob, entropy;

	/* alloc 4D memory */
	c_prob = (float ****)alloc_2d_array(num_group, num_pmodel, sizeof(float **));
	ptr1 = (float **)alloc_mem(num_group * num_pmodel * (1 << pm_accuracy) 
		* sizeof(float *));
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pmodel; i++) {
			c_prob[gr][i] = ptr1;
			ptr1 += (1 << pm_accuracy);
		}
	}
	ptr2 = (float *)alloc_mem(num_group * num_pmodel * (1 << pm_accuracy)
		* (maxval + 1) * sizeof(float));
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pmodel; i++) {
			for (j = 0; j < (1 << pm_accuracy); j++) {
				c_prob[gr][i][j] = ptr2;
				ptr2 += maxval + 1;
			}
		}
	}
	/* calc entropy of conditional probability */
	logfreq = (double *)alloc_mem(((maxval << 1) + 1) * sizeof(double));
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pmodel; i++) {
			pm_p = pmodels[gr][i];
			for (j = 0; j < (1 << pm_accuracy); j++) {
				pm = pm_p + j;
				freq_p = pm->freq;
				cumfreq_p = pm->cumfreq;
				for (k = 0; k < (maxval << 1) + 1; k++) {
					logfreq[k] = log(freq_p[k]);
				}
				for (k = 0; k < maxval + 1; k++) {
					total = cumfreq_p[k + maxval + 1] - cumfreq_p[k];
					logtotal = log(total);
					entropy = 0;
					for (l = k; l < k + maxval + 1; l++) {
						entropy += (float)(freq_p[l] * (logtotal - logfreq[l]));
					}
					entropy /= (float)(total * log2);
					c_prob[gr][i][j][k] = entropy;
				} 
			}
		}
	}
	return c_prob;
}

void calc_ratio_of_model_to_rate(ENCODER *enc)
{
	int y, x, prd, gr, e, base, frac, cl, totpel, *gr_pel;
	double *entropy, *cost, *ratio, totent, totcost, totratio;
	PMODEL *pm;
	static int calc_entropy_flag = 0;
	static float ****entropy_of_conditional_probability;


	if (calc_entropy_flag == 0) {
		entropy_of_conditional_probability =
			calc_entropy_of_conditional_probability(enc->pmodels, enc->num_group,
			enc->num_pmodel,enc->pm_accuracy,
			enc->maxval);
		calc_entropy_flag = 1;
	}
	entropy = (double *)alloc_mem(enc->num_group * sizeof(double));
	cost = (double *)alloc_mem(enc->num_group * sizeof(double));
	ratio = (double *)alloc_mem(enc->num_group * sizeof(double));
	gr_pel = (int *)alloc_mem(enc->num_group * sizeof(int));
	for (gr = 0; gr < enc->num_group; gr++) {
		entropy[gr] = 0;
		cost[gr] = 0;
		gr_pel[gr] = 0;
	}
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			cl = enc->class[y][x];
			gr = enc->group[y][x];
			e = enc->encval[y][x];
			prd = enc->prd[y][x];
			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost[gr] += pm->cost[base + e] + pm->subcost[base];
			entropy[gr] += entropy_of_conditional_probability[gr][pm->id][frac][base];
			gr_pel[gr]++;
		}
	}
	/* calc ratio */
	totpel = 0;
	totcost = totent = 0;
	for (gr = 0; gr < enc->num_group; gr++){
		if (entropy[gr] != 0.0)
			ratio[gr] = (1.0 - fabs(entropy[gr] - cost[gr]) / entropy[gr]);
		else
			ratio[gr] = 1.0;
		totent += entropy[gr];
		totcost += cost[gr];
		totpel += gr_pel[gr];
	}
	totratio = (1.0 - fabs(totent - totcost) / totent);
	/* print */
	printf("******* differences between entropy and rate *******\n");
	printf("(gr)  [shape]\tentropy\t\t|rate\t\t|fitness\t|pel\n");
	printf("------------------------------------------------------------------------------\n");
	for (gr = 0; gr < enc->num_group; gr++) {
		printf("(%2d)  [%.1f]\t%10d\t|%10d\t|   %.3f\t|%10d\n", gr,
			0.2 * (enc->pmlist[gr]->id + 1), (int)entropy[gr],
			(int)cost[gr], ratio[gr], gr_pel[gr]);
	}
	printf("------------------------------------------------------------------------------\n");
	printf("all         \t%10d\t|%10d\t|   %.3f\t|%10d\n", (int)totent, (int)totcost, totratio, totpel);
	free(entropy);
	free(cost);
	free(ratio);
	free(gr_pel);
}

IMAGE *read_pgm(char *filename)
{
	int i, j, width, height, maxval;
	char tmp[256];
	IMAGE *img;
	FILE *fp;

	fp = fileopen(filename, "rb");
	fgets(tmp, 256, fp);
	if (tmp[0] != 'P' || tmp[1] != '5') {
		fprintf(stderr, "Not a PGM file!\n");
		exit(1);
	}
	while (*(fgets(tmp, 256, fp)) == '#');
	sscanf(tmp, "%d %d", &width, &height);
	while (*(fgets(tmp, 256, fp)) == '#');
	sscanf(tmp, "%d", &maxval);

	if (maxval > 255) {
		fprintf(stderr, "Sorry, this version only supports 8bpp images!\n");
		exit(1);
	}
	img = alloc_image(width, height, maxval);
	for (i = 0; i < img->height; i++) {
		for (j = 0; j < img->width; j++) {
			img->val[i][j] = (img_t)fgetc(fp);
		}
	}
	fclose(fp);
	return (img);
}

int ***init_ref_offset(int height, int width, int prd_order)//予測器の範囲指定かと思う
{
	int ***roff, *ptr;
	int x, y, dx, dy, k;
	int order, min_dx, max_dx, min_dy;

	min_dx = max_dx = min_dy = 0;
	order = (prd_order > NUM_UPELS)? prd_order : NUM_UPELS;
	for (k = 0; k < order; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;
		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dx > max_dx) max_dx = dx;
	}
	roff = (int ***)alloc_2d_array(height, width, sizeof(int *));//3次元配列作成の例
	if (min_dy != 0) {
		ptr = (int *)alloc_mem((1 - min_dy) * (1 + max_dx - min_dx) * order * sizeof(int));
	}else {
		ptr = (int *)alloc_mem((2 + max_dx - min_dx) * order * sizeof(int));
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			if (y == 0) {
				if (x == 0) {
					roff[y][x] = ptr;
					dx = 0;
					dy = height;
					for (k = 0; k < order; k++) {
						*ptr++ = dy * width + dx;
					}
				} else if (x + min_dx <= 0 || x + max_dx >= width) {
					roff[y][x] = ptr;
					dy = 0;
					for (k = 0; k < order; k++) {
						dx = dyx[k].x;
						if (x + dx < 0) dx = -x;
						else if (dx >= 0) dx = -1;
						*ptr++ = dy * width + dx;
					}
				} else {
					roff[y][x] = roff[y][x - 1];
				}
				// for K = 1 and NUM_UPELS = 1
			} else if (min_dy == 0 && y == 1 && x == 0) {
				roff[y][x] = ptr;
				dy = -1;
				dx = 0;
				*ptr++ = dy * width + dx;
			} else if (y + min_dy <= 0) {
				if (x == 0) {
					roff[y][x] = ptr;
					for (k = 0; k < order; k++) {
						dy = dyx[k].y;
						if (y + dy < 0) dy = -y;
						else if (dy >= 0) dy = -1;
						dx = dyx[k].x;
						if (x + dx < 0) dx = -x;
						*ptr++ = dy * width + dx;
					}
				} else if (x + min_dx <= 0 || x + max_dx >= width) {
					roff[y][x] = ptr;
					for (k = 0; k < order; k++) {
						dy = dyx[k].y;
						if (y + dy < 0) dy = -y;
						dx = dyx[k].x;
						if (x + dx < 0) dx = -x;
						else if (x + dx >= width) {
							dx = width - x - 1;
						}
						*ptr++ = dy * width + dx;
					}
				} else {
					roff[y][x] = roff[y][x - 1];
				}
			} else {
				roff[y][x] = roff[y - 1][x];
			}
		}
	}
	return (roff);
}

ENCODER *init_encoder(IMAGE *img, int num_class, int num_group,
					  int prd_order, int coef_precision, int quadtree_depth,
					  int num_pmodel, int pm_accuracy)
{
	ENCODER *enc;
	int x, y, i, j;
	double c;
#if AUTO_PRD_ORDER
	cost_t *ptr;
	int k;
#endif

	enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
	enc->height = img->height;
	enc->width = img->width;
	enc->maxval = img->maxval;
	enc->num_class = num_class;
	enc->num_group = num_group;
#if AUTO_PRD_ORDER
	enc->prd_order = BASE_PRD_ORDER;
	enc->max_prd_order = MAX_PRD_ORDER;
#else
	enc->prd_order = prd_order;
	enc->max_prd_order = prd_order;
#endif
	enc->coef_precision = coef_precision;//6
	enc->max_coef = (2 << coef_precision);
	enc->quadtree_depth = quadtree_depth;
	enc->num_pmodel = num_pmodel;//16
	enc->pm_accuracy = pm_accuracy;//3
	enc->maxprd = enc->maxval << enc->coef_precision;//maxvalを64培したもの  
	enc->predictor = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order,
		sizeof(int));
	enc->num_nzcoef = (int *)alloc_mem(enc->num_class * sizeof(int));
#if AUTO_PRD_ORDER
	for (i = 0; i < enc->num_class; i++) {
		enc->num_nzcoef[i] = BASE_PRD_ORDER;
	}
#else
	for (i = 0; i < enc->num_class; i++) {
		enc->num_nzcoef[i] = prd_order;
	}
#endif
	enc->nzconv = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order, sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->nzconv[i][j] = j;
		}
	}
#if AUTO_PRD_ORDER
	enc->num_search = (int *)alloc_mem(enc->num_class * sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		enc->num_search[i] = 30;
	}
	enc->ord2mhd = (int *)alloc_mem(enc->max_prd_order * sizeof(int));
	for (i = 0; i < enc->max_prd_order; i++) {
		enc->ord2mhd[i] = (int)((sqrt(1 + 4 * i) - 1) / 2);
	}
	enc->prd_mhd = enc->ord2mhd[enc->max_prd_order - 1] + 1;
#endif
	enc->th = (int **)alloc_2d_array(enc->num_class, enc->num_group,
		sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->predictor[i][j] = 0;
		}
		for (j = 0; j < enc->num_group - 1; j++) {
			enc->th[i][j] = 0;
		}
		enc->th[i][enc->num_group - 1] = MAX_UPARA + 1;
	}

	enc->upara = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	enc->prd = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	enc->prd_b = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	enc->roff = init_ref_offset(enc->height, enc->width, enc->max_prd_order);
	enc->org = (int **)alloc_2d_array(enc->height+1, enc->width, sizeof(int));
	enc->err = (int **)alloc_2d_array(enc->height+1, enc->width, sizeof(int));
	if (enc->quadtree_depth > 0) {
		y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
		x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = enc->quadtree_depth - 1; i >= 0; i--) {
			enc->qtmap[i] = (char **)alloc_2d_array(y, x, sizeof(char));
			y <<= 1;
			x <<= 1;
		}
	}
	enc->ctx_weight = init_ctx_weight();
  enc->ctx_weight_double = init_ctx_weight_double();
	enc->class = (char **)alloc_2d_array(enc->height, enc->width,
		sizeof(char));
	enc->group = (char **)alloc_2d_array(enc->height, enc->width,
		sizeof(char));
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->group[y][x] = 0;
			enc->org[y][x] = img->val[y][x];
		}
	}
	enc->org[enc->height][0] = (enc->maxval + 1) >> 1;
	enc->err[enc->height][0] = (enc->maxval + 1) >> 2;
	enc->uquant = (char **)alloc_2d_array(enc->num_class, MAX_UPARA + 1,
		sizeof(char));
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j <= MAX_UPARA; j++) {
			enc->uquant[i][j] = enc->num_group - 1;
		}
	}
	enc->econv = (int **)alloc_2d_array(enc->maxval+1, (enc->maxval<<1)+1,
		sizeof(int));
	enc->bconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
	enc->fconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
	enc->pmlist = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));
	enc->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	enc->spm.cumfreq = &(enc->spm.freq[MAX_SYMBOL]);

	enc->sigma = sigma_a;

	enc->mtfbuf = (int *)alloc_mem(enc->num_class * sizeof(int));

	enc->coef_m = (int *)alloc_mem( enc->prd_order * sizeof(int));
	for (i = 0; i < enc->prd_order; i++) {
		enc->coef_m[i] = 0;
	}
#if AUTO_PRD_ORDER
	enc->zero_m = (int *)alloc_mem(enc->prd_mhd * sizeof(int));
	for (j = 0; j < enc->prd_mhd; j++) {
		enc->zero_m[j] = NUM_ZMODEL >> 1;
	}
	enc->zero_fr = (int *)alloc_mem(NUM_ZMODEL * sizeof(int));
	for (i = 0; i < NUM_ZMODEL; i++) {
		enc->zero_fr[i] = (int)(zerocoef_prob[i] * (double)TOT_ZEROFR);
	}
	enc->coef_cost = (cost_t ***)alloc_2d_array(NUM_ZMODEL, 16, sizeof(cost_t *));
	ptr = (cost_t *)alloc_mem(NUM_ZMODEL * 16 * (enc->max_coef + 1) * sizeof(cost_t));
	for (i = 0; i < NUM_ZMODEL; i++) {
		for (j = 0; j < 16; j++) {
			enc->coef_cost[i][j] = ptr;
			ptr += (enc->max_coef + 1);
		}
	}
	for (k = 0; k < NUM_ZMODEL; k++) {
		for (i = 0; i < 16; i++) {
#if OPT_SIDEINFO
			double p, zero, nonz;
			uint cumb = 0;
			nonz = log((double)TOT_ZEROFR / (double)enc->zero_fr[k]) / log(2.0);
			zero = log((double)TOT_ZEROFR / ((double)TOT_ZEROFR - (double)enc->zero_fr[k])) / log(2.0);
			set_spmodel(&enc->spm, enc->max_coef + 1, i);
			cumb = enc->spm.freq[0];
			p = log(enc->spm.cumfreq[enc->max_coef + 1] - cumb);
			for (j = 1; j <= enc->max_coef; j++) {
				enc->coef_cost[k][i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
				enc->coef_cost[k][i][j] += 1.0 + nonz;
			}
			enc->coef_cost[k][i][0] = zero;
#else		
			for (j = 0; j <= enc->max_coef; j++) {
				enc->coef_cost[k][i][j] = 0;
			}
#endif
		}
	}
#else
	enc->coef_cost = (cost_t **)alloc_2d_array(16, enc->max_coef + 1,
		sizeof(cost_t));
	for (i = 0; i < 16; i++) {
#if OPT_SIDEINFO
		double p;
		set_spmodel(&enc->spm, enc->max_coef + 1, i);
		p = log(enc->spm.cumfreq[enc->max_coef + 1]);
		for (j = 0; j <= enc->max_coef; j++) {
			enc->coef_cost[i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
			if (j > 0) enc->coef_cost[i][j] += 1.0;
		}
#else
		for (j = 0; j <= enc->max_coef; j++) {
			enc->coef_cost[i][j] = 0;
		}
#endif	
	}
#endif
	enc->th_cost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	for (i = 0; i < MAX_UPARA + 2; i++) {
		enc->th_cost[i] = 0;
	}
	enc->class_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));
	c = log((double)enc->num_class) / log(2.0);
	for (i = 0; i < enc->num_class; i++) {
		enc->class_cost[i] = c;
	}
	for (i = 0; i < (QUADTREE_DEPTH << 3); i++) {
		enc->qtflag_cost[i] = 1.0;
	}
#if AUTO_DEL_CL
	i =  (enc->width + ( MAX_BSIZE >> QUADTREE_DEPTH) - 1) / ( MAX_BSIZE >> QUADTREE_DEPTH) 
		* (enc->height + ((MAX_BSIZE >> QUADTREE_DEPTH) -1 ) / (MAX_BSIZE >> QUADTREE_DEPTH));
	enc->err_cost = (cost_t **)alloc_2d_array(i, enc->num_class, sizeof(cost_t));
#endif
	enc->cl_hist = (int *)alloc_mem(enc->num_class * sizeof(int));

/////////////////////////////////////
/////以下テンプレートマッチング用////
/////////////////////////////////////

  enc->temp_num = (int **)alloc_2d_array(  enc->height , enc->width ,sizeof(int));
  enc->point_th = (POINT *)alloc_mem(sizeof(POINT) * enc->height * enc->width);
//  enc->cost_save = (double **)alloc_2d_array(  enc->height + 1 , enc->width ,sizeof(double));
//  enc->cost_save[enc->height][0] = 0;
//  enc->tm = (TempleteM_S***)alloc_3d_array( enc->height , enc->width , Y_SIZE * X_SIZE * 2 + X_SIZE  , sizeof(TempleteM_S) );
  enc->TM_U = (int **)alloc_2d_array(enc->width , enc->height , sizeof(int));
  enc->TM_gr = (int **)alloc_2d_array(enc->width , enc->height , sizeof(int));
  enc->TM_cn = (int *)alloc_mem(sizeof(int) * enc->num_pmodel);

  enc->PRD_gr = (int **)alloc_2d_array(enc->width , enc->height , sizeof(int));
//  enc->PRD_cn = (int *)alloc_mem(sizeof(int) * enc->num_pmodel);
  enc->array = (int ***)alloc_3d_array(enc->height , enc->width ,MAX_DATA_SAVE_DOUBLE ,sizeof(int));


	return (enc);
}

void init_class(ENCODER *enc)//選択情報
{
	int k, x, y, ly, lx, i, j, v, cl, sum, num_block;
	int *var, *tmp, **ptr;

	num_block = ((enc->height + BASE_BSIZE - 1) / BASE_BSIZE)
		* ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE);//BASE_BSIZE = 8.blockの数.64*64=4096

	var = (int *)alloc_mem(num_block * sizeof(int));
	ptr = (int **)alloc_mem(num_block * sizeof(int *));
	for (k = 0; k < num_block; k++) {
		y = (k / ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;//0〜8〜16〜24〜などが64個ずつ入る

		x = (k % ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;//0〜8〜16〜24〜などブロックの始まる座標が入る

		ly = y + BASE_BSIZE;//最後の画像外に出てしまう場合の対処
		if (ly > enc->height)
			ly = enc->height;
		lx = x + BASE_BSIZE;
		if (lx > enc->width)
			lx = enc->width;
		var[k] = sum = 0;//初期化
		for (i = y; i < ly; i++) {
			for (j = x; j < lx; j++) {//あるブロック内での計算
				v = enc->org[i][j];//輝度値
				sum += v;//輝度値の合計
				var[k] += v * v;//あるブロック内輝度値の2乗の合計
			}
		}
		var[k] -= sum * sum / ((ly -y) * (lx - x));//あるブロック内での分散の計算.平均を引いた

		ptr[k] = &(var[k]);//pointer配列に.計算の高速化のため?
	}//fin k

	/* sort */
	for (i = num_block - 1; i > 0; i--) {
		for (j = 0; j < i; j++) {
			if (*ptr[j] > *ptr[j + 1]) {
				tmp = ptr[j];
				ptr[j] = ptr[j + 1];
				ptr[j + 1] = tmp;
			}
		}
	}
	for (k = 0; k < num_block; k++) {
		cl = (k * enc->num_class) / num_block;//64個分ずつ0〜40までの数.512*512なら41個.256なら21個.これは固定
		v = (int)(ptr[k] - var);//var[0] = 404.この計算ホントにあってるの?

		y = (v / ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;

		x = (v % ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;

//printf("===(%d,%d)===\n",x,y);
		ly = y + BASE_BSIZE;
		if (ly > enc->height) ly = enc->height;
		lx = x + BASE_BSIZE;
		if (lx > enc->width) lx = enc->width;
		for (i = y; i < ly; i++) {
			for (j = x; j < lx; j++) {
				enc->class[i][j] = cl;
			}
		}
	}
	free(ptr);
	free(var);
}

void set_cost_model(ENCODER *enc, int f_mmse)
{
	int gr, i, j, k;
	double a, b, c, var;
	int gauss_index;
	PMODEL *pm;

	/* temporary repairs */
	if( enc->num_pmodel == 1) {
		gauss_index = 0;
	}else if( enc->num_pmodel == 8 || enc->num_pmodel == 16 || enc->num_pmodel == 32 ) {
		gauss_index = ( 5 * enc->num_pmodel) / 8 - 1;//gauss_index = 9
	}else {
		fprintf(stderr, "Cannot build Gaussian Distribution!\n");//ガウス分布
		exit(1);
	}

	for (i = 0; i <= enc->maxval; i++) {
		for (j = 0; j <= (enc->maxval << 1); j++) {//enc->maxval * 2
			k = (i << 1) - j - 1;//i * 2 - j - 1.
			if (k < 0) k = -(k + 1);//511個分の整数になる.例えば0〜510,3〜0・0〜506とか
			enc->econv[i][j] = k;//全てのパターンが入る
		}
	}
	enc->encval = enc->err;//err[enc->height][0]には256/4 .即ち64が入っている・・・.輝度値か?
	for (gr = 0; gr < enc->num_group; gr++) {
		var = enc->sigma[gr] * enc->sigma[gr];//varはσ＾2
		if (f_mmse) {//1st
			a = 0;
			b = 1.0;
		} else {
			a = 0.5 * log(2 * M_PI * var) / log(2.0);//log2(0.5 * log(2 * M_PI * var)) おそらく論文式(2-12)第一項分母だと思われる
			b = 1.0 / (2.0 * log(2.0) * var);//式(2-12)第二項からe^2を除いたもの
		}
		enc->pmlist[gr] = pm = enc->pmodels[gr][gauss_index];//gauss_indexはcn?
		for (k = 0; k < pm->size; k++) {//pm->size = 256
			c = (double)k * 0.5 + 0.25; //誤差eそのもの
			pm->cost[k] = (float)(a + b * c * c);
		}
		pm->subcost[0] = 0.0;
	}
	for (k = 0; k <= enc->maxprd; k++) {//初期化?
		enc->bconv[k] = 0;
		enc->fconv[k] = 0;
	}
	return;
}

void set_cost_rate(ENCODER *enc)
{
	int gr, k, i, j, mask, shift, num_spm;
	double a, c;
	PMODEL *pm;

	enc->encval = enc->org;//輝度値の代入
	mask = (1 << enc->pm_accuracy) - 1;//8 - 1.pm_accuracyは量子化精度
	shift = enc->coef_precision - enc->pm_accuracy;//coef_precision = 6 なのでshift = 3. coef_precisionは市街地距離かな
	for (k = 0; k <= enc->maxprd; k++) {//maxprd = enc->maxval << enc->coef_precision なので255 * 2^6 = 255 * 64 .とかかな

    //予測値は1/64で量子化されている(なので64倍)

		i = (enc->maxprd - k + (1 << shift) / 2) >> shift;//maxprdからkずつ引き 0.5 を足す.maxprd には予測値が入る.
		enc->fconv[k] = (i & mask);//bitのAND. 7以下しか入らない.量子化1/8の際場所はどこかを表す.
		enc->bconv[k] = (i >> enc->pm_accuracy);//i / 8　(1/64だが予測誤差は1/8なので　1/64と1/8がどのように関係しているか)
	}
	num_spm = 1 << enc->pm_accuracy;//8.基本的にnum_spmは8分割分のうちどれを選ぶかの値

	a = 1.0 / log(2.0);//log2(e).もしくはlog2に直すための係数
	for (gr = 0; gr < enc->num_group; gr++) {//分散16個分
		for (i = 0; i < enc->num_pmodel; i++) {//cn16個分
			pm = enc->pmodels[gr][i];
			for (j = 0; j < num_spm; j++) {//量子化精度8個分.TMではいらないねj=0
				for (k = 0; k < pm->size; k++) {//pm->size = 511
					pm->cost[k] = (float)(-a * log(pm->freq[k]));//log2(pm->freq[]).pm->freqは確率モデルの高さ(というか面積).pm->costには全領域のコストが入る
				}
				for (k = 0; k <= enc->maxval; k++) {
					c = pm->cumfreq[k + enc->maxval + 1] - pm->cumfreq[k];//pm->cumfreqは確率モデルを左からみた場合の全ての高さの合計が入っている
					pm->subcost[k] = (float)(a * log(c));//論文式(2-8)下式第一項
				}
				pm++;
			}
		}
	}
}

void predict_region(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
	int x, y, k, l, cl, prd, org;
	int *coef_p, *nzc_p;
	int *prd_p;
	int *roff_p, **roff_pp;
	int *err_p, *org_p;
	char *class_p;

	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		org_p = &enc->org[y][tlx];
		roff_pp = &enc->roff[y][tlx];
		err_p = &enc->err[y][tlx];
		prd_p = &enc->prd[y][tlx];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			roff_p = *roff_pp++;
			coef_p = enc->predictor[cl];
			nzc_p = enc->nzconv[cl];
			prd = 0;
			for (k = 0; k < enc->num_nzcoef[cl]; k++) {
				l = nzc_p[k];
				prd += org_p[roff_p[l]] * (coef_p[l]);
			}
			org = *org_p++;
			*prd_p++ = prd;
			prd = CLIP(0, enc->maxprd, prd);
			prd >>= (enc->coef_precision - 1);
			*err_p++ = enc->econv[org][prd];
		}
	}
}

int calc_uenc(ENCODER *enc, int y, int x)
{
	int u, k, *err_p, *roff_p, *wt_p;
	err_p = &enc->err[y][x];
	roff_p = enc->roff[y][x];
	wt_p = enc->ctx_weight;

	u = 0;
	for (k =0; k < NUM_UPELS; k++) {
		u += err_p[*roff_p++] * (*wt_p++);
	}
	u >>= 6;
	if (u > MAX_UPARA) u = MAX_UPARA;
	return (u);
}

cost_t calc_cost(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
	cost_t cost;
	int x, y, u, cl, gr, prd, e, base, frac;
	int *upara_p, *prd_p, *encval_p;
	char *class_p, *group_p;
	PMODEL *pm;

	if (bry > enc->height) bry = enc->height;
	if (tlx < 0) tlx = 0;
	if (brx > enc->width) brx = enc->width;
	cost = 0;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		group_p = &enc->group[y][tlx];
		upara_p = &enc->upara[y][tlx];
		encval_p = &enc->encval[y][tlx];
		prd_p = &enc->prd[y][tlx];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			*upara_p++ = u = calc_uenc(enc, y, x);
			*group_p++ = gr = enc->uquant[cl][u];
			e = *encval_p++;
			prd = *prd_p++;
			prd = CLIP(0, enc->maxprd, prd);
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost += pm->cost[base + e] + pm->subcost[base];//コストの計算はここで行う.pm->costはfreq.pm->subcostはcumfreq
		}
	}
	return (cost);
}

cost_t design_predictor(ENCODER *enc, int f_mmse)
{
	double **mat, *weight, w, e, d, pivot;
	int x, y, i, j, k, cl, gr, pivpos, *index, *roff_p, *org_p, *nzc_p;

//やっていることは自己相関関数．ただし時間ではなく距離


	mat = (double **)alloc_2d_array(enc->prd_order, enc->prd_order + 1, sizeof(double));
	index = (int *)alloc_mem(sizeof(int) * enc->prd_order);
	weight = (double *)alloc_mem(sizeof(double) * enc->num_group);
	for (gr = 0; gr < enc->num_group; gr++) {//gr は分散の番号16種類
		if (f_mmse) {//f_mmseは基本0
			weight[gr] = 1.0;
		} else {
			weight[gr] = 1.0 / (enc->sigma[gr] * enc->sigma[gr]);//当該分散番号のσ^2で割っている
		}
	}
	for (cl = 0; cl < enc->num_class; cl++) {//クラスの数だけ
		nzc_p = enc->nzconv[cl];
		for (i = 0; i < enc->prd_order; i++) {
			for (j = 0; j <= enc->prd_order; j++) {
				mat[i][j] = 0.0;//初期化
			}
		}
		for (y = 0; y < enc->height; y++) {//参考テンプレ．縦に動いた分
			for (x = 0; x < enc->width; x++) {//横に動いた分
				if (enc->class[y][x] != cl) {
					x += BASE_BSIZE - 1;//classのためのブロックの設定
					continue;
				}
				gr = enc->group[y][x];//gr は分散の番号16種類
				roff_p = enc->roff[y][x];//init_ref_offsetが入っている．予測器の範囲指定と番号付け
				org_p = &enc->org[y][x];//2次元配列が入る．enc->org = (int **)alloc_2d_array(enc->height+1, enc->width, sizeof(int));
				for (i = 0; i < enc->prd_order; i++) {
					w = weight[gr] * org_p[roff_p[nzc_p[i]]];//予測器内の画素に何をかけるか?また係数の作成．とその場所の指定．nzc_pで場所の番号(市街地距離)．roffで番号位置の参照．orgで輝度値の抽出．nzc_p[i] = enc->nzconv = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order, sizeof(int));
					for (j = i; j < enc->prd_order; j++) {
						mat[i][j] += w * org_p[roff_p[nzc_p[j]]];//1stループでは初めに予測値を周りの画素としている
                                                     //とりあえず保存
					}
					mat[i][enc->prd_order] += w * org_p[0];//符号化対象画素に重み掛けてる
                                                 //この時点ですでに偏微分(cost/a_k)の計算中である
				}
			}
		}
		for (i = 0; i < enc->prd_order; i++) {
			index[i] = i;//予測器とインデックスの関連付け
			for (j = 0; j < i; j++) {
				mat[i][j] = mat[j][i];//対称行列にする
			}
		}
    //まず縦方向でpivotを見つけて横方向で数値をいじる
		for (i = 0; i < enc->prd_order; i++) {
			pivpos = i;//絶対値の一番大きい場所
			pivot = fabs(mat[index[i]][i]);//浮動小数点の絶対値化
			for (k = i + 1; k < enc->prd_order; k++) {
				if (fabs(mat[index[k]][i]) > pivot) {
					pivot = fabs(mat[index[k]][i]);
					pivpos = k;
				}
			}
			k = index[i];
			index[i] = index[pivpos];
			index[pivpos] = k;
			if (pivot > 1E-10) {
				d = mat[index[i]][i];//対象画素なので自己相関を1とするため
				for (j = i; j <= enc->prd_order; j++) {
					mat[index[i]][j] /= d;//自己相関を1にする．正規化
				}
				for (k = 0; k < enc->prd_order; k++) {
					if (k == i) continue;
					d = mat[index[k]][i];
					for (j = i; j <= enc->prd_order; j++) {
						mat[index[k]][j] -= d * mat[index[i]][j];
					}
				}
			}
		}
		w = (1 << enc->coef_precision);
		e = 0.0;
		for (i = 0; i < enc->prd_order; i++) {
			if (fabs(mat[index[i]][i]) > 1E-10) {
				d = mat[index[i]][enc->prd_order] * w;
			} else {
				d = 0.0;
			}
			k = (int)d;
			if (k > d) k--;
			if (k < -enc->max_coef) {
				d = k = -enc->max_coef;
			} else if (k > enc->max_coef) {
				d = k = enc->max_coef;
			}
			enc->predictor[cl][nzc_p[i]] = k;
			d -= k;
			e += d;
			mat[index[i]][enc->prd_order] = d;
		}
		/* minimize mean rounding errors */
		k = (int)(e + 0.5);
		for (;k > 0; k--) {
			d = 0;
			for (j = i = 0; i < enc->prd_order; i++) {
				if (mat[index[i]][enc->prd_order] > d) {
					d = mat[index[i]][enc->prd_order];
					j = i;
				}
			}
			if (enc->predictor[cl][nzc_p[j]] < enc->max_coef) enc->predictor[cl][nzc_p[j]]++;
			mat[index[j]][enc->prd_order] = 0;
		}
	}
	free(weight);
	free(index);
	free(mat);

	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}


cost_t optimize_group(ENCODER *enc)
{
	cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
	int x, y, th1, th0, k, u, cl, gr, prd, e, base, frac;
	int **trellis;
	PMODEL *pm, **pm_p;

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(cost_t));
	thc_p = enc->th_cost;
	for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
	/* Dynamic programming */
	for (cl = 0; cl < enc->num_class; cl++) {
		//if (enc->cl_hist[cl] == 0) continue;
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			for (u = 0; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] = 0;
			}
		}
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->class[y][x] == cl) {
					u = enc->upara[y][x] + 1;
					e = enc->encval[y][x];
					prd = enc->prd[y][x];
					prd = CLIP(0, enc->maxprd, prd);
					base = enc->bconv[prd];
					frac = enc->fconv[prd];
					pm_p = enc->pmlist;
					for (gr = 0; gr < enc->num_group; gr++) {
						pm = (*pm_p++) + frac;
						cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
					}
				}
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] += cbuf_p[u - 1];
			}
		}
		cbuf_p = cbuf[0];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u] + thc_p[u];
		}
		for (gr = 1; gr < enc->num_group - 1; gr++) {
			cbuf_p = cbuf[gr];
			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}
				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;
			}
		}

		cbuf_p = cbuf[gr];
		/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
		th1 = MAX_UPARA + 1;
		th0 = th1;
		min_cost = dpcost[th1] - cbuf_p[th1];
		for (k = 0; k < th1; k++) {
			cost = dpcost[k] - cbuf_p[k];
			if (cost < min_cost) {
				min_cost = cost;
				th0 = k;
			}
		}
		trellis[gr][th1] = th0;

		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
			enc->th[cl][gr - 1] = th1;
		}
	}
	/* set context quantizer */
	for (cl = 0; cl < enc->num_class; cl++) {
		u = 0;
		for (gr = 0; gr < enc->num_group; gr++) {
			for (; u < enc->th[cl][gr]; u++) {
				enc->uquant[cl][u] = gr;
			}
		}
	}
	/* renew groups */
	cost = 0;
	pm_p = enc->pmlist;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			cl = enc->class[y][x];
			u = enc->upara[y][x];
			enc->group[y][x] = gr = enc->uquant[cl][u];
			e = enc->encval[y][x];
			prd = enc->prd[y][x];
			prd = CLIP(0, enc->maxprd, prd);
			base = enc->bconv[prd];
			pm = pm_p[gr] + enc->fconv[prd];
			cost += pm->cost[base + e] + pm->subcost[base];
		}
	}
	/* optimize probability models */
	if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
		if (enc->num_pmodel > MAX_UPARA + 2) {
			free(cbuf);
			cbuf = (cost_t **)alloc_2d_array(enc->num_group, enc->num_pmodel,
				sizeof(cost_t));
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			for (k = 0; k < enc->num_pmodel; k++) {
				cbuf[gr][k] = 0;
			}
		}
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[y][x];
				e = enc->encval[y][x];
				prd = enc->prd[y][x];
				prd = CLIP(0, enc->maxprd, prd);
				base = enc->bconv[prd];
				frac = enc->fconv[prd];
				for (k = 0; k < enc->num_pmodel; k++) {
					pm = enc->pmodels[gr][k] + frac;
					cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
				}
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			pm = enc->pmodels[gr][0];
			cost = cbuf[gr][0];
			for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[gr][k];
				}
			}
			pm_p[gr] = pm;
		}
		cost = 0.0;
		for (gr = 0; gr < enc->num_group; gr++) {
			cost += cbuf[gr][pm_p[gr]->id];
		}
	}
	free(cbuf);
	free(dpcost);
	free(trellis);
	return (cost);
}

void set_prdbuf(ENCODER *enc, int **prdbuf, int **errbuf,
				int tly, int tlx, int bufsize)
{
	int x, y, brx, bry, cl, k, l, prd, *prdbuf_p, *errbuf_p, *coef_p;
	int buf_ptr, org, *org_p, *roff_p, *nzc_p;

	brx = (tlx + bufsize < enc->width) ? (tlx + bufsize) : enc->width;
	bry = (tly + bufsize < enc->height) ? (tly + bufsize) : enc->height;
	for (cl = 0; cl < enc->num_class; cl++) {
		buf_ptr = bufsize * (tly % bufsize) + tlx % bufsize;
		nzc_p = enc->nzconv[cl];
		for (y = tly; y < bry; y++) {
			prdbuf_p = &prdbuf[cl][buf_ptr];
			errbuf_p = &errbuf[cl][buf_ptr];
			buf_ptr += bufsize;
			org_p = &enc->org[y][tlx];
			for (x = tlx; x < brx; x++) {
				if (cl == enc->class[y][x]) {
					*prdbuf_p++ = enc->prd[y][x];
					*errbuf_p++ = enc->err[y][x];
					org_p++;
				} else {
					coef_p = enc->predictor[cl];
					roff_p = enc->roff[y][x];
					prd = 0;
					for (k = 0; k < enc->num_nzcoef[cl]; k++) {
						l = nzc_p[k];
						prd += org_p[roff_p[l]] * (coef_p[l]);
					}
					org = *org_p++;
					*prdbuf_p++ = prd;
					prd = CLIP(0, enc->maxprd, prd);
					prd >>= (enc->coef_precision - 1);
					*errbuf_p++ = enc->econv[org][prd];
				}
			}
		}
	}
}

int find_class(ENCODER *enc, int **prdbuf, int **errbuf,
			   int tly, int tlx, int bry, int brx, int bufsize, cost_t *err_cost)
{
	cost_t cost, min_cost;
	int x, y, bufptr, cl, min_cl;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	min_cost = 1E8;
	min_cl = 0;
	for	(cl = 0; cl < enc->num_class; cl++) {
		bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
		for (y = tly; y < bry; y++) {
			class_p = &enc->class[y][tlx];
			prd_p = &enc->prd[y][tlx];
			prdbuf_p = &prdbuf[cl][bufptr];
			err_p = &enc->err[y][tlx];
			errbuf_p = &errbuf[cl][bufptr];
			bufptr += bufsize;
			for (x = tlx; x < brx; x++) {
				*class_p++ = cl;
				*prd_p++ = *prdbuf_p++;
				*err_p++ = *errbuf_p++;
			}
		}
		cost = calc_cost(enc, tly, tlx, bry, brx);
		err_cost[cl] = cost;
		if( enc->optimize_loop == 2) {
			cost += enc->class_cost[enc->mtfbuf[cl]];
		}
		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}
	bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		prd_p = &enc->prd[y][tlx];
		prdbuf_p = &prdbuf[min_cl][bufptr];
		err_p = &enc->err[y][tlx];
		errbuf_p = &errbuf[min_cl][bufptr];
		bufptr += bufsize;
		for (x = tlx; x < brx; x++) {
			*class_p++ = min_cl;
			*prd_p++ = *prdbuf_p++;
			*err_p++ = *errbuf_p++;
		}
	}
	return (min_cl);
}

cost_t vbs_class(ENCODER *enc, int **prdbuf, int **errbuf, int tly, int tlx,
				 int blksize, int width, int level, int *blk)
{
	int y, x, k, bry, brx, cl, bufsize, bufptr, ctx, s_blk;
	int mtf_save[MAX_CLASS];
	char **qtmap;
	cost_t cost1, cost2, qtcost, *err_cost;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	s_blk = *blk;
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		bufsize = MAX_BSIZE;
	} else {
		bufsize = BASE_BSIZE;
	}
	brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
	bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;
	if (tlx >= brx || tly >= bry) return (0);
	for (k = 0; k < enc->num_class; k++) {
		mtf_save[k] = enc->mtfbuf[k];
	}
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx,
		blksize, width, enc->num_class);
	err_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));
	cl = find_class(enc, prdbuf, errbuf, tly, tlx, bry, brx, bufsize, err_cost);
	qtcost = enc->class_cost[enc->mtfbuf[cl]];
	if (level > 0) {
		/* context for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;
		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}
		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;
		ctx = ((level - 1) * 4 + ctx) << 1;
		/* Quad-tree partitioning */
		cost1 = calc_cost(enc, tly, tlx, bry, brx)
			+ enc->class_cost[enc->mtfbuf[cl]] + enc->qtflag_cost[ctx];
		blksize >>= 1;
		for (k = 0; k < enc->num_class; k++) {
			enc->mtfbuf[k] = mtf_save[k];
		}
		qtcost = enc->qtflag_cost[ctx + 1];
		qtcost += vbs_class(enc, prdbuf, errbuf, tly, tlx,
			blksize, width, level - 1, blk);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly, tlx+blksize,
			blksize, width, level - 1, blk);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly+blksize, tlx,
			blksize, width, level - 1, blk);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly+blksize, tlx+blksize,
			blksize, brx, level - 1, blk);
		cost2 = calc_cost(enc, tly, tlx, bry, brx) + qtcost;
		if (cost1 < cost2) {
			blksize <<= 1;
			for (k = 0; k < enc->num_class; k++) {
				enc->mtfbuf[k] = mtf_save[k];
			}
			mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx,
				blksize, width, enc->num_class);
			qtcost = enc->class_cost[enc->mtfbuf[cl]]
			+ enc->qtflag_cost[ctx];
			bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
			for (y = tly; y < bry; y++) {
				class_p = &enc->class[y][tlx];
				prd_p = &enc->prd[y][tlx];
				prdbuf_p = &prdbuf[cl][bufptr];
				err_p = &enc->err[y][tlx];
				errbuf_p = &errbuf[cl][bufptr];
				bufptr += bufsize;
				for (x = tlx; x < brx; x++) {
					*class_p++ = cl;
					*prd_p++ = *prdbuf_p++;
					*err_p++ = *errbuf_p++;
				}
			}
			tly = (tly / MIN_BSIZE) >> level;
			tlx = (tlx / MIN_BSIZE) >> level;
			bry = tly + 1;
			brx = tlx + 1;
			for (; level > 0; level--) {
				qtmap = enc->qtmap[level - 1];
				for (y = tly; y < bry; y++) {
					for (x = tlx; x < brx; x++) {
						qtmap[y][x] = 0;
					}
				}
				tly <<= 1;
				tlx <<= 1;
				bry <<= 1;
				brx <<= 1;
			}
#if AUTO_DEL_CL
			*blk = s_blk + 1;
			for (cl = 0; cl < enc->num_class; cl++) {
				enc->err_cost[s_blk][cl] = err_cost[cl];
			}
#endif
		} else {
			qtmap[y][x] = 1;
		}
	} else {
#if AUTO_DEL_CL
		*blk = s_blk + 1;
		for (cl = 0; cl < enc->num_class; cl++) {
			enc->err_cost[s_blk][cl] = err_cost[cl];
		}
#endif
	}
	free(err_cost);
	return (qtcost);
}

void count_cl(ENCODER *enc)
{
	int cl, y, x;

	for (cl = 0; cl < enc->num_class; cl++) enc->cl_hist[cl] = 0; 
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->cl_hist[(int)enc->class[y][x]]++;
		}
	}
}

cost_t optimize_class(ENCODER *enc)
{
	int y, x, i, blksize, level, blk = 0;
	int **prdbuf, **errbuf;

	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE;
	} else {
		level = 0;
		blksize = BASE_BSIZE;
	}
	for (i = 0; i < enc->num_class; i++) {
		enc->mtfbuf[i] = i;
	}
	prdbuf =(int **)alloc_2d_array(enc->num_class, blksize * blksize,
		sizeof(int));
	errbuf =(int **)alloc_2d_array(enc->num_class, blksize * blksize,
		sizeof(int));
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			set_prdbuf(enc, prdbuf, errbuf, y, x, blksize);
			vbs_class(enc, prdbuf, errbuf, y, x,
				blksize, enc->width, level, &blk);
		}
	}
	free(errbuf);
	free(prdbuf);
	count_cl(enc);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

#if AUTO_PRD_ORDER

/*****************************************************************************
enc->nzconv		: Non-Zero Coef to The Order (Index) table.
enc->num_nzcoef : Number of Non-Zero Coef
*****************************************************************************/
void set_prd_pels(ENCODER *enc)
{
	int cl, k, i;

	for (cl = 0; cl < enc->num_class; cl++) {
		k = 0;
		for (i = 0; i < enc->max_prd_order; i++) {
			if (enc->predictor[cl][i] != 0) {
				enc->nzconv[cl][k] = i;
				k++;
			}
		}
		for (i = k; i < enc->max_prd_order; i++) {
			enc->nzconv[cl][i] = (int)1E5;			//for bugfix
		}
		enc->num_nzcoef[cl] = k;
	}
}

/*****************************************************************************
Coef Modification
*****************************************************************************/
void optimize_coef(ENCODER *enc, int cl, int pos, int *num_eff)
{
#define S_RANGE 2		//Search Range  ex. 2 -> +-1
#define SUBS_RANGE 2	//Search Range of One Coef Modification
	cost_t *cbuf, *cbuf_p, *cbuf_p2, c_base, *coef_cost_p1, *coef_cost_p2;
	int i, j, k, l, x, y, df1, df2, df_f, *org_p, base, *roff_p;
	int coef_abs_pos, coef_abs, coef_pos;
	int num_swap_search, *swap_search;
	int prd, shift, maxprd, *coef_p, *nzc_p, prd_c;
	char *class_p, onoff;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;
	int diffy, diffx, *diff, *diff_p, *diff_p2, sgn;

	num_swap_search = enc->max_prd_order - enc->num_nzcoef[cl];
	cbuf = (cost_t *)alloc_mem((SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search) * sizeof(cost_t));
	diff = (int *) alloc_mem ((SUBS_RANGE + (enc->max_prd_order * S_RANGE)+ num_swap_search) * sizeof(int));
	//table of exchange coef
	swap_search = (int *)alloc_mem(num_swap_search * sizeof(int));
	nzc_p = enc->nzconv[cl];
	coef_p = enc->predictor[cl];
	coef_pos = coef_p[pos];
	coef_abs_pos = (coef_pos < 0) ? -coef_pos : coef_pos;
	cbuf_p = cbuf;
	cbuf_p2 = &cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
	diff_p = diff;
	diff_p2 = &diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
	i = enc->ord2mhd[pos];
	coef_cost_p1 = enc->coef_cost[enc->zero_m[i]][enc->coef_m[i]];
	for (i = 0; i < (SUBS_RANGE >> 1); i++) {
		df1 = coef_pos - (i + 1);
		for (j = 0; j < 2; j++) {
			y = df1;
			sgn = 0;
			if (y < 0) {
				y = -y;
				sgn = 1;
			}
			if (y > enc->max_coef) y = enc->max_coef;
			*cbuf_p++ = coef_cost_p1[y];
			*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;
			df1 += (i + 1) << 1;
		}
	}
	k = 0;
	for (i = 0; i < enc->max_prd_order; i++) {
		onoff = (coef_p[i] == 0) ? 1 : 0;
		j = enc->ord2mhd[i];
		coef_cost_p2 = enc->coef_cost[enc->zero_m[j]][enc->coef_m[j]];
		df2 = coef_p[i];
		coef_abs = (df2 < 0) ? -df2 : df2;
		c_base = coef_cost_p2[coef_abs];
		for (l = 0; l < (S_RANGE >> 1); l++) {
			df1 = coef_pos - (l + 1);
			df2 = coef_p[i] + (l + 1);
			for (j = 0; j < 2; j++) {		//change "+ or -" modification
				y = df1;
				x = df2;
				sgn = 0;
				if (y < 0) {
					y = -y;
					sgn = 1;
				}
				diffy = y - enc->max_coef;
				if (diffy < 0) diffy = 0;
				if (y > enc->max_coef) y = enc->max_coef;
				if (x < 0) x = -x;
				diffx = x - enc->max_coef;
				if (diffx < 0) diffx = 0;
				if (x > enc->max_coef) x = enc->max_coef;
				if (diffy > 0 || diffx > 0) {
					if (diffy > diffx) {
						x += (j) ? diffy - diffx : diffx - diffy;
						if (x < 0) {
							x = -x;
						}
					} else {
						y += (j) ? diffy - diffx : diffx - diffy;
						if (y < 0) {
							y = -y;
							sgn = (sgn) ? 0 : 1;
						}
					}
				}
				*cbuf_p++ = coef_cost_p1[y] + coef_cost_p2[x] - c_base;
				*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;
				df1 += (l + 1) << 1;
				df2 -= (l + 1) << 1;
			}
		}
		if (onoff == 1) {
			*cbuf_p2++ = coef_cost_p1[0] + coef_cost_p2[coef_abs_pos] - coef_cost_p2[0];
			*diff_p2++ = -coef_pos;
			swap_search[k++] = i;
		}
	}
	for (i = 0; i < S_RANGE; i++) {
		cbuf[SUBS_RANGE + pos * S_RANGE + i] = coef_cost_p1[coef_abs_pos];
		diff[SUBS_RANGE + pos * S_RANGE + i] = 0;
	}
	// before here -> the calculation of side cost
	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;
	for (y = 0; y < enc->height; y++) {
		class_p = enc->class[y];
		for (x = 0; x < enc->width; x++) {
			if (cl != *class_p++) continue;
			roff_p = enc->roff[y][x];
			prd = enc->prd[y][x];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos]];
			cbuf_p = cbuf;
			cbuf_p2 = &cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
			diff_p = diff;
			diff_p2 = &diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
			df_f = df1;
			//Only One Coef Modification
			for (i = 0; i < (SUBS_RANGE >> 1); i++) {
				for (j = 0; j < 2; j++) {	//change "+ or -" modification
					prd_c = prd + df_f * (*diff_p++);
					prd_c = CLIP(0, maxprd, prd_c);
					base = bconv_p[prd_c];
					pm = pm_p + fconv_p[prd_c];
					*cbuf_p++ += pm->cost[*org_p + base]
					+ pm->subcost[base];
				}
			}
			//Modification and Search The Other Coef
			for (i = 0; i < enc->max_prd_order; i++) {
				df2 = org_p[roff_p[i]];
				df_f = df1 - df2;
				for (l = 0; l < (S_RANGE >> 1); l++) {
					for (j = 0; j < 2; j++) {	//change "+ or -"
						prd_c = prd + df_f * (*diff_p++);
						prd_c = CLIP(0, maxprd, prd_c);
						base = bconv_p[prd_c];
						pm = pm_p + fconv_p[prd_c];
						*cbuf_p++ += pm->cost[*org_p + base]
						+ pm->subcost[base];
					}
				}
			}
			//exchange two coefs
			for (j = 0; j < num_swap_search; j++) {
				k = swap_search[j];
				prd_c = prd + (df1 - org_p[roff_p[k]]) * (*diff_p2++);
				prd_c = CLIP(0, maxprd, prd_c);
				base = bconv_p[prd_c];
				pm = pm_p + fconv_p[prd_c];
				*cbuf_p2++ += pm->cost[*org_p + base]
				+ pm->subcost[base];
			}
		}
	}
	j = SUBS_RANGE + pos * S_RANGE;
	for (i = 0; i < SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search; i++) {
		if (cbuf[i] < cbuf[j]) {
			j = i;
		}
	}
	free(cbuf);
	if (j == SUBS_RANGE + pos * S_RANGE) {
		free(swap_search);
		free(diff);
		return;
	}
	if (j < SUBS_RANGE) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					org_p = &enc->org[y][x];
					roff_p = enc->roff[y][x];
					enc->prd[y][x] += org_p[roff_p[pos]] * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
	} else {
		if (j < SUBS_RANGE + (enc->max_prd_order * S_RANGE)) {
			i = j - SUBS_RANGE;
			i /= S_RANGE;
		} else {
			i = j - SUBS_RANGE - (enc->max_prd_order * S_RANGE);
			i = swap_search[i];
		}
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					org_p = &enc->org[y][x];
					roff_p = enc->roff[y][x];
					enc->prd[y][x] += (org_p[roff_p[pos]] - org_p[roff_p[i]]) * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
		coef_p[i] -= diff[j];
	}
	free(swap_search);
	free(diff);
	if (diff[j] != 0) (*num_eff)++;
}

cost_t optimize_predictor(ENCODER *enc)
{
	int cl, pos, k, num_nzc, num_eff;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif

	for (cl = 0; cl < enc->num_class; cl++) {
		num_nzc = enc->num_nzcoef[cl];
		num_eff = 0;
		if (enc->cl_hist[cl] == 0) continue;
		for (k = 0; k < enc->num_search[cl]; k++) {
			if (enc->num_nzcoef[cl] == 0) continue;
			pos = (int)(((double)rand() * enc->num_nzcoef[cl]) / (RAND_MAX+1.0));
			pos = enc->nzconv[cl][pos];
			optimize_coef(enc, cl, pos, &num_eff);
			set_prd_pels(enc);
		}
		enc->num_search[cl] = num_eff + 3;
	}
	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

#else

void optimize_coef(ENCODER *enc, int cl, int pos1, int pos2)
{
#define SEARCH_RANGE 11
#define SUBSEARCH_RANGE 3
	cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
	int i, j, k, x, y, df1, df2, base;
	int prd, prd_f, shift, maxprd, *coef_p, *roff_p, *org_p;
	char *class_p;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;

	cbuf_p = cbuf;
	coef_p = enc->predictor[cl];
	k = 0;
	for (i = 0; i < SEARCH_RANGE; i++) {
		y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);
		if (y < 0) y = -y;
		if (y > enc->max_coef) y = enc->max_coef;
		for (j = 0; j < SUBSEARCH_RANGE; j++) {
			x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1))
				- (j - (SUBSEARCH_RANGE >> 1));
			if (x < 0) x = -x;
			if (x > enc->max_coef) x = enc->max_coef;

			cbuf_p[k++] = enc->coef_cost[enc->coef_m[pos1]][y]
			+ enc->coef_cost[enc->coef_m[pos2]][x];

		}
	}
	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;
	for (y = 0; y < enc->height; y++) {
		class_p = enc->class[y];
		for (x = 0; x < enc->width; x++) {
			if (cl != *class_p++) continue;
			roff_p = enc->roff[y][x];
			prd = enc->prd[y][x];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos1]];
			df2 = org_p[roff_p[pos2]];
			prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1)
				+ df2 * (SUBSEARCH_RANGE >> 1);
			cbuf_p = cbuf;
			for (i = 0; i < SEARCH_RANGE; i++) {
				for (j = 0; j < SUBSEARCH_RANGE; j++) {
					prd = prd_f;
					prd = CLIP(0, maxprd, prd);
					base = bconv_p[prd];
					pm = pm_p + fconv_p[prd];
					(*cbuf_p++) += pm->cost[*org_p + base]
					+ pm->subcost[base];
					prd_f -= df2;
				}
				prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
			}
		}
	}
	cbuf_p = cbuf;
	j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;
	for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
		if (cbuf_p[i] < cbuf_p[j]) {
			j = i;
		}
	}
	i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
	j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
	y = coef_p[pos1] + i;
	x = coef_p[pos2] - i - j;
	if (y < -enc->max_coef) y = -enc->max_coef;
	else if (y > enc->max_coef) y = enc->max_coef;
	if (x < -enc->max_coef) x = -enc->max_coef;
	else if (x > enc->max_coef) x = enc->max_coef;
	i = y - coef_p[pos1];
	j = x - coef_p[pos2];
	if (i != 0 || j != 0) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					roff_p = enc->roff[y][x];
					org_p = &enc->org[y][x];
					enc->prd[y][x] += org_p[roff_p[pos1]] * i
						+ org_p[roff_p[pos2]] * j;
				}
			}
		}
		coef_p[pos1] += i;
		coef_p[pos2] += j;
	}
}

cost_t optimize_predictor(ENCODER *enc)
{
	int cl, k, pos1, pos2;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif

	for (cl = 0; cl < enc->num_class; cl++) {
		//if (enc->cl_hist[cl] == 0) continue;
		for (k = 0; k < enc->max_prd_order; k++) {
retry:
			pos1 = (int)(((double)rand() * enc->max_prd_order) / (RAND_MAX+1.0));
			pos2 = (int)(((double)rand() * enc->max_prd_order) / (RAND_MAX+1.0));
			if (pos1 == pos2) goto retry;
			optimize_coef(enc, cl, pos1, pos2);
		}
	}
	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

#endif

int putbits(FILE *fp, int n, uint x)
{
	static int bitpos = 8;
	static uint bitbuf = 0;
	int bits;

	bits = n;
	if (bits <= 0) return (0);
	while (n >= bitpos) {
		n -= bitpos;
		if (n < 32) {
			bitbuf |= ((x >> n) & (0xff >> (8 - bitpos)));
		}
		putc(bitbuf, fp);
		bitbuf = 0;
		bitpos = 8;
	}
	bitpos -= n;
	bitbuf |= ((x & (0xff >> (8 - n))) << bitpos);
	return (bits);
}

void remove_emptyclass(ENCODER *enc)
{
	int cl, i, k, x, y;

	for (cl = 0; cl < enc->num_class; cl++) {
		enc->mtfbuf[cl] = 0;
	}
	for (y = 0; y < enc->height; y += MIN_BSIZE) {
		for (x = 0; x < enc->width; x += MIN_BSIZE) {
			cl = enc->class[y][x];
			enc->mtfbuf[cl]++;
		}
	}
	for (i = cl = 0; i < enc->num_class; i++) {
		if (enc->mtfbuf[i] == 0) {
			enc->mtfbuf[i] = -1;
		} else {
			enc->mtfbuf[i] = cl++;
		}
	}
	if (cl == enc->num_class) return;	/* no empty class */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			i = enc->class[y][x];
			enc->class[y][x] = enc->mtfbuf[i];
		}
	}
	for (i = cl = 0; i < enc->num_class; i++) {
		if (enc->mtfbuf[i] < 0) continue;
		if (cl != i) {
			for (k = 0; k < enc->max_prd_order; k++) {
				enc->predictor[cl][k] = enc->predictor[i][k];
			}
			for (k = 0; k < enc->num_group - 1; k++) {
				enc->th[cl][k] = enc->th[i][k];
			}
		}
		cl++;
	}
	printf("M = %d\n", cl);
	enc->num_class = cl;
}

int write_header(ENCODER *enc, FILE *fp)
{
	int bits;

	bits = putbits(fp, 16, MAGIC_NUMBER);
	bits += putbits(fp, 8, VERSION);
	bits += putbits(fp, 16, enc->width);
	bits += putbits(fp, 16, enc->height);
	bits += putbits(fp, 16, enc->maxval);
//	bits += putbits(fp, 6, enc->num_class);
	bits += putbits(fp, 6, enc->num_group);
	bits += putbits(fp, 7, enc->max_prd_order);
	bits += putbits(fp, 6, enc->num_pmodel - 1);
	bits += putbits(fp, 4, enc->coef_precision - 1);
	bits += putbits(fp, 3, enc->pm_accuracy);
//	bits += putbits(fp, 1, (enc->quadtree_depth < 0)? 0 : 1);
	return (bits);
}

void set_qtindex(ENCODER *enc, int *index, uint *hist, int *numidx,
				 int tly, int tlx, int blksize, int width, int level)
{
	int i, cl, x, y, ctx;
	char **qtmap;

	if (tly >= enc->height || tlx >= enc->width) return;//範囲外なら演算しない
	if (level > 0) {//四分木を使用していないのでスルー
		/* context modeling for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[level - 1];
		y = ((tly + MIN_BSIZE - 1) / MIN_BSIZE) >> level;//levelによってブロックがどこまでか決める
		x = ((tlx + MIN_BSIZE - 1) / MIN_BSIZE) >> level;
		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (tlx + blksize < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}
		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;
		ctx = ((level - 1) * 4 + ctx) << 1;
		if (qtmap[y][x] == 1) {
			ctx++;
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[ctx]++;
			blksize >>= 1;
			set_qtindex(enc, index, hist, numidx, tly, tlx,
				blksize, width, level - 1);
			set_qtindex(enc, index, hist, numidx, tly, tlx + blksize,
				blksize, width, level - 1);
			set_qtindex(enc, index, hist, numidx, tly + blksize, tlx,
				blksize, width, level - 1);
			width = tlx + blksize * 2;
			if (width >= enc->width) width = enc->width;
			set_qtindex(enc, index, hist, numidx, tly + blksize, tlx + blksize,
				blksize, width, level - 1);
			return;
		} else {
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[ctx]++;
		}
	}//level fin
	cl = enc->class[tly][tlx];
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx, blksize, width, enc->num_class);
	i = enc->mtfbuf[cl];//iには順番が入る
	index[(*numidx)++] = i;//indexに順番を保存
	hist[i]++;//ヒストグラムの作成
	return;
}

int encode_class(FILE *fp, ENCODER *enc, int flag)
{
	int i, j, k, numidx, blksize, level, x, y, ctx, bits, *index;
	uint *hist;
	cost_t cost;
	PMODEL *pm;
	double p, c;
	int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];
	cost_t qtflag_cost[QUADTREE_DEPTH << 3], *class_cost;
//#if (!OPT_SIDEINFO)
//	if (fp == NULL) return(0);
//#endif
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;//四文木の情報
		blksize = MAX_BSIZE;//32
		numidx = 0;

		y = (enc->height + MIN_BSIZE - 1) / MIN_BSIZE;//MIN_BSIZE = 4
		x = (enc->width + MIN_BSIZE - 1) / MIN_BSIZE;
		for (k = 0; k <= level; k++) {
			numidx += y * x;
			y = (y + 1) >> 1;
			x = (x + 1) >> 1;
		}
		for (k = 0; k < QUADTREE_DEPTH << 3; k++) {
			enc->qtctx[k] = 0;
		}
	} else {//基本型.四分木は使用してないのでこちらを使用
		level = 0;

		blksize = BASE_BSIZE;
		numidx = ((enc->height + (BASE_BSIZE - 1)) / BASE_BSIZE)
			* ((enc->width + (BASE_BSIZE - 1)) / BASE_BSIZE);//numidx = 4096.when 512*512

	}
	hist = (uint *)alloc_mem(enc->num_class * sizeof(uint));
	index = (int *)alloc_mem(numidx * sizeof(int));
	for (i = 0; i < enc->num_class; i++) {//0〜41
		hist[i] = 0;
		enc->mtfbuf[i] = i;
	}
	numidx = 0;
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			set_qtindex(enc, index, hist, &numidx, y, x,
				blksize, enc->width, level);//ここでmtfの計算を行っていた.indexに値が保存されている
		}
	}
	bits = 0;
	/* Arithmetic */
	class_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));
	/* context modeling for quad-tree flag */
	if (level > 0) {//ここは使用していない

		for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
			cost = INT_MAX;
			for (i = k = 0; i < 7; i++) {
				p = qtree_prob[i];
				c = -log(p) * (cost_t)enc->qtctx[(ctx << 1) + 1] 
				-log(1.0 - p) * (cost_t)enc->qtctx[ctx << 1];
				if (c < cost) {
					k = i;
					cost = c;
				}
			}
			p = qtree_prob[qtree_code[ctx] = k];
			qtflag_cost[(ctx << 1) + 1] = -log(p) / log(2.0);
			qtflag_cost[ctx << 1] = -log(1.0 - p) / log(2.0);
		}
	}

	/* quantization of log-transformed probability */
	c = 0.0;
	for (i = 0; i < enc->num_class; i++) {
		c += (double)hist[i];

	}
	for (i = 0; i < enc->num_class; i++) {//Move To Front
		p = (double)hist[i] / c;//c = 4096.hist[i]は順番のヒストグラムかと思われる

		if (p > 0.0) {
			mtf_code[i] = (int)(-log(p) / log(2.0) * (PMCLASS_LEVEL / PMCLASS_MAX));//PMCLASS_MAX = 16
			if (mtf_code[i] >= PMCLASS_LEVEL) {//PMCLASS_LEVEL = 32
				mtf_code[i] = PMCLASS_LEVEL - 1;
			}
		} else {
			mtf_code[i] = PMCLASS_LEVEL - 1;
		}//mtf_codeは頻度分布を再構成するために必要
    
		p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
			* PMCLASS_MAX / PMCLASS_LEVEL);//確率モデルの山原型かな

		class_cost[i] = -log(p) / log(2.0);
		hist[i] = (uint)(p * (1 << 10));

		if (hist[i] <= 0) hist[i] = 1;//確率モデルの1補完
	}
	if (fp == NULL) {
		cost = 0.0;
		for (j = 0; j < numidx; j++) {
			i = index[j];
			if (i < 0) {
				ctx = -(i + 1);
				cost += qtflag_cost[ctx];
			} else {
				cost += class_cost[i];
			}
		}
		if (flag == 1) {
			for (i = 0; i < (QUADTREE_DEPTH << 3); i++)
				enc->qtflag_cost[i] = qtflag_cost[i];
			for (i = 0; i < enc->num_class; i++)
				enc->class_cost[i] = class_cost[i];
		}
		bits = (int)cost;
	} else {	/* actually encode */
		PMODEL cpm[1];
		/* additional info. */
		pm = &enc->spm;
		if (level > 0) {//四分木を使用していないためスルー
			set_spmodel(pm, 7, -1);
			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				i = qtree_code[ctx];
				rc_encode(fp, enc->rc, pm->cumfreq[i], pm->freq[i],
					pm->cumfreq[pm->size],0);
			}
		}//fin level
		set_spmodel(pm, PMCLASS_LEVEL, -1);
		for (i = 0; i < enc->num_class; i++) {
			j = mtf_code[i];//頻度分布を再生するために必要

			rc_encode(fp, enc->rc, pm->cumfreq[j], pm->freq[j],
				pm->cumfreq[pm->size],0);//ここで画像毎の頻度分布を送っている
			if (pm->cumfreq[pm->size] < MAX_TOTFREQ) {
				for (j = 0; j < pm->size; j++) {
					if (j < mtf_code[i]) {
						pm->freq[j] /= 2;
					} else {
						pm->freq[j] *= 2;
					}
					if (pm->freq[j] <= 0) pm->freq[j] = 1;
					pm->cumfreq[j + 1] = pm->cumfreq[j] + pm->freq[j];
				}
			}
		}
		/* set prob. models */
		if (level > 0) {
			pm->size = QUADTREE_DEPTH << 3;
			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				i = qtree_code[ctx];
				p = qtree_prob[i];
				pm->freq[(ctx << 1) + 1] = (uint)(p * (1 << 10));
				p = 1.0 - p;
				pm->freq[(ctx << 1)] = (uint)(p * (1 << 10));
			}
			for (i = 0; i < pm->size; i++) {
				pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
			}
		}//fin level

//ヒストグラムをもとに確率モデルを作成することで最高効率を出す

		cpm->size = enc->num_class;
		cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
		cpm->cumfreq = &(cpm->freq[cpm->size]);
		cpm->cumfreq[0] = 0;
		for (i = 0; i < enc->num_class; i++) {
			cpm->freq[i] = hist[i];
			cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
		}
		for (j = 0; j < numidx; j++) {//4096.when 512*512
			i = index[j];//mtfの順番が保存されている
			if (i < 0) {//マイナスのものは調整
				i = -(i + 1);
				ctx = i & (~1);
				rc_encode(fp, enc->rc,
					pm->cumfreq[i] - pm->cumfreq[ctx], pm->freq[i],
					pm->cumfreq[ctx + 2] - pm->cumfreq[ctx],0);
			} else {
				rc_encode(fp, enc->rc, cpm->cumfreq[i], cpm->freq[i],
					cpm->cumfreq[cpm->size],0);
			}
		}
		bits += (int)enc->rc->code;
		enc->rc->code = 0;
	}
	free(index);
	free(hist);
	return (bits);
}

#if AUTO_PRD_ORDER

int encode_predictor(FILE *fp, ENCODER *enc, int flag)
{
	int cl, coef, sgn, k, m, min_m, bits, d;
	cost_t cost, min_cost, t_cost;
	PMODEL *pm;	
	uint cumb;
	int zrfreq, nzfreq;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	t_cost = 0.0;
	for (d = 0; d < enc->prd_mhd; d++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 8; m++) {
			cost = 0.0;
			for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
				for (cl = 0; cl < enc->num_class; cl++) {
					coef = enc->predictor[cl][k];
					if (coef < 0) coef = -coef;
					cost += enc->coef_cost[enc->zero_m[d]][m][coef];
				}
			}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		if (flag) enc->coef_m[d] = min_m;
		min_cost = INT_MAX;
		for (m = min_m = 0; m < NUM_ZMODEL; m++) {
			cost = 0.0;
			for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
				for (cl = 0; cl < enc->num_class; cl++) {
					coef = enc->predictor[cl][k];
					if (coef < 0) coef = -coef;
					cost += enc->coef_cost[m][enc->coef_m[d]][coef];
				}
			}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		if (flag) enc->zero_m[d] = min_m;
		t_cost += min_cost;
	}
	bits = (int)t_cost;
	/* Arithmetic */
	if (fp != NULL) {
		bits = 0;
		//		PMODEL *pm;
		//		uint cumb;
		//		int zrfreq, nzfreq;
		pm = &enc->spm;
		for (d = 0; d < enc->prd_mhd; d++) {
			rc_encode(fp, enc->rc, enc->zero_m[d], 1, NUM_ZMODEL,0);
			rc_encode(fp, enc->rc, enc->coef_m[d], 1, 8,0);
			nzfreq = enc->zero_fr[enc->zero_m[d]];
			zrfreq = TOT_ZEROFR - nzfreq;
			set_spmodel(pm, enc->max_coef + 1, enc->coef_m[d]);
			cumb = pm->freq[0];
			for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
				for (cl = 0; cl < enc->num_class; cl++) {
					coef = enc->predictor[cl][k];
					if (coef == 0) {
						rc_encode(fp, enc->rc, 0, zrfreq, TOT_ZEROFR,0);
					} else {
						rc_encode(fp, enc->rc, zrfreq, nzfreq, TOT_ZEROFR,0);
						sgn = (coef < 0)? 1 : 0;
						if (coef < 0) coef = -coef;
						rc_encode(fp, enc->rc, pm->cumfreq[coef] - cumb,  pm->freq[coef],
							pm->cumfreq[pm->size] - cumb,0);
						rc_encode(fp, enc->rc, sgn, 1, 2,0);
					}
				}
			}
		}
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}

#else

/* change pmodel each position */
int encode_predictor(FILE *fp, ENCODER *enc, int flag)
{
	int cl, coef, sgn, k, m, min_m, bits;
	cost_t cost, min_cost, t_cost;
	PMODEL *pm;
	pm = &enc->spm;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	t_cost = 0.0;
	for (k = 0; k < enc->prd_order; k++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 16; m++) {
			cost = 0.0;
			for (cl = 0; cl < enc->num_class; cl++) {
				coef = enc->predictor[cl][k];
				if (coef < 0) coef = -coef;
				cost += enc->coef_cost[m][coef];
			}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		t_cost += min_cost;
		if (flag) enc->coef_m[k] = min_m;
	}
	bits = (int)t_cost;
	if (fp != NULL) {
		bits = 0;
		for (k = 0; k < enc->prd_order; k++) {
			set_spmodel(pm, enc->max_coef + 1, enc->coef_m[k]);
			rc_encode(fp, enc->rc, enc->coef_m[k], 1, 16,0);
			for (cl = 0; cl < enc->num_class; cl++) {
				coef = enc->predictor[cl][k];
				sgn = (coef < 0)? 1 : 0;
				if (coef < 0) coef = -coef;
				rc_encode(fp, enc->rc, pm->cumfreq[coef],  pm->freq[coef],
					pm->cumfreq[pm->size],0);
				if (coef > 0) {
					rc_encode(fp, enc->rc, sgn, 1, 2,0);
				}
			}
		}
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}

#endif

int encode_threshold(FILE *fp, ENCODER *enc, int flag)
{
	int cl, gr, i, k, m, min_m, bits;
	cost_t cost, min_cost;
	PMODEL *pm;
	double p;
	pm = &enc->spm;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	/* Arithmetic */
	min_cost = INT_MAX;
	for (m = min_m = 0; m < 16; m++) {
		set_spmodel(pm, MAX_UPARA + 2, m);
		cost = 0.0;
		for (cl = 0; cl < enc->num_class; cl++) {
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				i = enc->th[cl][gr - 1] - k;
				p = (double)pm->freq[i]
				/ (pm->cumfreq[pm->size - k]);
				cost += -log(p);
				k += i;
				if (k > MAX_UPARA) break;
			}
		}
		cost /= log(2.0);
		if (cost < min_cost) {
			min_cost = cost;
			min_m = m;
		}
	}
	set_spmodel(pm, MAX_UPARA + 2, min_m);
	p = log(pm->cumfreq[MAX_UPARA + 2]);
	if (fp == NULL) {
		if (flag == 1){
			for (i = 0; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);
			}
		}
		bits = (int)min_cost;
	} else {
		rc_encode(fp, enc->rc, min_m, 1, 16,1);
		for (cl = 0; cl < enc->num_class; cl++) {
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				i = enc->th[cl][gr - 1] - k;
				rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i],
					pm->cumfreq[pm->size - k],1);
				k += i;
				if (k > MAX_UPARA) break;
			}
		}
		if (enc->num_pmodel > 1) {
			for (gr = 0; gr < enc->num_group; gr++) {
				pm = enc->pmlist[gr];
				rc_encode(fp, enc->rc, pm->id, 1, enc->num_pmodel,1);
			}
		}
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}



int encode_threshold_temp(FILE *fp, ENCODER *enc, int flag)
{
	int gr, i, k, bits;
	cost_t cost, min_cost;
	PMODEL *pm;
	double p;
	pm = &enc->spm;
  int count = 0;


#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	/* Arithmetic */
	min_cost = INT_MAX;
		set_spmodel(pm, MAX_UPARA + 2, 0);
		cost = 0.0;
			k = 0;
			for (gr = 0; gr < enc->num_group - 1; gr++) {
				i = enc->threshold[gr];
				p = (double)pm->freq[i]
				/ (pm->cumfreq[pm->size]);
				cost += -log(p);
				if (i > MAX_UPARA) break;
        count++;
			}
		cost /= log(2.0);
	set_spmodel(pm, MAX_UPARA + 2, 0);
	p = log(pm->cumfreq[MAX_UPARA + 2]);
	if (fp == NULL) {
		if (flag == 1){
			for (i = 0; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);
			}
		}
		bits = (int)min_cost;
	} else {
    printf("th_start_range = %llu\n",enc->rc->range);
		rc_encode(fp, enc->rc, count, 1, 16,1);
    printf("th_m_range = %llu\n",enc->rc->range);
    printf("count = %d\n",count);
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				i = enc->threshold[gr];
				if (i > MAX_UPARA) break;
        printf("th[%d] = %d\n",gr,i);
				rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i],
				pm->cumfreq[pm->size],1);
        printf("th_range = %llu\n",enc->rc->range);
        }
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
		}
	return (bits);
}

int encode_prd_threshold(FILE *fp, ENCODER *enc, int flag)
{
	int gr, i, k, bits;
	cost_t cost, min_cost;
	PMODEL *pm;
	double p;
	pm = &enc->spm;
  int count = 0;


#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	/* Arithmetic */
	min_cost = INT_MAX;
		set_spmodel(pm, MAX_UPARA + 2, 0);
		cost = 0.0;
			k = 0;
			for (gr = 0; gr < enc->num_group - 1; gr++) {
				i = enc->prd_threshold[gr];
				p = (double)pm->freq[i]
				/ (pm->cumfreq[pm->size]);
				cost += -log(p);
				if (i > MAX_UPARA) break;
        count++;
			}
		cost /= log(2.0);
	set_spmodel(pm, MAX_UPARA + 2, 0);
	p = log(pm->cumfreq[MAX_UPARA + 2]);
	if (fp == NULL) {
		if (flag == 1){
			for (i = 0; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);
			}
		}
		bits = (int)min_cost;
	} else {
    printf("th_start_range = %llu\n",enc->rc->range);
		rc_encode(fp, enc->rc, count, 1, 16,1);
    printf("th_m_range = %llu\n",enc->rc->range);
    printf("count = %d\n",count);
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				i = enc->prd_threshold[gr];
				if (i > MAX_UPARA) break;
        printf("prd_th[%d] = %d\n",gr,i);
				rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i],
				pm->cumfreq[pm->size],1);
        printf("th_range = %llu\n",enc->rc->range);
        }
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
		}
	return (bits);
}



struct DIFF
{
  int value;
  int id;
  
};


int select_nearest_var(int var){

//既存のσと実際の分散とが一番近い物を選ぶプログラム

int i,j;
struct DIFF diff[16],temp;


for(i = 0; i < 16 ; i++){

    diff[i].value = abs(var - sigma_a[i]);
    diff[i].id = i;
  }

 for (i = 0; i < 16 - 1; i++) {
        for (j = 16 - 1; j > i; j--) {
            if (diff[j - 1].value > diff[j].value) {  // 前の要素の方が大きかったら 
                temp = diff[j];        // 交換する 
                diff[j] = diff[j - 1];
                diff[j - 1]= temp;
            }
        } 
        
    }


return(diff[0].id);
}



double lngamma_w(double xx)//Γ関数
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




double calc_ggprob_w(double beta, double shape, double h, double x)
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

void set_freqtable_w(int *freq, int size, int ssize,
				   double shape, double sigma, double h, double off)
{
	double beta, norm;
	int i;


	/* Generalized Gaussian distribution */
	beta = exp(0.5*(lngamma_w(3.0/shape)-lngamma_w(1.0/shape))) / sigma;//一般化ガウス関数.ηのみ?.論文式(2-4).shape = cn
	norm = 0.0;//初期化
	for (i = 0; i < size; i++) {
		norm += calc_ggprob_w(beta, shape, h, i - off);//式(2-4)上式のexpとその中身
	}
	norm = (double)(MAX_TOTFREQ - size * MIN_FREQ) / norm;//MAX_TOTFREQ = (1 << 14) . MIN_FREQ = 1.
	norm += 1E-8;	/* to avoid machine dependent rounding errors */ //量子化するとき誤差が出ないようにするため
  for(i = 0 ; i < size;i++){	
	freq[i] = (uint)(norm * calc_ggprob_w(beta, shape, h, i - off ) + MIN_FREQ);//size分のfreq
	}
  return;
}

void w_newone_pmodels(int num_pmodel, int i,  double sigma, int size ,int *freq)
{
//	PMODEL ***pmodels, *pmbuf, *pm;
  int pm_accuracy = 3;
	int j = 0, num_subpm, ssize;
	double delta_c, c, s, sw, off;

  
	num_subpm = 1 << pm_accuracy;//num_subpm is 8.pm_accuracy is 3
	ssize = size;//size is 256
	size = size + ssize - 1;//size = 511
	sw = 1.0 / (double)num_subpm;// 1/8の量子化

	delta_c = 3.2 / (double)num_pmodel;//cnは0.2ずつ動くので. num_pmodel is 16 . delta_c = 0.2 
	off = (double)(ssize - 1);//off = 255
		s = sigma;//分散を代入
	  c = delta_c * (double)(i + 1);//常にこれを使用する　idとcnが配列の使用によってずれてしまうためidに1を足す
		set_freqtable_w(freq, size, ssize, c, s, sw/2.0, off - sw * j);//sは実際の分散値
	
  return;
}



double continuous_GGF(ENCODER *enc, double e,int w_gr){

int lngamma(double);
int cn,num_pmodel = enc->num_pmodel;
double sigma,delta_c,shape,eta,p;
double accuracy = 1 / (double)NAS_ACCURACY;
  cn = WEIGHT_CN; 
  sigma = enc->sigma[w_gr];

	delta_c = 3.2 / (double)num_pmodel;//cnは0.2ずつ動くので. num_pmodel is 16 . delta_c = 0.2 
	shape = delta_c * (double)(cn);//cnを実際の数値に変換

	eta = exp(0.5*(lgamma(3.0/shape)-lgamma(1.0/shape))) / sigma;//一般化ガウス関数.ηのみ

if(e <= accuracy){
  p = 1;
}else{
	p = exp(-pow(eta * (e), shape));
}
  
return(p);
}

void composite_MMF(ENCODER *enc, PMODEL *pm,int x,int y,int u,
 uint *freq_shift_array,uint *cumfreq_shift_array,int limit){

int i,prd;
int start = enc->coef_start,end = enc->coef_end;
double coef = 0;
uint *freq_array,*array_original;

  if(u > limit){

  if( (u < start) || (end < u) ){
    if( u < start ){return;;
    }else if( end < u ){coef = enc->prd_coef[end - start - 1];
    }else{printf("error - accidental U\n");}
  }else{
  coef = enc->prd_coef[u - start];
  }
  if(coef <= 0){return;}

freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));
array_original =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));

  for(i = 0;i < pm->size ; i++){
    array_original[i]  = freq_shift_array[i];
  }
  for(i = 0;i < pm->size ; i++){
    freq_array[i]  = 1;
  }
  prd = enc->prd_b[y][x];
  prd = CLIP(0, enc->maxprd, prd);
  prd = prd >> 6;
  for(i = 0;i < pm->size ; i++){
    freq_array[i + prd]  = pm->freq[i];
  }

  int buf1,buf2;  
  double buf;
  buf1 = 0;
  buf2 = 0;
  for(i = 0;i < pm->size ; i++){
    if(buf1 < array_original[i]){buf1 = array_original[i];}
    if(buf2 < freq_array[i]){buf2 = freq_array[i];}
  }
  if(buf1 > buf2){
  buf = ((double)buf1 / (double)buf2);
  }else{
  buf = ((double)buf2 / (double)buf1);
  }
  for(i = 0;i < pm->size ; i++){
    freq_array[i] = (int)((double)freq_array[i] / buf * coef); 
    array_original[i]  = array_original[i];
  }
  for(i = 0;i < pm->size ; i++){
    freq_array[i] = (freq_array[i] + array_original[i]) ; 
    if(freq_array[i] == 0){
      freq_array[i] = 1;
    }
  }
  for(i = 0;i < pm->size ; i++){
    freq_shift_array[i] = freq_array[i]; 
  }
  cumfreq_shift_array[0] = 0;

	for (i = 0; i < pm->size; i++) {
		cumfreq_shift_array[i + 1] = cumfreq_shift_array[i] + freq_shift_array[i];
	} 
  
free(freq_array);
free(array_original);
  }
return;
}


int Shift_freq(PMODEL *pm, int x , int y, ENCODER *enc,int ***a
,uint **freq_shift_array_save,uint *freq_shift_array,uint *cumfreq_shift_array,
unsigned long long int *sum_freq_shift_array,int multimodal,int w_gr,int w_func,int flg,int flgi){  //稲村flgi追加

  int i,j,by[multimodal],bx[multimodal],e[multimodal];
  int flag2 = 0;
  double er_sum = 0.0,er_coef = 0.0;
  int w[multimodal];
//  int multi = multimodal;

int multi = multimodal - 1; //稲村変更



  
  int mc[multimodal];
  double mc_w[multimodal];
  int flg_min = 0;
  

//fundamental用
int *org_p,*roff_p;
int area1_sum;
double area1_ave;

double mc_double;

if(flg == 0){
  for(i = 0 ; i < multi ; i++){//いくつの山にするか
    by[i] = a[y][x][i*4+1];
    bx[i] = a[y][x][i*4+2];  
    e[i] = 0;
    for(j = 0;j < pm->size ; j++){
       freq_shift_array[j+e[i]]  = pm->freq[j];
    }
  }
}

else{

  if( (x < 5) && (y == 0) ){
    multi = x;
    flag2 = 1;
    flg_min = 1;
  }

if((flg_min != 1) && (enc->temp_num[y][x] < multi) ){
multi = enc->temp_num[y][x];
flag2 = 0;
}

for(j = 0;j < pm->size  ; j++){//初期化
  for(i = 0;i < multi ; i++){
    freq_shift_array_save[i][j] = 1;
  }
  freq_shift_array[j] = 1;
  sum_freq_shift_array[j] = 0;
}

	  roff_p = enc->roff[y][x];
		org_p = &enc->org[y][x];
    area1_sum = 0;
    area1_ave = 0;
		for(i=0;i < AREA; i++){//市街地距離AREA個分
			area1_sum += org_p[roff_p[i]];
	  }
    area1_ave = area1_sum / AREA;

for(i = 0 ; i < multi ; i++){//いくつの山にするか
    by[i] = a[y][x][i*4+1];
    bx[i] = a[y][x][i*4+2];  

    e[i] = (int)(enc->encval[by[i]][bx[i]] - enc->array[y][x][i] + area1_ave);//TM最小の時の輝度値
    
    if(e[i] < 0 || e[i] > enc->maxval){e[i] = area1_ave;}
    mc[i] = a[y][x][i*4+3];


 for(j = 0;j < pm->size ; j++){
    freq_shift_array_save[i][j+e[i]]  = pm->freq[j]; //freq_arrayに値を代入していく.
  }
}//multi fin











if(flgi == 1){  //稲村


if (y == 0 || x == 0){


 by[multi] = a[y][x][multi*4+1];
    bx[multi] = a[y][x][multi*4+2];  

    e[multi] = (int)(enc->encval[by[multi]][bx[multi]] - enc->array[y][x][multi] + area1_ave);//TM最小の時の輝度値
    
    if(e[multi] < 0 || e[multi] > enc->maxval){e[multi] = area1_ave;}
    mc[multi] = a[y][x][multi*4+3];


 for(j = 0;j < pm->size ; j++){
    freq_shift_array_save[multi][j+e[multi]]  = pm->freq[j]; //freq_arrayに値を代入していく.

}




multi++;


}

else{


    by[multi] = a[y][x][multi*4+1];
    bx[multi] = a[y][x][multi*4+2];  

    e[multi] = (int)((enc->org[y - 1][x] + enc->org[y][x - 1] +  enc->org[y - 1][x - 1]) / 3);//稲村変更


//printf("%d,%d,%lf\n", enc->encval[by[multi]][bx[multi]],enc->array[y][x][multi],area1_ave);


    
    if(e[multi] < 0 || e[multi] > enc->maxval){e[multi] = area1_ave;}
    mc[multi] = a[y][x][multi*4+3];//稲村変更


 for(j = 0;j < pm->size ; j++){
    freq_shift_array_save[multi][j+e[multi]]  = pm->freq[j]; //freq_arrayに値を代入していく.

}

multi++;   //稲村追加

}

}














///////////////////////////
//重み付け，正規化処理/////
///////////////////////////

if(flag2 == 1){

for(i = 0 ; i < multi ; i++){
  for(j = 0;j < pm->size ; j++){
    sum_freq_shift_array[j] += freq_shift_array_save[i][j];//変形させたものをsumに入れていく
  }
}

for(j = 0;j < pm->size ; j++){
  sum_freq_shift_array[j] = sum_freq_shift_array[j] / multi;
}

}else{

/////////////////////////////////////
///一般化ガウス関数による重み付け////
/////////////////////////////////////


if(w_func == 1){

for(i = 0 ; i < multi ; i++){
mc_double = 0;
if(mc[i] != 0){
mc_double = (double)mc[i] / NAS_ACCURACY; 
}
mc_w[i] = continuous_GGF(enc, mc_double,w_gr);
er_sum = er_sum + mc_w[i];
}

if(er_sum == 0){er_sum = 1;}

er_coef = (100) / er_sum;//er_coefが0になるのを未然に防ぐ
for(i = 0 ; i < multi ; i++){
  w[i] = (int)(mc_w[i] * er_coef);
}

/////////////////////////
///逆数による重み付け////
/////////////////////////

}else if(w_func == 0){

  for(i = 0 ; i < multi ; i++){//いくつの山にするか
    if(mc[i] == 0){mc[i] = 1;}
    mc_w[i] = 1 / ((double)mc[i] / NAS_ACCURACY);

    er_sum = er_sum + mc_w[i];
  }//multi fin

er_coef = (100) / er_sum;//er_coefが0になるのを未然に防ぐ
  for(i = 0 ; i < multi ; i++){
    if(mc_w[i] * er_coef < 1.0){
    w[i] = 0;
    }else{
    w[i] = (int)(mc_w[i] * er_coef);
    }
  }
}

for(i = 0 ; i < multi ; i++){
  for(j = 0;j < pm->size ; j++){
     freq_shift_array_save[i][j] = freq_shift_array_save[i][j] * w[i];//全てに重みをかけて保存
  }
  
  for(j = 0;j < pm->size ; j++){
    sum_freq_shift_array[j] += freq_shift_array_save[i][j];//重みをつけたものをsumに入れていく
  }
}//multimodal fin

for(j = 0;j < pm->size ; j++){
    sum_freq_shift_array[j] = sum_freq_shift_array[j] / 100;
    if(sum_freq_shift_array[j] == 0){
      sum_freq_shift_array[j] = 1;
    }
   }
}//else fin

for(j = 0;j < pm->size  ; j++){
  freq_shift_array[j] = sum_freq_shift_array[j];
  }
}//else fin
    cumfreq_shift_array[0] = 0;

	for (i = 0; i < pm->size; i++) {
		cumfreq_shift_array[i + 1] = cumfreq_shift_array[i] + freq_shift_array[i];
	}




return(0);
}

void calc_prd(ENCODER *enc)
{
int *org_p,*roff_p,**prd_m;
int k,x,y,prd = 0,prd_buf;
int prd_max = enc->maxval << 6;

prd_m = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));

  for(y = 0; y < enc->height;y++){
    for(x = 0; x < enc->width;x++){
      roff_p = enc->roff[y][x];
      org_p = &enc->org[y][x];
      prd_buf = prd;
      prd = 0;
      for (k = 0; k < PRED_AREA; k++) {
         prd += org_p[roff_p[k]] * enc->temp_predictor[y][x][k];
      }

			//prd = CLIP(0, enc->maxprd, prd);
      if(prd < 0 || prd_max < prd ){prd = prd_buf;}
      prd_m[y][x] = prd;
    }
  }
enc->prd_b = prd_m;

return;
}


struct point{
  int x;
  int y;
};

struct point get_ref_pels(int x, int y, int num){
//市街地距離３の範囲の輝度を返す(当該の左から，時計回り)
struct point ref[TRAIN_AREA] = {
  {0,0},
  {-1,0},{0,-1},
  {-2,0},{-1,-1},{0,-2},{1,-1},
  {-3,0},{-2,-1},{-1,-2},{0,-3},{1,-2},{2,-1},
  {-4,0},{-3,-1},{-2,-2},{-1,-3},{0,-4},{1,-3},{2,-2},{3,-1},
  {-5,0},{-4,-1},{-3,-2},{-2,-3},{-1,-4},{0,-5},{1,-4},{2,-3},{3,-2},{4,-1},
  {-6,0},{-5,-1},{-4,-2},{-3,-3},{-2,-4},{-1,-5},{0,-6},{1,-5},{2,-4},{3,-3},{4,-2},{5,-1},
  {-7,0},{-6,-1},{-5,-2},{-4,-3},{-3,-4},{-2,-5},{-1,-6},{0,-7},{1,-6},{2,-5},{3,-4},{4,-3},{5,-2},{6,-1},
  {-8,0},{-7,-1},{-6,-2},{-5,-3},{-4,-4},{-3,-5},{-2,-6},{-1,-7},{0,-8},{1,-7},{2,-6},{3,-5},{4,-4,},{5,-3},{6,-2},{7,-1},
  {-9,0},{-8,-1},{-7,-2},{-6,-3},{-5,-4},{-4,-5},{-3,-6},{-2,-7},{-1,-8},{0,-9},{1,-8},{2,-7},{3,-6},{4,-5},{5,-4},{6,-3},{7,-2},{8,-1},
  {-10,0},{-9,-1},{-8,-2},{-7,-3},{-6,-4},{-5,-5},{-4,-6},{-3,-7},{-2,-8},{-1,-9},{0,-10},{1,-9},{2,-8},{3,-7},{4,-6},{5,-5},{6,-4},{7,-3},{8,-2},{9,-1}
};
/*
x += ref[num].x;
y += ref[num].y;
*/
return(ref[num]);
}



 void design_predictor_for_temp(ENCODER *enc, int f_mmse)
{
  double **mat, *weight, w, e, d, pivot;
  int x, y, i, j, k,u, gr, pivpos, *index, *roff_p, *org_p;
  int ***predictor;
  int *th_p;
  struct point pos;
  int num,bx,by;

  predictor = (int***)alloc_3d_array(enc->height , enc->width ,PRED_AREA ,sizeof(int));

  mat = (double **)alloc_2d_array(PRED_AREA, PRED_AREA + 1, sizeof(double));
  index = (int *)alloc_mem(sizeof(int) * PRED_AREA);
  weight = (double *)alloc_mem(sizeof(double) * enc->num_group);

  for (gr = 0; gr < enc->num_group; gr++) {//gr は分散の番号16種類
    if (f_mmse) {//f_mmseは基本0
      weight[gr] = 1.0;
    } else {
      weight[gr] = 1.0 / (enc->sigma[gr] * enc->sigma[gr]);//当該分散番号のσ^2で割っている
    }
  }
    for (i = 0; i < PRED_AREA; i++) {
      for (j = 0; j <= PRED_AREA; j++) {
        mat[i][j] = 0.0;//初期化
      }
    }

//    th_p = enc->threshold;
    for (y = 0; y < enc->height; y++) {
      for (x = 0; x < enc->width; x++) {
        if(f_mmse){
          gr = 5;//二乗誤差最小では関係ないので適当
        }else{

          th_p = enc->threshold;
          u = enc->TM_U[y][x];      
          for (gr = 0; gr < enc->num_group; gr++) {
            if (u < th_p[gr]) break;
          }
        }
        for (i = 0; i < PRED_AREA; i++) {
          for (j = 0; j <= PRED_AREA; j++) {
            mat[i][j] = 0.0;
          }
        }

        for(num = 1; num < TRAIN_AREA; num++){//13画素分　num:0〜12
          pos = get_ref_pels( x, y, num);
          bx = x + pos.x;
          by = y + pos.y;
          if((by < 0) || (by > enc->height))continue;
          if((bx < 0) || (bx > enc->width))continue;

          roff_p = enc->roff[by][bx];
          org_p = &enc->org[by][bx];
          for (i = 0; i < PRED_AREA; i++) {
            w = weight[gr] * org_p[roff_p[i]];
              for (j = i; j < PRED_AREA; j++) {
                mat[i][j] += w * org_p[roff_p[j]];
              }
            mat[i][PRED_AREA] += w * org_p[0];//当該の答え
          }
        }

        for (i = 0; i < PRED_AREA; i++) {
          index[i] = i;//予測器とインデックスの関連付け
          for (j = 0; j < i; j++) {
            mat[i][j] = mat[j][i];//対称行列にする
          }
        }
        for (i = 0; i < PRED_AREA; i++) {
          pivpos = i;
          pivot = fabs(mat[index[i]][i]);
          for (k = i + 1; k < PRED_AREA; k++) {
            if (fabs(mat[index[k]][i]) > pivot) {
              pivot = fabs(mat[index[k]][i]);
              pivpos = k;
            }
          }
          k = index[i];
          index[i] = index[pivpos];
          index[pivpos] = k;
          if (pivot > 1E-10) {
            d = mat[index[i]][i];
            for (j = i; j <= PRED_AREA; j++) {
              mat[index[i]][j] /= d;
            }
            for (k = 0; k < PRED_AREA; k++) {
              if (k == i) continue;
              d = mat[index[k]][i];
              for (j = i; j <= PRED_AREA; j++) {
                mat[index[k]][j] -= d * mat[index[i]][j];
              }
            }
          }
        }
        w = (1 << enc->coef_precision);
        e = 0.0;
        for (i = 0; i < PRED_AREA; i++) {
          if (fabs(mat[index[i]][i]) > 1E-10) {
            d = mat[index[i]][PRED_AREA] * w;
          } else {
            d = 0.0;
          }
          k = (int)d;
          if (k > d) k--;
          if (k < -enc->max_coef) {
            d = k = -enc->max_coef;
          } else if (k > enc->max_coef) {
            d = k = enc->max_coef;
          }
          predictor[y][x][i] = k;
          d -= k;
          e += d;
          mat[index[i]][PRED_AREA] = d;
        }
      }//x fin
    }//y fin

  enc->temp_predictor = predictor;


  calc_prd(enc);


  free(weight);
  free(index);
  free(mat);

}


//roff出力値をx,yに変換するプログラム。使えるかどうかは知らない

void convert_roff_value(ENCODER *enc,int r ,int *roff_x,int *roff_y){
int roff;
int quotient,scs,k;
int width = enc->width;

roff = r;
quotient = roff / width;
scs = roff % width;
k = scs + width;
if((width / 2 < k) && (width >= k)){k = k - width;}
*roff_y = quotient;
*roff_x = k;

return;
}





struct EI_Member
{
  double cost;
  int gr;
  int cn;
  
};

struct EI_Member_M
{
  double cost;
  int gr;
  int cn;
  int multi;
};

void make_cost_u(ENCODER *enc , int ***array){

  printf("Making U...\n");

	int x, y, e,k,c, base, bits, gr,cn, cumbase,flg;
  double cost,subcost,tot_cost,a;
  int multi=0;
	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;

  int **temp_u;
  double area_cost;
  int **area_cost_qt;
  int *roff_p;
  int *org_p;
  double *wt_p;
  int count;
  int u;
  area_cost_qt = enc->TM_U;
  temp_u =  (int **)alloc_2d_array(enc->height ,enc->width , sizeof(int));

freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));



int **cost_save;
  cost_save = (int **)alloc_2d_array(  enc->height+1 , enc->width ,sizeof(int));
	cost_save[enc->height][0] = 0;
double a_cost;

int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

gr = enc->gr = 10;
cn = enc->cn;
multi = enc->multi;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
			pm = enc->pmodels[gr][cn];

for(y = 0; y < enc->height ; y++){
		for (x = 0; x < enc->width; x++) {
    cost_save[y][x] = 0;
    }
}

	wt_p = enc->ctx_weight_double;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      if(y == 0 && x == 0){
      u = 0;
      }else{
	

				roff_p = enc->roff[y][x];
        org_p = &cost_save[y][x];
        area_cost = 0;        
        count = 0;
        for(k=0;k < U_AREA; k++){

          a_cost = (double)org_p[roff_p[k]];
          a_cost = a_cost / COST_ACCURACY;
          if(a_cost == 0){
          a_cost = 0.1 * (wt_p[k]);
          }else{
          a_cost = a_cost * (wt_p[k]) ;
          }

 				  area_cost += a_cost; 

          if(a_cost != 0){
          count++;
          }
        }

        area_cost_qt[y][x] = (int)(area_cost / count );
        if(area_cost_qt[y][x] > MAX_UPARA){
          area_cost_qt[y][x] = MAX_UPARA;
        }
      }//else fin      


        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入
      base = 255;
      
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      cost_save[y][x] = (int)((cost + subcost) * COST_ACCURACY);//Uのためにコストを保存
      tot_cost += cost + subcost;
		
		}//x fin
	}//y fin


//free(cost_save);
return;
}



void make_cost_u_2nd(ENCODER *enc , int ***array){


  printf("Making U...\n");
	int x, y, e,k,c, base, bits, gr,cn, cumbase,flg;
  double cost,subcost,tot_cost,a;
  int multi=0;
	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;

  double area_cost;
  int **area_cost_qt;
  int *roff_p;
  int *org_p;
  double *wt_p;
  int count,u;
  int *th_p;
  area_cost_qt = enc->TM_U;

freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

cn = enc->cn;
multi = enc->multi;
int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

double a_cost;
int **cost_save;
  cost_save = (int **)alloc_2d_array(  enc->height+1 , enc->width ,sizeof(int));
	cost_save[enc->height][0] = 0;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
for(y = 0; y < enc->height ; y++){
		for (x = 0; x < enc->width; x++) {
    cost_save[y][x] = 0;
    }
}
	wt_p = enc->ctx_weight_double;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

      th_p = enc->threshold;
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      //cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入
      base = 255;
      
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      cost_save[y][x] = (int)((cost + subcost) * COST_ACCURACY);//Uのためにコストを保存
      tot_cost += cost + subcost;
		
      if(y == 0 && x == 0){
      u = 0;
      }else{
	

				roff_p = enc->roff[y][x];
        org_p = &cost_save[y][x];
        area_cost = 0;        
        count = 0;
        for(k=0;k < U_AREA; k++){

          a_cost = (double)org_p[roff_p[k]];
          a_cost = a_cost / COST_ACCURACY;
          if(a_cost == 0){
          a_cost = 0.1 * (wt_p[k]);
          }else{
          a_cost = a_cost * (wt_p[k]) ;
          }

 				  area_cost += a_cost; 

          if(a_cost != 0){
          count++;
          }
        }

        area_cost_qt[y][x] = (int)(area_cost / count );
        if(area_cost_qt[y][x] > MAX_UPARA){
          area_cost_qt[y][x] = MAX_UPARA;
        }
      }      
		}//x fin
	}//y fin

return;
}

void make_prd_th(ENCODER *enc ){

  cost_t min_cost, *dpcost, *cbuf_p, *thc_p;
	int th1, th0;;
	int **trellis;

	int x, y, e,k,u, base, bits, gr,cn;
  double cost,tot_cost,a;
  int *threshold;
	PMODEL *pm;
  int prd; 
  cost_t **cbuf,**tempbuf;

threshold = (int *)alloc_mem((enc->num_group) * sizeof(int));
cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));
tempbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));

cn = enc->cn;

printf("Making th\n");


bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
		thc_p = enc->th_cost;

	for (k = 0; k < MAX_UPARA + 2; k++){ trellis[0][k] = 0;}
	/* Dynamic programming */
		for (gr = 0; gr < enc->num_group; gr++) {
			for (u = 0; u < MAX_UPARA + 2; u++) {
			  cbuf[gr][u] = 0;
			}
		}

	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

      u = enc->TM_U[y][x];
			prd = enc->prd_b[y][x];
			prd = CLIP(0, enc->maxprd, prd);
			base = enc->bconv[prd];
			e = enc->encval[y][x];//輝度値代入

      for(gr = 0 ; gr < enc->num_group ;gr++){
			  pm = enc->pmodels[gr][cn];
			  cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
      }

		}//x fin
	}//y fin

		for (gr = 0; gr < enc->num_group; gr++) {
      cbuf[gr][0] = 0;
			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf[gr][u] += cbuf[gr][u - 1];
			}
		}

		cbuf_p = cbuf[0];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u] + thc_p[u];
		}
		for (gr = 1; gr < enc->num_group - 1; gr++) {
			cbuf_p = cbuf[gr];
			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
        
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}
				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;

			}
		}

		cbuf_p = cbuf[gr];
		th1 = MAX_UPARA + 1;
		th0 = th1;
		min_cost = dpcost[th1] - cbuf_p[th1];
		for (k = 0; k < th1; k++) {
			cost = dpcost[k] - cbuf_p[k];
			if (cost < min_cost) {
				min_cost = cost;
				th0 = k;
			}
		}
		trellis[gr][th1] = th0;


		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
      threshold[gr - 1] = th1;
		}
    threshold[enc->num_group - 1] = MAX_UPARA + 1;
    threshold[0] = 0;
    for(gr = 0;gr < enc->num_group;gr++){printf("prd_threshold=%d\n",threshold[gr]);}

      enc->prd_threshold = threshold;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < threshold[gr]) break;
      }
      enc->PRD_gr[y][x] = gr;
    }
  }

return;
}


void make_prd_th_2nd(ENCODER *enc , int ***array){

  cost_t min_cost, *dpcost, *cbuf_p, *thc_p;
	int th1, th0;;
	int **trellis;

	int x, y, e,c,k,u, base, bits, gr,cn, cumbase,flg,*th_p;
  int multi=0;
  double cost,subcost,tot_cost,a;
  int *threshold;
	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
  cost_t **cbuf,**tempbuf;

threshold = (int *)alloc_mem((enc->num_group) * sizeof(int));
cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));
tempbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));

freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

  uint *freq_array_buf,*cumfreq_array_buf;
freq_array_buf =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));
cumfreq_array_buf =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));


cn = enc->cn;
multi = enc->multi;
int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Making th\n");


enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
		thc_p = enc->th_cost;
	for (k = 0; k < MAX_UPARA + 2; k++){ trellis[0][k] = 0;}
	/* Dynamic programming */
		for (gr = 0; gr < enc->num_group; gr++) {
			for (u = 0; u < MAX_UPARA + 2; u++) {
			  cbuf[gr][u] = 0;
			}
		}

	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

      u = enc->TM_U[y][x];

      th_p = enc->threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

      for(k = 0 ; k < pm->size + 1; k++){
        freq_array_buf[k] = freq_array[k];
        cumfreq_array_buf[k] = cumfreq_array[k];
      }


      for(gr = 0 ; gr < enc->num_group ;gr++){

      cn = enc->PRD_cn;
			pm = enc->pmodels[gr][cn];

      composite_MMF(enc,pm, x, y,u,freq_array,cumfreq_array,enc->limit_PRD);




			e = enc->encval[y][x];//輝度値代入
      base = 255;
      
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      cbuf[gr][u] += cost + subcost;//しきい値のために分散ごとのコストを保存

      for(k = 0 ; k < 512 + 2; k++){
        freq_array[k] = freq_array_buf[k];
        cumfreq_array[k] = cumfreq_array_buf[k];
      }


      }

		}//x fin
	}//y fin

		for (gr = 0; gr < enc->num_group; gr++) {
      cbuf[gr][0] = 0;
			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf[gr][u] += cbuf[gr][u - 1];
			}
		}

		cbuf_p = cbuf[0];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u] + thc_p[u];
		}
		for (gr = 1; gr < enc->num_group - 1; gr++) {
			cbuf_p = cbuf[gr];
			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
        
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}
				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;

			}
		}

		cbuf_p = cbuf[gr];
		/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
		th1 = MAX_UPARA + 1;
		th0 = th1;
		min_cost = dpcost[th1] - cbuf_p[th1];
		for (k = 0; k < th1; k++) {
			cost = dpcost[k] - cbuf_p[k];
			if (cost < min_cost) {
				min_cost = cost;
				th0 = k;
			}
		}
		trellis[gr][th1] = th0;


		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
      threshold[gr - 1] = th1;
		}
    threshold[enc->num_group - 1] = MAX_UPARA + 1;
    threshold[0] = 0;
    for(gr = 0;gr < enc->num_group;gr++){printf("prd_threshold=%d\n",threshold[gr]);}
      enc->prd_threshold = threshold;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < threshold[gr]) break;
      }
      enc->PRD_gr[y][x] = gr;
    }
  }

return;
}




void make_th(ENCODER *enc , int ***array){

  cost_t min_cost, *dpcost, *cbuf_p, *thc_p;
	int th1, th0;;
	int **trellis;

	int x, y, e,c,k,u, base, bits, gr,cn, cumbase,flg;
  int multi=0;
  double cost,subcost,tot_cost,a;
  int *threshold;
	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
  cost_t **cbuf,**tempbuf;

threshold = (int *)alloc_mem((enc->num_group) * sizeof(int));
cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));
tempbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(cost_t));

freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

cn = enc->cn;
multi = enc->multi;
int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Making th\n");


enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
		thc_p = enc->th_cost;
	for (k = 0; k < MAX_UPARA + 2; k++){ trellis[0][k] = 0;}
	/* Dynamic programming */
		//if (enc->cl_hist[cl] == 0) continue;
		for (gr = 0; gr < enc->num_group; gr++) {
			for (u = 0; u < MAX_UPARA + 2; u++) {
			  cbuf[gr][u] = 0;
			}
		}

	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

      u = enc->TM_U[y][x];

      for(gr = 0 ; gr < enc->num_group ;gr++){
      //cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入
      base = 255;
      
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      cbuf[gr][u] += cost + subcost;//しきい値のために分散ごとのコストを保存

      }

		}//x fin
	}//y fin

		for (gr = 0; gr < enc->num_group; gr++) {
      cbuf[gr][0] = 0;
			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf[gr][u] += cbuf[gr][u - 1];
			}
		}

		cbuf_p = cbuf[0];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u] + thc_p[u];
		}
		for (gr = 1; gr < enc->num_group - 1; gr++) {
			cbuf_p = cbuf[gr];
			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
        
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}
				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;

			}
		}

		cbuf_p = cbuf[gr];
		/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
		th1 = MAX_UPARA + 1;
		th0 = th1;
		min_cost = dpcost[th1] - cbuf_p[th1];
		for (k = 0; k < th1; k++) {
			cost = dpcost[k] - cbuf_p[k];
			if (cost < min_cost) {
				min_cost = cost;
				th0 = k;
			}
		}
		trellis[gr][th1] = th0;


		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
      threshold[gr - 1] = th1;
      printf("%d,%d\n",gr - 1,threshold[gr - 1]); 
			//enc->th[cl][gr - 1] = th1;
		}
    threshold[enc->num_group - 1] = MAX_UPARA + 1;
    threshold[0] = 0;
    for(gr = 0;gr < enc->num_group;gr++){printf("threshold=%d\n",threshold[gr]);}
      enc->threshold = threshold;


	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < threshold[gr]) break;
      }
      enc->TM_gr[y][x] = gr;
    }
  }


return;
}


void w_gr_search(ENCODER *enc , int ***array, int flag){

	int x, y, e,i,l,j,k,c, base, bits, gr,cn, cumbase,flg;
  int multi,th = 0;
  double cost,subcost,tot_cost,a;
  struct EI_Member ei[enc->num_group ];
  struct EI_Member temp_ei;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

int u,*th_p;
int w_gr;

int w_func = 1;

printf("Set W_gr...\n");
setbuf(stderr,NULL);

	bits = 0;
  bzero(&ei, sizeof(ei));
  bzero(&temp_ei, sizeof(temp_ei));

multi = enc->multi;
cn = enc->cn;
th = enc->tm_th;

  k = 0;
  l = 0;

  for(w_gr = 0 ; w_gr < enc->num_group ;w_gr++){

bits = 0;
tot_cost = 0;
      th_p = enc->threshold;

a = 1.0 / log(2.0);

	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      //cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入

      base = 255;
      
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;

		}//x fin
	}//y fin

ei[k].cost = tot_cost;
ei[k].gr = w_gr;
ei[k].cn = cn;
k++;

fprintf(stderr,"%3.1d / 16\r",l+1);
l++;

}//w_gr fin
fprintf(stderr,"...OK      \r");
printf("\n");


 for (i = 0; i <(  enc->num_group  ) - 1; i++) {
        for (j = ( enc->num_group  ) - 1; j > i; j--) {
            if (ei[j - 1].cost > ei[j].cost) {  // 前の要素の方が大きかったら 
                temp_ei = ei[j];        // 交換する 
                ei[j] = ei[j - 1];
                ei[j - 1]= temp_ei;
            }
        } 
    }




printf("w_gr = %d  \n",ei[0].gr);

enc->tm_cost = ei[0].cost;
enc->w_gr = ei[0].gr;



}

void set_cn_old(ENCODER *enc , int ***array){

	int x, y, e,i,l,j,k,c, base, bits,cn,gr = 0, cumbase,flg = 0;
  int multi;
  double cost,subcost,tot_cost,a;
  struct EI_Member ei[enc->num_pmodel ];
  struct EI_Member temp_ei;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

int w_gr;
w_gr = enc->w_gr;
int w_func = 1;
int u,*th_p;

printf("Set Cn...\n");
setbuf(stderr,NULL);

	bits = 0;

  bzero(&ei, sizeof(ei));
  bzero(&temp_ei, sizeof(temp_ei));

multi = enc->multi;

  k = 0;
  l = 0;
      
      th_p = enc->threshold;
    for(cn = 0 ; cn < enc->num_pmodel ;cn++){


bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入
      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;
		
		}//x fin
	}//y fin

ei[k].cost = tot_cost;
ei[k].gr = gr;
ei[k].cn = cn;
//printf("cost[%d] = %f\n",k,tot_cost);
k++;


fprintf(stderr,"%3.1d / 16\r",l+1);
l++;

}//cn fin
fprintf(stderr,"...OK      \r");
printf("\n");

 for (i = 0; i <(  enc->num_pmodel  ) - 1; i++) {
        for (j = (  enc->num_pmodel  ) - 1; j > i; j--) {
            if (ei[j - 1].cost > ei[j].cost) {  // 前の要素の方が大きかったら 
                temp_ei = ei[j];        // 交換する 
                ei[j] = ei[j - 1];
                ei[j - 1]= temp_ei;
            }
        } 
    }

printf("| cn = %d | multi = %d |\n",ei[0].cn,enc->multi);

enc->tm_cost = ei[0].cost;
enc->PRD_cn = cn = enc->cn = ei[0].cn;

}

void set_cn_prd(ENCODER *enc , int ***array){

	int x, y, e,i,l,j,k,c, base, bits,cn,gr=0, cumbase;
  int multi;
  double cost,subcost,tot_cost,a;
  struct EI_Member ei[enc->num_pmodel ];
  struct EI_Member temp_ei;

  int prd_buf,prd = 0;
  int *roff_p;
  int *org_p;
  int prd_max = enc->maxval << 6;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));


int w_gr;
w_gr = enc->w_gr;

printf("Set Cn...\n");
setbuf(stderr,NULL);

	bits = 0;

  bzero(&ei, sizeof(ei));
  bzero(&temp_ei, sizeof(temp_ei));

multi = enc->multi;

  k = 0;
  l = 0;
      
    for(cn = 0 ; cn < enc->num_pmodel ;cn++){


  k = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      gr = enc->PRD_gr[y][x];
 			pm = enc->pmodels[gr][cn];
      prd_buf = prd;
      prd = 0;
    	roff_p = enc->roff[y][x];
      org_p = &enc->org[y][x];
      prd = enc->prd_b[y][x];
      if(prd < 0 || prd_max < prd ){prd = prd_buf;}

			base = enc->bconv[prd];

			e = enc->encval[y][x];//輝度値代入
			cumbase = pm->cumfreq[base];//baseまでの高さの和

      cost = (float)(-a * log(pm->freq[base + e]));
      c = pm->cumfreq[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;
//      }
		}//x fin
	}//y fin

ei[cn].cost = tot_cost;
ei[cn].cn = cn;
ei[cn].gr = gr;
printf("gr= %d cn= %d cost= %f\n",gr,cn,tot_cost);
  }//cn

 for (i = 0; i <(  enc->num_pmodel  ) - 1; i++) {
        for (j = (  enc->num_pmodel  ) - 1; j > i; j--) {
            if (ei[j - 1].cost > ei[j].cost) {  // 前の要素の方が大きかったら 
                temp_ei = ei[j];        // 交換する 
                ei[j] = ei[j - 1];
                ei[j - 1]= temp_ei;
            }
        } 
    }
for(l = 0 ; l < enc->num_pmodel ;l++){
printf("|cost= %f,cn= %d|\n",ei[l].cost,ei[l].cn);
}


enc->PRD_cn = ei[0].cn;
printf("|TM_cn[%d]= %d|\n",gr,ei[0].cn);


fprintf(stderr,"...OK      \r");
printf("\n");



}

void set_cn_prd_second(ENCODER *enc , int ***array){

	int x, y, e,i,l,j,k,c, base, bits,cn,gr=0, cumbase,flg;
  int multi;
  double cost,subcost,tot_cost,a;
  struct EI_Member ei[enc->num_pmodel ];
  struct EI_Member temp_ei;


	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));


int w_gr;
w_gr = enc->w_gr;
int w_func = 1;
int u,*th_p;

printf("Set Cn...\n");
setbuf(stderr,NULL);

	bits = 0;

  bzero(&ei, sizeof(ei));
  bzero(&temp_ei, sizeof(temp_ei));

multi = enc->multi;

  k = 0;
  l = 0;
      
    for(cn = 0 ; cn < enc->num_pmodel ;cn++){


  k = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      u = enc->TM_U[y][x];
      gr = enc->TM_gr[y][x];
 			pm = enc->pmodels[gr][enc->TM_cn[gr]];


        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
      }else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }


      th_p = enc->prd_threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
			pm = enc->pmodels[gr][cn];

      composite_MMF(enc,pm, x, y,u,freq_array,cumfreq_array,enc->limit_PRD);

      e = enc->encval[y][x];
      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和
      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;
		}//x fin
	}//y fin

ei[cn].cost = tot_cost;
ei[cn].cn = cn;
  }//cn

 for (i = 0; i <(  enc->num_pmodel  ) - 1; i++) {
        for (j = (  enc->num_pmodel  ) - 1; j > i; j--) {
            if (ei[j - 1].cost > ei[j].cost) {  // 前の要素の方が大きかったら 
                temp_ei = ei[j];        // 交換する 
                ei[j] = ei[j - 1];
                ei[j - 1]= temp_ei;
            }
        } 
    }
for(l = 0 ; l < enc->num_pmodel ;l++){
printf("|cost= %f,cn= %d|\n",ei[l].cost,ei[l].cn);
}


enc->PRD_cn = ei[0].cn;
printf("|PRD_cn = %d|\n",ei[0].cn);


fprintf(stderr,"...OK      \r");
printf("\n");



}


void set_cn(ENCODER *enc , int ***array){

	int x, y, e,i,l,j,k,c, base, bits,cn,gr, cumbase,flg;
  int multi;
  double cost,subcost,tot_cost,a;
  struct EI_Member ei[enc->num_pmodel ];
  struct EI_Member temp_ei;
  int *TM_cn;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

TM_cn = (int *)alloc_mem(sizeof(int) * enc->num_pmodel);

int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Set Cn...\n");
setbuf(stderr,NULL);

	bits = 0;

  bzero(&ei, sizeof(ei));
  bzero(&temp_ei, sizeof(temp_ei));

multi = enc->multi;

  k = 0;
 for (gr = 0; gr < enc->num_group; gr++) {
  l = 0;
      
    for(cn = 0 ; cn < enc->num_pmodel ;cn++){


  k = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			if (enc->TM_gr[y][x] == gr) {

 			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入
      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;
		
      }
		}//x fin
	}//y fin

ei[cn].cost = tot_cost;
ei[cn].cn = cn;
ei[cn].gr = gr;
printf("gr= %d cn= %d cost= %f\n",gr,cn,tot_cost);
  }//cn

 for (i = 0; i <(  enc->num_pmodel  ) - 1; i++) {
        for (j = (  enc->num_pmodel  ) - 1; j > i; j--) {
            if (ei[j - 1].cost > ei[j].cost) {  // 前の要素の方が大きかったら 
                temp_ei = ei[j];        // 交換する 
                ei[j] = ei[j - 1];
                ei[j - 1]= temp_ei;
            }
        } 
    }
for(l = 0 ; l < enc->num_pmodel ;l++){
printf("|cost= %f,cn= %d|\n",ei[l].cost,ei[l].cn);
}


TM_cn[gr] = ei[0].cn;
if((ei[0].cn == 0) && (ei[1].cn == 0)){TM_cn[gr] = ( 5 * enc->num_pmodel) / 8 - 1;}
printf("|TM_cn[%d]= %d|\n",gr,ei[0].cn);

}//gr

fprintf(stderr,"...OK      \r");
printf("\n");

enc->TM_cn = TM_cn;

for(gr = 0 ; gr < enc->num_group ;gr++){
  printf("|gr[%d]= %d|\n",gr,enc->TM_cn[gr]);
}


}


void calculate_cost_all(ENCODER *enc , int ***array, int flag){
////////基本使用しない。テスト用なので動くかも微妙
	int x, y, e,i,l,j,k,c, base, bits, gr,cn, cumbase,flg;
  int multi,th = 0;
  double cost,subcost,tot_cost,a;
  struct EI_Member ei[enc->num_pmodel * enc->num_group ];
  struct EI_Member temp_ei;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

int w_gr;
w_gr = enc->w_gr;
int w_func = 0;


printf("Start Calculate ALL...\n");
setbuf(stderr,NULL);

	bits = 0;

  bzero(&ei, sizeof(ei));
  bzero(&temp_ei, sizeof(temp_ei));

multi = enc->multi;
th = enc->tm_th;

  k = 0;
  l = 0;
      
  for(gr = 0 ; gr < enc->num_group ;gr++){
    for(cn = 0 ; cn < enc->num_pmodel ;cn++){


bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);
			pm = enc->pmodels[gr][cn];
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入
      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));


      tot_cost += cost + subcost;
		

		}//x fin
	}//y fin

ei[k].cost = tot_cost;
ei[k].gr = gr;
ei[k].cn = cn;
printf("cost[%d] = %f\n",k,tot_cost);
k++;


fprintf(stderr,"%3.1d / 256\r",l+1);
l++;

}//cn fin
}//gr fin
fprintf(stderr,"...OK      \r");
printf("\n");

 for (i = 0; i <(  enc->num_pmodel * enc->num_group  ) - 1; i++) {
        for (j = (  enc->num_pmodel * enc->num_group  ) - 1; j > i; j--) {
            if (ei[j - 1].cost > ei[j].cost) {  // 前の要素の方が大きかったら 
                temp_ei = ei[j];        // 交換する 
                ei[j] = ei[j - 1];
                ei[j - 1]= temp_ei;
            }
        } 
    }

printf("gr = %d  | cn = %d | multi = %d | th = %d\n",ei[0].gr,ei[0].cn,enc->multi,enc->tm_th);

enc->tm_cost = ei[0].cost;
gr = enc->gr = ei[0].gr;
cn = enc->cn = ei[0].cn;

return;
}


void calculate_cost(ENCODER *enc , int ***array, int flag){

	int x, y, e,l,k,c, base, bits, gr,cn, cumbase,flg;
  int multi,th = 0;
  double cost,subcost,tot_cost,a;
  int *th_p;
  int u;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Start Calculate Cost...\n");

	bits = 0;

multi = enc->multi;
cn = enc->cn;
th = enc->tm_th;

  k = 0;
  l = 0;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

      th_p = enc->threshold;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入

      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;

		}//x fin
	}//y fin

enc->tm_cost = tot_cost;
printf("%f bits/pel\n",tot_cost / (enc->width * enc->height));

return;
}

void calc_cn_select(ENCODER *enc , int ***array, int flag){

	int x, y, e,l,k,c, base, bits, gr,cn, cumbase,flg;
  int multi,th = 0;
  double cost,subcost,tot_cost,a;
  int *th_p;
  int u;

  double cn_cost[2];

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));
int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Select Cn...\n");

	bits = 0;

multi = enc->multi;
cn = enc->cn;
th = enc->tm_th;

  k = 0;
  l = 0;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);


for(l = 0;l < 2 ;l++){
tot_cost = 0;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      th_p = enc->threshold;
 
      if(y == 0 && x == 0){
      u = 0;
      }
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      if(l == 0){
      cn = enc->cn;//old
      }else if(l == 1){
      cn = enc->TM_cn[gr];//select
      }
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
      }else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

      e = enc->encval[y][x];//輝度値代入


      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;

		}//x fin
	}//y fin
cn_cost[l] = tot_cost;
printf("cn_cost= %f\n",tot_cost);
}
if(cn_cost[0] < cn_cost[1]){
  for(l = 0;l < enc->num_pmodel;l++){
    enc->TM_cn[l] = enc->cn;
  }
}
enc->tm_cost = tot_cost;

return;
}

void calc_hist_of_temp_and_prd(ENCODER *enc , int ***array){
//中嶋論文第5章ヒストグラムを求め近似曲線を算出する
	int x, y, e,l,k,c, base, bits, gr,cn, cumbase,flg;
  int multi,th = 0;
  double cost,subcost,tot_cost,a;
  int *th_p;
  int u;
  int temp;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

  int **pred_rate,**temp_rate;
pred_rate =  (int **)alloc_2d_array(enc->height ,enc->width , sizeof(int));
temp_rate =  (int **)alloc_2d_array(enc->height ,enc->width , sizeof(int));


int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Start Calculate Histgram...\n");

	bits = 0;

multi = enc->multi;
cn = enc->cn;
th = enc->tm_th;

  k = 0;
  l = 0;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

int num,prd;
int *roff_p;
int *org_p;
 
int *u_prd,*u_temp;
int *u_prd_res,*u_temp_res;
u_prd = (int *)alloc_mem( (MAX_UPARA + 2) * sizeof(int )); 
u_temp = (int *)alloc_mem( (MAX_UPARA + 2) * sizeof(int )); 

u_prd_res = (int *)alloc_mem( (MAX_UPARA + 2) * sizeof(int )); 
u_temp_res = (int *)alloc_mem( (MAX_UPARA + 2) * sizeof(int )); 



  for(k = 0;k < MAX_UPARA + 2 ; k++){
    u_temp[k] = u_prd[k] = 0;
  }

      th_p = enc->threshold;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      
      th_p = enc->threshold;
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入

      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;
      temp_rate[y][x] = (int)((cost + subcost) * COST_ACCURACY);

//////////
      th_p = enc->prd_threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
			pm = enc->pmodels[gr][cn];

    	roff_p = enc->roff[y][x];
      org_p = &enc->org[y][x];

      prd = 0;
      for (k = 0; k < PRED_AREA; k++) {
         prd += org_p[roff_p[k]] * enc->temp_predictor[y][x][k];
      }
        
	    base = enc->bconv[prd];
      num = enc->encval[y][x];

      pred_rate[y][x] = (int)((pm->cost[base + num] + pm->subcost[base]) * COST_ACCURACY);


/////////////////hist
      if(pred_rate[y][x] < temp_rate[y][x]){
        u_prd[u]++;
      }else{
        u_temp[u]++;
      }

		}//x fin
	}//y fin

  for(k = 0;k < MAX_UPARA + 2 ; k++){
    u_temp_res[k] = u_temp[k];
    u_prd_res[k] = u_prd[k];
  }
  for (k = 0; k < MAX_UPARA; k++) {
    for (l = MAX_UPARA; l > k; l--) {
      if (u_temp_res[l - 1] < u_temp_res[l]) {
        temp = u_temp_res[l];
        u_temp_res[l] = u_temp_res[l - 1];
        u_temp_res[l - 1] = temp;
      }
    } 
  }
  for (k = 0; k < MAX_UPARA; k++) {
    for (l = MAX_UPARA; l > k; l--) {
      if (u_prd_res[l - 1] < u_prd_res[l]) {
        temp = u_prd_res[l];
        u_prd_res[l] = u_prd_res[l - 1];
        u_prd_res[l - 1] = temp;
      }
    } 
  }
////各項目の端を決定

  for(l = 0; l < MAX_UPARA; l++){
    if(u_temp_res[0] == u_temp[l]){
      break;
    }
  }
int start_temp,end_temp,start_prd,end_prd;
  for(k = l;k > 0 ; k--){
    if(u_temp[k] == 0){
      break;
    }
  }
  start_temp = k + 1;
  for(k = l;k < MAX_UPARA ; k++){
    if(u_temp[k] == 0){
      break;
    }
  }
  end_temp = k - 1;
    //printf("start= %d,end= %d\n",start_temp,end_temp);

  for(l = 0; l < MAX_UPARA; l++){
    if(u_prd_res[0] == u_prd[l]){
      break;
    }
  }

  for(k = l;k > 0 ; k--){
    if(u_prd[k] == 0){
      break;
    }
  }
  start_prd = k + 1;
  for(k = l;k < MAX_UPARA ; k++){
    if(u_prd[k] == 0){
      break;
    }
  }
  end_prd = k - 1;

////最終的な端を決定
int start,end;
  if(start_temp < start_prd){
    start = start_prd;
  }else{
    start = start_temp;
  }  

  if(end_temp > end_prd){
    end = end_prd;
  }else{
    end = end_temp;
  }

////重みとなる係数を計算
int num_data = end - start;
double *fra;
int count;
fra = (double *)alloc_mem( (num_data + 1 ) * sizeof(double )); 
  count = 0;
  for(k = start; k <= end ;k++){
    fra[count] = (double)u_prd[k] / (double)u_temp[k];
    count++;
  } 
 
////近似曲線の算出
float tl_a = 0,tl_b = 0;
double sum_xy = 0, sum_x = 0, sum_y = 0, sum_x2 = 0;
  for(k = 0; k < num_data; k++){
    sum_xy += (double)k * fra[k];
    sum_x += (double)k;
    sum_y += fra[k];
    sum_x2 += (double)(k * k);
  }

  tl_a = (float)((num_data * sum_xy - sum_x * sum_y) / (num_data * sum_x2 - sum_x * sum_x));
  tl_b = (float)((sum_x2 * sum_y - sum_xy * sum_x) / (num_data * sum_x2 - sum_x * sum_x));

enc->trendline_a = tl_a;
enc->trendline_b = tl_b;
  printf("//////////////////\n");
  printf("a= %f\n",tl_a);
  printf("b= %f\n",tl_b);
  printf("//////////////////\n");

////近似曲線の結果保存
double *prd_coef;
prd_coef = (double *)alloc_mem( (num_data + 1 ) * sizeof(double )); 

  for(k = 0; k < num_data; k++){
    prd_coef[k] = tl_a * (double)k + tl_b;
  }
enc->prd_coef = prd_coef;
enc->coef_start = start;
enc->coef_end = end;
  printf("//////////////////\n");
  printf("start= %d\n",start);
  printf("end= %d\n",end);
  printf("//////////////////\n");


return;
}

void limit_PRD(ENCODER *enc , int ***array, int flag){
//論文未記入。Uの値0〜512からしきい値として一個数値を送り，しきい値以下なら第5章の予測の混合は行わない．
//これが無ければ符号化レートは悪くなる．実質近似曲線の範囲制限みたいなもん．
	int x, y, e,l,k,c, base, bits, gr,cn, cumbase,flg;
  int multi;
  double cost,subcost,tot_cost,a;
  int *th_p;
  int u;
  int g,h;
  double temp;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

  uint *freq_array_buf,*cumfreq_array_buf;
freq_array_buf =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));
cumfreq_array_buf =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));



int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Start Calculate Cost...\n");

	bits = 0;

multi = enc->multi;
cn = enc->cn;

  k = 0;
  l = 0;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

 
int i;
double *tcost,*tcost_buf;
tcost = (double *)alloc_mem( (MAX_UPARA + 1 ) * sizeof(double )); 
tcost_buf = (double *)alloc_mem( (MAX_UPARA + 1 ) * sizeof(double )); 

for(i = 0; i < MAX_UPARA; i++ ){
  tcost[i] = 0;
}

	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      
      th_p = enc->threshold;
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

      for(k = 0 ; k < pm->size + 1; k++){
        freq_array_buf[k] = freq_array[k];
        cumfreq_array_buf[k] = cumfreq_array[k];
      }
      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

			e = enc->encval[y][x];//輝度値代入
      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));
      tot_cost += cost + subcost;
 
      th_p = enc->prd_threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->cn;
			pm = enc->pmodels[gr][cn];

      for(i = 0 ; i < MAX_UPARA ; i++){
      composite_MMF(enc,pm, x, y,u,freq_array,cumfreq_array,i);

      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tcost[i] += cost + subcost;
      for(k = 0 ; k < 512 + 2; k++){
        freq_array[k] = freq_array_buf[k];
        cumfreq_array[k] = cumfreq_array_buf[k];
      }

      }
		}//x fin
	}//y fin
printf("%f bits/pel(original)\n",tot_cost / (enc->width * enc->height));
for(i = 0; i < MAX_UPARA; i++ ){
      tcost_buf[i] = tcost[i];
printf("tcost= %f\n",tcost[i]);
}

for (g = 0; g < MAX_UPARA - 1; g++) {
        for (h = MAX_UPARA - 1; h > g; h--) {
            if (tcost[h - 1] > tcost[h]) {  /* 前の要素の方が大きかったら */
                temp = tcost[h];        /* 交換する */
                tcost[h] = tcost[h - 1];
                tcost[h - 1] = temp;
            }
        } 
    }

for(i = 0; i < MAX_UPARA; i++ ){
  if(tcost[0] == tcost_buf[i]) break;
}

enc->limit_PRD = i;
printf("limit= %d\n",enc->limit_PRD);

return;
}




void calc_threshold_bypix(ENCODER *enc , int ***array, int flag){
//論文未記入．結果は研究経過報告で．Uの値0〜512の値をしきい値として一個伝送し，事例ベースと予測を画素毎に切り替える．
	int x, y, e,l,k,c, base, bits, gr,cn, cumbase,flg;
  int multi,th = 0;
  double cost,subcost,tot_cost,a;
  int *th_p;
  int u;
  int g,h,temp;
  int prd_buf,prd = 0;
  int prd_max = enc->maxval << 6;

	PMODEL *pm;
  uint **freq_array_save,*freq_array,*cumfreq_array;//*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

  int **pred_rate,**temp_rate;
  int *tot_rate,*buf_rate;
pred_rate =  (int **)alloc_2d_array(enc->height ,enc->width , sizeof(int));
temp_rate =  (int **)alloc_2d_array(enc->height ,enc->width , sizeof(int));

tot_rate =  (int *)alloc_mem((MAX_UPARA + 2) * sizeof(int));
buf_rate =  (int *)alloc_mem((MAX_UPARA + 2) * sizeof(int));

int w_gr;
w_gr = enc->w_gr;
int w_func = 1;

printf("Start Calculate Cost...\n");

	bits = 0;

multi = enc->multi;
cn = enc->cn;
th = enc->tm_th;

  k = 0;
  l = 0;

enc->rc->range = (range_t) -1;
enc->rc->code = 0;
bits = 0;
tot_cost = 0;

a = 1.0 / log(2.0);

int num;
int *roff_p;
int *org_p;
 


      th_p = enc->threshold;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      
      th_p = enc->threshold;
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      //cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);}
        else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

			e = enc->encval[y][x];//輝度値代入

      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;
      temp_rate[y][x] = (int)((cost + subcost) * COST_ACCURACY);

//////////
      th_p = enc->prd_threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->PRD_cn;
			pm = enc->pmodels[gr][cn];

    	roff_p = enc->roff[y][x];
      org_p = &enc->org[y][x];

      prd = 0;
      for (k = 0; k < PRED_AREA; k++) {
         prd += org_p[roff_p[k]] * enc->temp_predictor[y][x][k];
      }
        
	    base = enc->bconv[prd];
      num = enc->encval[y][x];

      pred_rate[y][x] = (int)((pm->cost[base + num] + pm->subcost[base]) * COST_ACCURACY);


		}//x fin
	}//y fin



for (th = 0; th < MAX_UPARA ; th++) {
  tot_rate[th] = 0;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
  
      u = enc->TM_U[y][x];      
      if(u < th){
      tot_rate[th] += temp_rate[y][x];
      }else{
      tot_rate[th] += pred_rate[y][x];
      }
        
    }
  }
}

for (th = 0; th < MAX_UPARA + 2; th++) {
  buf_rate[th] = tot_rate[th];
} 

for (g = 0; g < MAX_UPARA - 1; g++) {
        for (h = MAX_UPARA - 1; h > g; h--) {
            if (buf_rate[h - 1] > buf_rate[h]) {  /* 前の要素の方が大きかったら */
                temp = buf_rate[h];        /* 交換する */
                buf_rate[h] = buf_rate[h - 1];
                buf_rate[h - 1] = temp;
            }
        } 
    }

for(th = 0; th < MAX_UPARA + 2; th++ ){
  if(buf_rate[0] == tot_rate[th]) break;
}

enc->picsel_th = th;

tot_cost = 0;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      th_p = enc->threshold;
 
      if(y == 0 && x == 0){
      u = 0;
      }
      u = enc->TM_U[y][x];      
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->cn;
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
      }else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,0);
        }

      e = enc->encval[y][x];//輝度値代入

      if(u < enc->picsel_th){

      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      tot_cost += cost + subcost;

      }else{

      th_p = enc->prd_threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->PRD_cn;
			pm = enc->pmodels[gr][cn];

      prd_buf = prd;
      prd = 0;
      prd = enc->prd_b[y][x];
      if(prd < 0 || prd_max < prd ){prd = prd_buf;}
       
	    base = enc->bconv[prd];
      num = enc->encval[y][x];

      tot_cost += pm->cost[base + num] + pm->subcost[base];


      }
		}//x fin
	}//y fin

enc->tm_cost = tot_cost;

return;
}




struct TM_th
{
  double cost;
  int th;

};

struct th_and_multi
{
  double cost;
  int multi;
  int th;
  int number;
};


void calculate_encode(ENCODER *enc,int ***array )
{
  int rep;
  int gr;
  int *tm_count;
  uint **freq_array_save,*freq_array,*cumfreq_array;//,*sum_freq_array;
  
freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
tm_count =  (int *)alloc_mem(2 << TH_BITS * sizeof(int));//一次元配列

unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

setbuf(stderr,NULL);

  printf("Strating Calculate Cost\n");

double before_cost=0;
int count=0,i;
int threshold_save[enc->num_group];
int **TM_U_save;

TM_U_save = (int **)alloc_2d_array(enc->height ,enc->width ,  sizeof(int));


//最適化処理は基本的にここで行う．関数名にprdがあるものは画素適応予測のために使用している．
//予測をしたくなかったらprd関連をはずした後，encode_image内の書き換え(削除)を行えばいいはず．

for(gr = 0; gr < enc->num_group; gr++ ){
enc->PRD_cn = enc->cn = enc->TM_cn[gr] = ( 5 * enc->num_pmodel) / 8 - 1;//gauss_index(cn) = 9
}

// design_predictor_for_temp(enc, 1);
 enc->multi = MAX_MULTIMODAL;
 make_cost_u(enc , array);
 make_th(enc , array);
// make_prd_th(enc );
 set_cn_old(enc , array);
 w_gr_search(enc , array, 0);
 calculate_cost(enc , array, 0);
 before_cost = enc->tm_cost;
for(rep = 0 ; rep < REPETITION ; rep++){
 count++;

// design_predictor_for_temp(enc, 1);
 make_cost_u_2nd(enc , array);
 make_th(enc , array);
// make_prd_th(enc );
 set_cn_old(enc , array);
 w_gr_search(enc , array, 0);
 calculate_cost(enc , array, 0);
// calc_hist_of_temp_and_prd(enc , array);
 if((count != 1) && (enc->tm_cost >= before_cost)){

 for(i = 0;i < enc->num_group ;i++){
    threshold_save[i] = enc->threshold[i];
 }

 break;
 }
 before_cost = enc->tm_cost;
 }
 set_cn(enc , array);
// set_cn_prd(enc , array);

 calc_cn_select(enc , array, 0);

// limit_PRD(enc , array, 0);
// set_cn_prd_second(enc , array);
// make_prd_th_2nd(enc , array);

return;
}


int encode_image(FILE *fp, ENCODER *enc,int ***array )
{
	int x, y, e,k, base, bits, gr,cn, cumbase,flg;
  int multi;
  int *tm_count;
	PMODEL *pm ;
  uint **freq_array_save,*freq_array,*cumfreq_array;//,*sum_freq_array;
  
  int *th_p;
  int u;


freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
freq_array_save =  (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
//sum_freq_array =  (uint *)alloc_mem((512) * sizeof(uint));
tm_count =  (int *)alloc_mem(2 << TH_BITS * sizeof(int));//一次元配列

unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

setbuf(stderr,NULL);

////////////////////////////////////
/////////ファイルに書き込み/////////
////////////////////////////////////
int w_gr;
w_gr = enc->w_gr;
int w_func = 1;
cn = enc->cn;
multi = enc->multi = MAX_MULTIMODAL;

  printf("Strating Encode Image\n");

 int *roff_p;
  int *org_p;
  int c;
  double a,cost,subcost;
int area_cost_qt,count;
double area_cost;
int prd;
int **cost_save;
  cost_save = (int **)alloc_2d_array(  enc->height+1 , enc->width ,sizeof(int));
	cost_save[enc->height][0] = 0;

  double *wt_p;
  double a_cost;
a = 1.0 / log(2.0);
	wt_p = enc->ctx_weight_double;
bits = 0;
k = 0;
int prd_buf;
prd_buf = 0;
prd = 0;
double tot_cost;

	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      th_p = enc->threshold;

//特徴量Uの再算出 
      if(y == 0 && x == 0){
      u = 0;
      }else{
	
			roff_p = enc->roff[y][x];
      org_p = &cost_save[y][x];
      area_cost = 0;        
      count = 0;
      for(k=0;k < U_AREA; k++){
          a_cost = (double)org_p[roff_p[k]];
          a_cost = a_cost / COST_ACCURACY;
          if(a_cost == 0){
          a_cost = 0.1 * (wt_p[k]);
          }else{
          a_cost = a_cost * (wt_p[k]) ;
          }
 				  area_cost += a_cost; 
          if(a_cost != 0){
          count++;
          }
      }
        
      area_cost_qt = (int)(area_cost / count );
        if(area_cost_qt > MAX_UPARA){
          area_cost_qt = MAX_UPARA;
        }
       u = area_cost_qt;
      }

//確率モデリング
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,1);
      }else{
      flg = 1;
      Shift_freq(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg,1);
        }

//画素適応予測
/*
      th_p = enc->prd_threshold;
      for (gr = 0; gr < enc->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = enc->PRD_cn;
			pm = enc->pmodels[gr][cn];

      composite_MMF(enc,pm, x, y,u,freq_array,cumfreq_array,enc->limit_PRD);
*/
//符号化
      e = enc->encval[y][x];//輝度値代入



//printf("%d,%d,%d\n",x,y,e);








      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和

			rc_encode(fp, enc->rc,
				cumfreq_array[base + e] - cumbase,
				freq_array[base + e],
				cumfreq_array[base + enc->maxval + 1] - cumbase,
        1);//256個分に制限．rangecoder

//見積もり符号量の算出
			cumbase = cumfreq_array[base];//baseまでの高さの和
      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

      cost_save[y][x] = (int)((cost + subcost) * COST_ACCURACY);//Uのためにコストを保存

      tot_cost += cost + subcost;
		}//x fin
	}//y fin
	rc_finishenc(fp, enc->rc,1);
  
	bits += (int)enc->rc->code;







  free(freq_array_save);
  free(freq_array);
  free(cumfreq_array);
  free(sum_freq_array);

	return (bits);
}

#if AUTO_DEL_CL
void opcl_sub(ENCODER *enc, int k, int *blk, int tly, int tlx,
			  int blksize, int width, int level)
{
	int cl, min_cl, x, y, bry, brx;
	char **qtmap;
	cost_t *err_cost, cost, min_cost;

	brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
	bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;
	if (tly >= bry || tlx >= brx) return;
	if (level > 0) {
		qtmap = enc->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;
		if (qtmap[y][x] == 1) {
			blksize >>= 1;
			opcl_sub(enc, k, blk, tly, tlx, blksize, width, level - 1);
			opcl_sub(enc, k, blk, tly, tlx + blksize, blksize, width, level - 1);
			opcl_sub(enc, k, blk, tly + blksize, tlx, blksize, width, level - 1);
			opcl_sub(enc, k, blk, tly + blksize, tlx + blksize, blksize, brx, level - 1);
			return;
		}
	}
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx,
		blksize, width, enc->num_class);
	err_cost = enc->err_cost[*blk]; 
	min_cost = 1E5;
	min_cl = 0;
	for (cl = 0; cl < enc->num_class; cl++) {
		if (cl == k) continue;
		cost = err_cost[cl];
		cost += enc->class_cost[enc->mtfbuf[cl]];
		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}
	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			enc->class[y][x] = min_cl;
		}
	}
	(*blk)++;
}

cost_t opcl(ENCODER *enc, int k, int *blk, int restore)
{
	int x, y, i, j, *th_s, *prd_s, blksize, level;
	char *uq_s;
	cost_t cost;

	for (i = 0; i < enc->num_class; i++) enc->mtfbuf[i] = i;
	if( enc->optimize_loop > 1 && enc->quadtree_depth >= 0) {
		blksize = MAX_BSIZE;
	}else {
		blksize = BASE_BSIZE;
	}
	level =  enc->quadtree_depth;
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			opcl_sub(enc, k, blk, y, x, blksize, enc->width, level);
		}
	}
	th_s = (int *)alloc_mem(enc->num_group * sizeof(int));
	prd_s = (int *)alloc_mem(enc->max_prd_order * sizeof(int));
	uq_s = (char *)alloc_mem((MAX_UPARA + 1) * sizeof(char));
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			if (enc->class[y][x] > k) {
				enc->class[y][x]--;
			}
		}
	}
	for (i = 0; i < enc->num_group; i++) {
		th_s[i] = enc->th[k][i];
	}
	for (i = 0; i < enc->max_prd_order; i++) {
		prd_s[i] = enc->predictor[k][i];
	}
	for (i = 0; i < MAX_UPARA + 1; i++) {
		uq_s[i] = enc->uquant[k][i];
	}
	for (i = k; i < enc->num_class - 1; i++) {
		for (j = 0; j < enc->num_group; j++) {
			enc->th[i][j] = enc->th[i + 1][j];
		}
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->predictor[i][j] = enc->predictor[i + 1][j];
		}
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->nzconv[i][j] = enc->nzconv[i + 1][j];
		}
		enc->num_nzcoef[i] = enc->num_nzcoef[i + 1];
		for (j = 0; j < MAX_UPARA + 1; j++) {
			enc->uquant[i][j] = enc->uquant[i + 1][j];
		}
	}
#if AUTO_PRD_ORDER
	set_prd_pels(enc);
#endif
	enc->num_class--;
	predict_region(enc, 0, 0, enc->height, enc->width);
	cost = calc_cost(enc, 0, 0, enc->height, enc->width);
	cost += encode_class(NULL, enc, 0);
	cost += encode_predictor(NULL, enc, 0);
	cost += encode_threshold(NULL, enc, 0);
	if (restore == 1) {
		enc->num_class++;
		for (y = 0; y < enc->height; y ++) {
			for (x = 0; x < enc->width; x ++) {
				if (enc->class[y][x] >= k) {
					enc->class[y][x]++;
				}
			}
		}
		for (i = enc->num_class - 2; i >= k; i--) {
			for (j = 0; j < enc->num_group; j++) {
				enc->th[i + 1][j] = enc->th[i][j];
			}
			for (j = 0; j < enc->max_prd_order; j++) {
				enc->predictor[i + 1][j] = enc->predictor[i][j];
			}
			for (j = 0; j < MAX_UPARA + 1; j++) {
				enc->uquant[i + 1][j] = enc->uquant[i][j];
			}
		}
		for (i = 0; i < enc->num_group; i++) {
			enc->th[k][i] = th_s[i];
		}
		for (i = 0; i < enc->max_prd_order; i++) {
			enc->predictor[k][i] = prd_s[i];
		}
		for (i = 0; i < MAX_UPARA + 1; i++) {
			enc->uquant[k][i] = uq_s[i];
		}
#if AUTO_PRD_ORDER
		set_prd_pels(enc);
#endif
	}
	free(th_s);
	free(prd_s);
	free(uq_s);
	return(cost);
}

cost_t auto_del_class(ENCODER *enc, cost_t pre_cost)
{
	int x, y, i, j, k, del_cl, blk;
	cost_t cost, min_cost;
	char **class;

	class = (char **)alloc_2d_array(enc->height, enc->width,
		sizeof(char));
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			class[y][x] = enc->class[y][x];
		}
	}
	min_cost = 1E10;
	del_cl = 0;
	for (k = 0; k < enc->num_class; k++) {
		blk = 0;
		cost = opcl(enc, k, &blk, 1);
		if (cost < min_cost) {
			min_cost = cost;
			del_cl = k;
		}
	}
	if (pre_cost > min_cost) {
		blk = 0;
		cost = opcl(enc, del_cl, &blk, 0);
		for (i = del_cl; i < enc->num_class; i++) {
			for (j = 0; j < blk; j++) {
				enc->err_cost[j][i] = enc->err_cost[j][i + 1];
			}
		}
		printf("D");
		optimize_class(enc);
	} else {
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->class[y][x] = class[y][x];
			}
		}
	}
	predict_region(enc, 0, 0, enc->height, enc->width);
	cost = calc_cost(enc, 0, 0, enc->height, enc->width)
		+ encode_class(NULL, enc, 1) + encode_predictor(NULL, enc, 1) + encode_threshold(NULL, enc, 1);
	free(class);
	return(cost);
}
#endif


/* 比較関数 */


struct TM_Member   //データ
{
    int id;    //ID
    int by;   // 位置情報
    int bx;  // 位置情報
    int sum;    // TMの値
    int ave_o;

} ;


int comp( const void *a , const void *b ) {



const struct TM_Member* pmem1 = (struct TM_Member*)a;
const struct TM_Member* pmem2 = (struct TM_Member*)b;  
  
  int diff;

 diff = pmem1->sum - pmem2->sum;
        return diff;


  /* 引数はvoid*型と規定されているのでint型にcastする */
}
int*** TempleteM(ENCODER *enc,int ***array){
		int x, y, bx, by ,g,h,i, j, k,count,*area1 ,*area_o ,*tm_array , sum , m  , *roff_p, *org_p;
    
    clock_t start = 0,end = 0;
    start = clock();

    int x_size = X_SIZE;


//相関係数用
    int sum1,sum_o;
    double ave1,ave_o;
 
//nas用
double nas;

    

/////////////////////////
///////メモリ確保////////
/////////////////////////

	area1 = (int *)alloc_mem(AREA * sizeof(int  ));//一次元配列メモリ確保.初期化しようとしたらエラーが出た!
	area_o = (int *)alloc_mem(AREA * sizeof(int ));//一次元配列メモリ確保;
  tm_array = (int *)alloc_mem((Y_SIZE * X_SIZE * 2 + X_SIZE) * 4 * sizeof(int)) ;

  struct TM_Member tm[Y_SIZE * X_SIZE * 2 + X_SIZE ];

///////////////////////////
////////画像の走査/////////
///////////////////////////

  printf("     Start Calculating Templete Matching\n");
  for(y = 0 ; y < enc->height ; y++){
				for (x = 0; x < enc->width; x++){
  
        bzero(&tm, sizeof(tm));

				roff_p = enc->roff[y][x];//init_ref_offsetが入っている．予測器の範囲指定と番号付け
				org_p = &enc->org[y][x];

			  for(i=0;i < AREA; i++){//市街地距離AREA個分
          area1[i] = 0; 
				  area1[i] = org_p[roff_p[i]];
	      }
				
///////////////////////////
//テンプレートマッチング///
///////////////////////////

	j = 0;

if(y == 0 || y == 1 || y == 2){x_size = 50;
}else{
x_size = X_SIZE;
}

	for (by = y - Y_SIZE ; by <= y ; by++) {
		if((by < 0) || (by > enc->height))continue;
			for (bx = x - x_size ; bx <= x + x_size - 1; bx++) {
        
        if((by == y) && (bx == x) )break;
				if((bx < 0) || (bx > enc->width))continue;
				
      	sum = 0;
				roff_p = enc->roff[by][bx];
				org_p = &enc->org[by][bx];

				for(k=0;k < AREA; k++){//市街地距離AREA個分 
          area_o[k] = 0;
				  area_o[k] = org_p[roff_p[k]];
				}

        sum1 = 0;
        sum_o = 0;
 				for(m = 0; m < AREA ; m++){//平均の計算			
        sum1 += area1[m];
        sum_o += area_o[m];
      	}
        ave1 = (double)sum1 / AREA;
        ave_o = (double)sum_o / AREA;

nas = 0;  

				for(m = 0; m < AREA ; m++){//テンプレートマッチングの計算
				  nas += fabs( ((double)area1[m] - ave1) - ((double)area_o[m] - ave_o) );
				//	a = ( area1[m] - area_o[m] ) * ( area1[m] - area_o[m]);
				//	sum = sum + a;
				}
        tm[j].id = j;
        tm[j].by = by; 
        tm[j].bx = bx; 
        tm[j].ave_o = (int)ave_o;
        tm[j].sum = (int)(nas * NAS_ACCURACY);
       
		   			j++;	
			}//bx fin
    }//by fin
/////////////////////////
///////ソートの実行//////
/////////////////////////

enc->temp_num[y][x] = j;
 struct TM_Member temp;
for (g = 0; g < j - 1; g++) {
        for (h = j - 1; h > g; h--) {
            if (tm[h - 1].sum > tm[h].sum) {  /* 前の要素の方が大きかったら */
                temp = tm[h];        /* 交換する */
                tm[h] = tm[h - 1];
                tm[h - 1] = temp;
            }
        } 
    }

for(k = 0 ; k < Y_SIZE * X_SIZE * 2 + X_SIZE  ; k++){
  count = 0;
  tm_array[k * 4 + count] = 0;
  tm_array[k * 4 + count] = tm[k].id;
  count++;
  tm_array[k * 4 + count] = 0;
  tm_array[k * 4 + count] = tm[k].by; 
  count++;
  tm_array[k * 4 + count] = 0;
  tm_array[k * 4 + count] = tm[k].bx;
  count++;
  tm_array[k * 4 + count] = 0;
  tm_array[k * 4 + count] = tm[k].sum;
}


for(k = 0 ; k < MAX_DATA_SAVE ; k++){
array[y][x][k] = tm_array[k];

}

for(k = 0 ; k < MAX_DATA_SAVE_DOUBLE ; k++){
enc->array[y][x][k] = tm[k].ave_o;
}



}//x fin
}//y fin

end = clock();

printf("number of hours worked:%lf[s]\n",(float)(end - start)/CLOCKS_PER_SEC);


/////////////////////////////
////////メモリ解放///////////
////////////////////////////
/*
free(area1);
free(area_o);
free(tm_array);
*/
return(array);
}

/*
int encode_cn_prd(FILE *fp,ENCODER *enc){

int i,bits = 0;
  for(i = 0; i < enc->num_group ; i++){
    bits += putbits(fp,4,enc->PRD_cn[i]);
  }
return(bits);
}
*/
union floatUnion
{
  float f;
  char bytes[4];
};

int putfloat(FILE *fp, float x){
int i,bits = 0;
union floatUnion fu;

  fu.f = x;
  for(i = 0; i < 4; i++){
    putbits(fp,8,(uint)fu.bytes[i]);
    //putc(fu.bytes[i],fp);
    printf("%d\n",(uint)fu.bytes[i]);
    bits += 8;
  }
return(bits);
}

int encode_prd_factor(FILE *fp,ENCODER *enc){

int bits = 0;
    bits += putbits(fp,9,enc->coef_start);
    bits += putbits(fp,9,enc->coef_end);
    bits += putbits(fp,9,enc->limit_PRD);
    bits += putfloat(fp,enc->trendline_a);
    bits += putfloat(fp,enc->trendline_b);
printf("start %d\n",enc->coef_start);
printf("end %d\n",enc->coef_end);
printf("limit %d\n",enc->limit_PRD);
printf("trendline_a %f\n",enc->trendline_a);
printf("trendline_b %f\n",enc->trendline_b);
 
return(bits);
}


int encode_cn(FILE *fp,ENCODER *enc){

int i,bits = 0;
  for(i = 0; i < enc->num_group ; i++){
    bits += putbits(fp,4,enc->TM_cn[i]);
  }
return(bits);
}


void calc_cost_test(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
	cost_t cost;
	int x, y, u, cl, gr, prd, e, base, frac;
	int *upara_p, *prd_p, *encval_p;
	char *class_p, *group_p;
	PMODEL *pm;

	if (bry > enc->height) bry = enc->height;
	if (tlx < 0) tlx = 0;
	if (brx > enc->width) brx = enc->width;
	cost = 0;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		group_p = &enc->group[y][tlx];
		upara_p = &enc->upara[y][tlx];
		encval_p = &enc->encval[y][tlx];
		prd_p = &enc->prd[y][tlx];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			*upara_p++ = u = calc_uenc(enc, y, x);
			*group_p++ = gr = enc->uquant[cl][u];
			e = *encval_p++;
			prd = *prd_p++;
			prd = CLIP(0, enc->maxprd, prd);
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost += pm->cost[base + e] + pm->subcost[base];//コストの計算はここで行う.pm->costはfreq.pm->subcostはcumfreq
 //     printf("(%d,%d) u= %d cost= %f\n",x,y,u,pm->cost[base + e] + pm->subcost[base]);
		}
	}
}



int main(int argc, char **argv)
{
	cost_t cost, min_cost;

	int i,  k,/* x, y, xx, yy, cl,*/ bits, sw;
  int ***area_temp,***array;
//  double *bits_map;
	int header_info, class_info = 0, pred_info = 0, th_info = 0, err_info;
//	char **qtmap_save[QUADTREE_DEPTH];
//	int num_class_save;
//	char **class_save;
//	char **qtmap_save[QUADTREE_DEPTH];
	double rate;
	IMAGE *img;
	ENCODER *enc;
	PMODEL **pmlist_save;
	double elapse = 0.0;
	int f_mmse = 0;
	int f_optpred = 0;
	int quadtree_depth = QUADTREE_DEPTH;
	int num_class = NUM_CLASS;
	int num_group = NUM_GROUP;
	int prd_order = PRD_ORDER;
	int coef_precision = COEF_PRECISION;
	int num_pmodel = NUM_PMODEL;
	int pm_accuracy = PM_ACCURACY;
	int max_iteration = MAX_ITERATION;
	char *infile, *outfile;
	FILE *fp;
  
//自分用


	
  clock_t start = 0,end = 0;
  start = clock();

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
		case 'M':
			num_class = atoi(argv[++i]);
			if (num_class <= 0 || num_class > MAX_CLASS) {
				num_class = NUM_CLASS;
			}
			break;
		case 'K':
			prd_order = atoi(argv[++i]);
			if (prd_order <= 0 || prd_order > MAX_PRD_ORDER) {
				prd_order = PRD_ORDER;
			}
			break;
		case 'P':
			coef_precision = atoi(argv[++i]);
			if (coef_precision <= 0 || coef_precision > 16) {
				coef_precision = COEF_PRECISION;
			}
			break;
		case 'V':
			num_pmodel = atoi(argv[++i]);
			if (num_pmodel <= 0 || num_pmodel > 64) {
				num_pmodel = NUM_PMODEL;
			}
			break;
		case 'A':
			pm_accuracy = atoi(argv[++i]);
			if (pm_accuracy < 0 || pm_accuracy > 6) {
				pm_accuracy = PM_ACCURACY;
			}
			break;
		case 'I':
			max_iteration = atoi(argv[++i]);
			if (max_iteration <= 0) {
				max_iteration = MAX_ITERATION;
			}
			break;
		case 'm':
			f_mmse = 1;
			break;
		case 'o':
			f_optpred = 1;
			break;
		case 'f':
			quadtree_depth = -1;
			break;
		default:
			fprintf(stderr, "Unknown option: %s!\n", argv[i]);
			exit (1);
			}
		} else {
			if (infile == NULL) {
				infile = argv[i];
			} else {
				outfile = argv[i];
			}
		}
	}

	if (pm_accuracy > coef_precision) pm_accuracy = coef_precision;
	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", 0.01 * VERSION);
		printf("usage: encmrp [options] infile outfile\n");
		printf("options:\n");
		printf("    -M num  Number of predictors [%d]\n", num_class);
		printf("    -K num  Prediction order [%d]\n", prd_order);
		printf("    -P num  Precision of prediction coefficients (fractional bits) [%d]\n", coef_precision);
		printf("    -V num  Number of probability models [%d]\n", num_pmodel);
		printf("    -A num  Accuracy of probability models [%d]\n", pm_accuracy);
		printf("    -I num  Maximum number of iterations [%d]\n", max_iteration);
		printf("    -m      Use MMSE predictors\n");
		printf("    -f      Fixed block-size for adaptive prediction\n");
		printf("    -o      Further optimization of predictors (experimental)\n");
		printf("infile:     Input file (must be in a raw PGM format)\n");
		printf("outfile:    Output file\n");
		exit(0);
	}

	img = read_pgm(infile);

	fp = fileopen(outfile, "wb");
	

	k = img->width * img->height;
	if (num_class < 0) {
		num_class = (int)(10.4E-5 * k + 13.8);
		if (num_class > MAX_CLASS) num_class = MAX_CLASS;
	}
	if (prd_order < 0) {
		prd_order = (int)(12.0E-5 * k + 17.2);
		for (i = 1; i < 8; i++) {
			if (prd_order < (i + 1) * (i+1)) {
				prd_order = i * (i+1);
				break;
			}
		}
		if (i >= 8) prd_order = 72;
	}
#if AUTO_DEL_CL
	num_class = 63;
#endif
	printf("%s -> %s (%dx%d)\n", infile, outfile, img->width, img->height);
#if AUTO_PRD_ORDER
	printf("M = %d, K = %d, P = %d, V = %d, A = %d\n",
		num_class, MAX_PRD_ORDER, coef_precision, num_pmodel, pm_accuracy);
#else
	printf("M = %d, K = %d, P = %d, V = %d, A = %d\n",
		num_class, prd_order, coef_precision, num_pmodel, pm_accuracy);
#endif

	enc = init_encoder(img, num_class, num_group, prd_order, coef_precision,
		quadtree_depth, num_pmodel, pm_accuracy);
	enc->pmodels = init_pmodels(enc->num_group, enc->num_pmodel,
		enc->pm_accuracy, NULL, enc->sigma,
		enc->maxval + 1);//ガウス関数の初期化

#if defined(_WIN32) && (LOG_LIST_MODE || LOG_PUT_OUT_ENC)
	if (set_directory()) {
		exit(1);
	}
#endif

#if LOG_LIST_MODE
	init_log_sheet(enc, infile);
#endif

	set_cost_model(enc, f_mmse);


	
	init_class(enc);
	count_cl(enc);
/*
	num_class_save = num_class;
	prd_save = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order, sizeof(int));
	th_save = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
	class_save = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));
	if (enc->quadtree_depth > 0) {
		y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
		x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = enc->quadtree_depth - 1; i >= 0; i--) {
			qtmap_save[i] = (char **)alloc_2d_array(y, x, sizeof(char));
			y <<= 1;
			x <<= 1;
		}
	}
 */ 
	pmlist_save = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));

	/* 1st loop */
	enc->optimize_loop = 1;
	min_cost = INT_MAX;
/*
	for (i = j = 0; i < max_iteration; i++) {
		printf("[%2d] cost =", i);
		cost = design_predictor(enc, f_mmse);
		printf(" %d ->", (int)cost);
		cost = optimize_group(enc);
		printf(" %d ->", (int)cost);
		cost = optimize_class(enc);
		printf(" %d", (int)cost);
		if (cost < min_cost) {
			printf(" *\n");
			min_cost = cost;
			j = i;
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					class_save[y][x] = enc->class[y][x];
				}
			}
			for (cl = 0; cl < enc->num_class; cl++) {
				for (k= 0; k < enc->max_prd_order; k++) {
					prd_save[cl][k] = enc->predictor[cl][k];
				}
				for (k= 0; k < enc->num_group; k++) {
					th_save[cl][k] = enc->th[cl][k];
				}
			}
		} else {
			printf("\n");
		}
		if (i - j >= EXTRA_ITERATION) break;
		elapse += cpu_time();
	}
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->class[y][x] = class_save[y][x];
		}
	}
	for (cl = 0; cl < enc->num_class; cl++) {
		for (k= 0; k < enc->max_prd_order; k++) {
			enc->predictor[cl][k] = prd_save[cl][k];
		}
		i = 0;
		for (k= 0; k < enc->num_group; k++) {
			enc->th[cl][k] = th_save[cl][k];
			for (; i < enc->th[cl][k]; i++) {
				enc->uquant[cl][i] = k;
			}
		}
	}
*/
	set_cost_rate(enc);
#if AUTO_PRD_ORDER
//	set_prd_pels(enc);
#endif
	predict_region(enc, 0, 0, enc->height, enc->width);
	cost = calc_cost(enc, 0, 0, enc->height, enc->width);

  calc_cost_test(enc, 0, 0, enc->height,enc->width );
	enc->optimize_loop = 2;
	min_cost = INT_MAX;
	sw = 0;
/*
		printf("(%2d) cost =", i);
		if (f_optpred) {
//			cost = optimize_predictor(enc);
			printf(" %d", (int)cost);
		}
		side_cost = sc = encode_predictor(NULL, enc, 1);
		printf("[%d] ->", (int)sc);
//		cost = optimize_group(enc);
		side_cost += sc = encode_threshold(NULL, enc, 1);
		printf(" %d[%d] ->", (int)cost, (int)sc);
//		cost = optimize_class(enc);
		side_cost += sc = encode_class(NULL, enc, 1);
		printf(" %d[%d] (%d)", (int)cost, (int)sc, (int)side_cost);
		cost += side_cost;

	printf("cost = %d\n", (int)cost);
*/
	/* 2nd loop */
	enc->optimize_loop = 2;
	min_cost = INT_MAX;
	sw = 0;
/*	
	for (i = j = 0; i < max_iteration; i++) {
		printf("(%2d) cost =", i);
		if (f_optpred) {
			cost = optimize_predictor(enc);
			printf(" %d", (int)cost);
		}
		side_cost = sc = encode_predictor(NULL, enc, 1);
		printf("[%d] ->", (int)sc);
		cost = optimize_group(enc);
		side_cost += sc = encode_threshold(NULL, enc, 1);
		printf(" %d[%d] ->", (int)cost, (int)sc);
		cost = optimize_class(enc);
		side_cost += sc = encode_class(NULL, enc, 1);
		printf(" %d[%d] (%d)", (int)cost, (int)sc, (int)side_cost);
		cost += side_cost;
#if AUTO_DEL_CL
		if (sw != 0) {
			if( enc->num_class > 1) {
				sw = enc->num_class;
				cost = auto_del_class(enc, cost);
				while (sw != enc->num_class) {
					sw = enc->num_class;
					cost = auto_del_class(enc, cost);
				}
				printf(" -> %d[%d]", (int)cost, enc->num_class);
			}
		}
#endif
		if (cost < min_cost) {
			printf(" *\n");
			min_cost = cost;
			j = i;
			if (f_optpred) {
				num_class_save = enc->num_class;
				for (y = 0; y < enc->height; y++) {
					for (x = 0; x < enc->width; x++) {
						class_save[y][x] = enc->class[y][x];
					}
				}
				if (enc->quadtree_depth > 0) {
					y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
					x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
					for (k = enc->quadtree_depth - 1; k >= 0; k--) {
						for (yy = 0; yy < y; yy++) {
							for (xx = 0; xx < x; xx++) {
								qtmap_save[k][yy][xx] = enc->qtmap[k][yy][xx];
							}
						}
						y <<= 1;
						x <<= 1;
					}
				}
				for (gr = 0; gr < enc->num_group; gr++) {
					pmlist_save[gr] = enc->pmlist[gr];
				}
				for (cl = 0; cl < enc->num_class; cl++) {
					for (k= 0; k < enc->max_prd_order; k++) {
						prd_save[cl][k] = enc->predictor[cl][k];
					}
					for (k= 0; k < enc->num_group; k++) {
						th_save[cl][k] = enc->th[cl][k];
					}
				}
			}
		} else {
			sw = 1;
			printf("\n");
		}
		if (f_optpred) {
			if (i - j >= EXTRA_ITERATION) break;
		} else {
			if (i > j) break;
		}
		elapse += cpu_time();
	}
	if (f_optpred) {
		enc->num_class = num_class_save;
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->class[y][x] = class_save[y][x];
			}
		}
		if (enc->quadtree_depth > 0) {
			y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
			x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
			for (k = enc->quadtree_depth - 1; k >= 0; k--) {
				for (yy = 0; yy < y; yy++) {
					for (xx = 0; xx < x; xx++) {
						enc->qtmap[k][yy][xx] = qtmap_save[k][yy][xx];
					}
				}
				y <<= 1;
				x <<= 1;
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			enc->pmlist[gr] = pmlist_save[gr];
		}
		for (cl = 0; cl < enc->num_class; cl++) {
			for (k= 0; k < enc->max_prd_order; k++) {
				enc->predictor[cl][k] = prd_save[cl][k];
			}
			i = 0;
			for (k= 0; k < enc->num_group; k++) {
				enc->th[cl][k] = th_save[cl][k];
				for (; i < enc->th[cl][k]; i++) {
					enc->uquant[cl][i] = k;
				}
			}
		}
#if AUTO_PRD_ORDER
		set_prd_pels(enc);
#endif
		predict_region(enc, 0, 0, enc->height, enc->width);
		calc_cost(enc, 0, 0, enc->height, enc->width);
	}
	remove_emptyclass(enc);

#if AUTO_PRD_ORDER
	set_prd_pels(enc);	
#endif
*/

	bits = header_info = write_header(enc, fp);
	printf("header info.\t:%10d bits\n", header_info);

	enc->rc = rc_init();

//  encode_class(NULL,enc,0);
  area_temp = (int***)alloc_3d_array(enc->height , enc->width ,MAX_DATA_SAVE ,sizeof(int));
  array = TempleteM(enc ,area_temp);



  calculate_encode(enc,array );


//  bits += encode_prd_factor(fp,enc);
  bits += encode_cn(fp,enc);
//  putbits(fp,4,enc->PRD_cn);
  putbits(fp,4,enc->w_gr);
  putbits(fp, 7,0 );//rc_startdec用
  bits += encode_threshold_temp(fp, enc, 0);
//  bits += encode_prd_threshold(fp, enc, 0);
	bits += err_info = encode_image(fp, enc,array);//輝度値符号化に直す
	printf("TM_pred\t:%10d bits\n", err_info);
	printf("------------------------------\n");
	printf("total\t\t:%10d bits\n", bits);
	rate = (double)bits / (enc->height * enc->width);
	printf("coding rate\t:%10.5f b/p\n", rate);
	fclose(fp);
	elapse += cpu_time();
	printf("cpu time: %.2f sec.\n", elapse);


  end = clock();
  printf("Total worked times:%lf[s]\n",(float)(end - start)/CLOCKS_PER_SEC);

	enc->optimize_loop = -1;	// Coding Done
	calc_ratio_of_model_to_rate(enc);

#if LOG_LIST_MODE
	finish_log_sheet(enc, header_info, class_info, pred_info, th_info, err_info, bits, rate);
#endif

	print_rate_map(enc,array, outfile);
//	print_temp_map(enc,array, outfile);
//  print_multimodal(enc,array,outfile);
//    print_etc(enc,array, outfile);
#if LOG_PUT_OUT_ENC
//	print_predictor(enc->predictor, enc->max_prd_order, enc->num_class, enc->max_coef, outfile);
//	print_threshold(enc->th, enc->num_group, enc->num_class, enc->pmlist, NULL, outfile);
	//	print_class(enc->class, enc->num_class, enc->height, enc->width, outfile);
//	print_class_color(enc->class, enc->num_class, enc->height, enc->width, outfile);
//	print_block_size(enc->org, enc->qtmap, enc->quadtree_depth, enc->height, enc->width, outfile);
//	print_class_and_block(enc->class, enc->num_class, enc->qtmap, enc->quadtree_depth, enc->height, enc->width, outfile);
//	print_amp_chara(enc->predictor, enc->max_prd_order, enc->num_class, enc->height, enc->width, outfile);
//	calc_var_upara(enc, outfile);
#endif
	return (0);
}
