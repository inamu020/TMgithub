/***** Decoder *****/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mrp.h"

extern CPOINT dyx[];
extern double sigma_a[];
extern double qtree_prob[];
extern double zerocoef_prob[];

uint getbits(FILE *fp, int n)
{
	static int bitpos = 0;
	static uint bitbuf = 0;
	int x = 0;

	if (n <= 0) return (0);
	while (n > bitpos) {
		n -= bitpos;
		x = (x << bitpos) | bitbuf;
		bitbuf = getc(fp) & 0xff;
		bitpos = 8;
	}
	bitpos -= n;
	x = (x << n) | (bitbuf >> bitpos);
	bitbuf &= ((1 << bitpos) - 1);
	return (x);
}

int ***init_ref_offset(int height, int width, int prd_order)
///////////////////////////////
//roffのreturn値は座標ではありません．正確には何画素前かをマイナスの値で表します
//例えば当該の横なら-1，その上なら-512となります(例は512×512)
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








DECODER *init_decoder(FILE *fp)
{
	DECODER *dec;
	int i;

	dec = (DECODER *)alloc_mem(sizeof(DECODER));
	if (getbits(fp, 16) != MAGIC_NUMBER) {
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}
	dec->version = getbits(fp, 8);
	dec->width = getbits(fp, 16);
	dec->height = getbits(fp, 16);
	dec->maxval = getbits(fp, 16);
//	dec->num_class = getbits(fp, 6);
	dec->num_group = getbits(fp, 6);
	dec->max_prd_order = getbits(fp, 7);
	dec->num_pmodel = getbits(fp, 6) + 1;
	dec->coef_precision = getbits(fp, 4) + 1;
	dec->max_coef = (2 << dec->coef_precision);
	dec->pm_accuracy = getbits(fp, 3);
//	dec->quadtree_depth = (getbits(fp, 1))? QUADTREE_DEPTH : -1;
	dec->maxprd = dec->maxval << dec->coef_precision;

/*
printf("%d\n",MAGIC_NUMBER);
printf("%d\n",dec->version);
printf("%d\n",dec->width);
printf("%d\n",dec->height);
printf("%d\n",dec->maxval);
printf("%d\n",dec->num_group);
printf("%d\n",dec->num_pmodel);
printf("%d\n",dec->pm_accuracy);
*/

	dec->roff = init_ref_offset(dec->height, dec->width, dec->max_prd_order);

	dec->org = (int **)alloc_2d_array(dec->height+1, dec->width, sizeof(int));

	dec->org[dec->height][0] = (dec->maxval + 1) >> 1;
	dec->predictor = (int **)alloc_2d_array(dec->num_class, dec->max_prd_order,
		sizeof(int));
	dec->num_nzcoef = (int *)alloc_mem(dec->num_class * sizeof(int));
	dec->nzconv = (int **)alloc_2d_array(dec->num_class, dec->max_prd_order, sizeof(int));
	dec->th = (int **)alloc_2d_array(dec->num_class, dec->num_group - 1,
		sizeof(int));
	dec->err = (int **)alloc_2d_array(dec->height, dec->width, sizeof(int));
	dec->ctx_weight = init_ctx_weight();
	if (dec->quadtree_depth > 0) {
		int x, y, xx, yy;
		yy = (dec->height + MAX_BSIZE - 1) / MAX_BSIZE;
		xx = (dec->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = dec->quadtree_depth - 1; i >= 0; i--) {
			dec->qtmap[i] = (char **)alloc_2d_array(yy, xx, sizeof(char));
			for (y = 0; y < yy; y++) {
				for (x = 0; x < xx; x++) {
					dec->qtmap[i][y][x] = 0;
				}
			}
			yy <<= 1;
			xx <<= 1;
		}
	}
	dec->class = (char **)alloc_2d_array(dec->height, dec->width,
		sizeof(char));
	if (dec->num_pmodel > 1) {
		dec->pm_idx = (int *)alloc_mem(dec->num_group * sizeof(int));
	} else {
		dec->pm_idx = NULL;
	}
	dec->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	dec->spm.cumfreq = &(dec->spm.freq[MAX_SYMBOL]);

	dec->sigma = sigma_a;

	dec->mtfbuf = (int *)alloc_mem(dec->num_class * sizeof(int));
#if AUTO_PRD_ORDER
	dec->ord2mhd = (int *)alloc_mem(dec->max_prd_order * sizeof(int));
	for (i = 0; i < dec->max_prd_order; i++) {
		dec->ord2mhd[i] = (int)((sqrt(1 + 4 * i) - 1) / 2);
	}
	dec->prd_mhd = dec->ord2mhd[dec->max_prd_order - 1] + 1;
	dec->zero_fr = (int *)alloc_mem(NUM_ZMODEL * sizeof(int));
	for (i = 0; i < NUM_ZMODEL; i++) {
		dec->zero_fr[i] = (int)(zerocoef_prob[i] * (double)TOT_ZEROFR);
	}
#endif

  
dec->ctx_weight_double = init_ctx_weight_double();
dec->threshold = (int *)alloc_mem((dec->num_group) * sizeof(int));
dec->prd_threshold = (int *)alloc_mem((dec->num_group) * sizeof(int));
dec->TM_cn = (int *)alloc_mem((dec->num_group) * sizeof(int));
dec->temp_num = (int**)alloc_2d_array(  dec->height , dec->width ,sizeof(int));
dec->cost_save = (double**)alloc_2d_array(  dec->height , dec->width ,sizeof(double));	
dec->array = (int *)alloc_mem((MAX_DATA_SAVE_DOUBLE) * sizeof(int));
return (dec);
}
/*
#if AUTO_PRD_ORDER

void decode_predictor(FILE *fp, DECODER *dec)
{
	int k, cl, coef, sgn, d, zero_m, coef_m;
	PMODEL *pm;
	int b;

	pm = &dec->spm;
	pm->size = dec->max_coef + NUM_ZMODEL + 5;
	pm->cumfreq[dec->max_coef + 5] = 0;
	for(k = dec->max_coef + 5; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	b = dec->max_coef + 2;
	for (d = 0; d < dec->prd_mhd; d++) {
		zero_m = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + NUM_ZMODEL + 5)
			- (dec->max_coef + 5);
		coef_m = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + 13)
			- (dec->max_coef + 5);
		pm->cumfreq[b] = 0;
		pm->freq[b] = TOT_ZEROFR - dec->zero_fr[zero_m];
		pm->freq[b + 1] = dec->zero_fr[zero_m];
		pm->cumfreq[b + 1] = pm->freq[b];
		pm->cumfreq[b + 2] = TOT_ZEROFR;
		set_spmodel(pm, dec->max_coef + 1, coef_m);
		for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
			for (cl = 0; cl < dec->num_class; cl++) {
				coef = rc_decode(fp, dec->rc, pm, dec->max_coef + 2, dec->max_coef + 4)
					- (dec->max_coef + 2);
				if (coef == 1) {
					coef = rc_decode(fp, dec->rc, pm, 1, dec->max_coef + 1);
					sgn = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + 7)
						- (dec->max_coef + 5);
					if (sgn) {
						coef = -coef;
					}
					dec->predictor[cl][k] = coef;
				} else {
					dec->predictor[cl][k] = 0;
				}
			}
		}
	}
	for (cl = 0; cl < dec->num_class; cl++) {
		d = 0;
		for (k = 0; k < dec->max_prd_order; k++) {
			if (dec->predictor[cl][k] != 0) {
				dec->nzconv[cl][d++] = k;
			}
		}
		dec->num_nzcoef[cl] = d;
	}
	return;
}

#else

void decode_predictor(FILE *fp, DECODER *dec)
{
	int k, m, cl, coef, sgn, d;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = dec->max_coef + 18;
	pm->cumfreq[dec->max_coef + 2] = 0;
	for(k = dec->max_coef + 2; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	for (k = 0; k < dec->max_prd_order; k++) {
		m = rc_decode(fp, dec->rc, pm, dec->max_coef + 2, dec->max_coef + 18) - (dec->max_coef + 2);
		set_spmodel(pm, dec->max_coef + 1, m);
		for (cl = 0; cl < dec->num_class; cl++) {
			coef = rc_decode(fp, dec->rc, pm, 0, dec->max_coef + 1);
			if (coef > 0) {
				sgn = rc_decode(fp, dec->rc, pm, dec->max_coef+2, dec->max_coef+4)
					- (dec->max_coef + 2);
				if (sgn) {
					coef = -coef;
				}
			}
			dec->predictor[cl][k] = coef;
		}
	}
	for (cl = 0; cl < dec->num_class; cl++) {
		d = 0;
		for (k = 0; k < dec->max_prd_order; k++) {
			if (dec->predictor[cl][k] != 0) {
				dec->nzconv[cl][d++] = k;
			}
		}
		dec->num_nzcoef[cl] = d;
	}
	return;
}

#endif

void decode_threshold(FILE *fp, DECODER *dec)
{
	int cl, gr, m, k;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = 16;
	pm->cumfreq[0] = 0;
	for (k = 0; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	m = rc_decode(fp, dec->rc, pm, 0, pm->size);
	set_spmodel(pm, MAX_UPARA + 2, m);
	for (cl = 0; cl < dec->num_class; cl++) {
		k = 0;
		for (gr = 1; gr < dec->num_group; gr++) {
			if (k <= MAX_UPARA) {
				k += rc_decode(fp, dec->rc, pm, 0, pm->size - k);
			}
			dec->th[cl][gr - 1] = k;
		}
	}

	if (dec->num_pmodel > 1) {
		pm->size = dec->num_pmodel;
		pm->freq[0] = 0;
		for (k = 0; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}
		for (gr = 0; gr < dec->num_group; gr++) {
			dec->pm_idx[gr] = rc_decode(fp, dec->rc, pm, 0, pm->size);
		}
	}
	return;
}

void decode_qtindex(FILE *fp, DECODER *dec, PMODEL *cpm,
					int tly, int tlx, int blksize, int width, int level)
{
	int i, cl, y, x, bry, brx, ctx;
	char **qtmap;
	PMODEL *pm;

	brx = (tlx + blksize < dec->width) ? (tlx + blksize) : dec->width;
	bry = (tly + blksize < dec->height) ? (tly + blksize) : dec->height;
	if (tlx >= brx || tly >= bry) return;
	if (level > 0) {
		ctx = 0;
		qtmap = dec->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;
		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}
		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;
		ctx = ((level - 1) * 4 + ctx) << 1;

		pm = &dec->spm;
		i = rc_decode(fp, dec->rc, pm, ctx, ctx + 2) - ctx;

		if (i == 1) {
			qtmap[y][x] = 1;
			blksize >>= 1;
			decode_qtindex(fp, dec, cpm, tly, tlx,
				blksize, width, level - 1);
			decode_qtindex(fp, dec, cpm, tly, tlx + blksize,
				blksize, width, level - 1);
			decode_qtindex(fp, dec, cpm, tly + blksize, tlx,
				blksize, width, level - 1);
			decode_qtindex(fp, dec, cpm, tly + blksize, tlx + blksize,
				blksize, brx, level - 1);
			return;
		}
	}
	i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);

	mtf_classlabel(dec->class, dec->mtfbuf, tly, tlx,
		blksize, width, dec->num_class);
	for (cl = 0; cl < dec->num_class; cl++) {
		if (dec->mtfbuf[cl] == i) break;
	}
	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			dec->class[y][x] = cl;
		}
	}
	return;
}

void decode_class(FILE *fp, DECODER *dec)
{
	int i, j, x, y, blksize, level;
	PMODEL *pm, cpm[1];
	double p;
	int ctx;
	int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];

	if (dec->quadtree_depth >= 0) {
		level = dec->quadtree_depth;
		blksize = MAX_BSIZE;
	} else {
		level = 0;
		blksize = BASE_BSIZE;
	}

	pm = &dec->spm;
	if (dec->quadtree_depth > 0) {
		set_spmodel(pm, 7, -1);
		for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
			qtree_code[ctx] = rc_decode(fp, dec->rc, pm, 0, pm->size);
		}
	}
	set_spmodel(pm, PMCLASS_LEVEL, -1);
	for (i = 0; i < dec->num_class; i++) {
		mtf_code[i] = rc_decode(fp, dec->rc, pm, 0, pm->size);
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
	// set prob. models 
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
	}
	cpm->size = dec->num_class;
	cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
	cpm->cumfreq = &cpm->freq[cpm->size];
	cpm->cumfreq[0] = 0;
	for (i = 0; i < dec->num_class; i++) {
		p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
			* PMCLASS_MAX/PMCLASS_LEVEL);
		cpm->freq[i] = (uint)(p * (1 << 10));
		if (cpm->freq[i] <= 0) cpm->freq[i] = 1;
		cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
	}

	for (i = 0; i < dec->num_class; i++) {
		dec->mtfbuf[i] = i;
	}
	for (y = 0; y < dec->height; y += blksize) {
		for (x = 0; x < dec->width; x += blksize) {
			decode_qtindex(fp, dec, cpm, y, x,
				blksize, dec->width, level);
		}
	}
	return;
}


int calc_udec(DECODER *dec, int y, int x)
{
	int rx, ry, k, u;
	int **err, *wt_p;

	u = 0;
	err = dec->err;
	wt_p = dec->ctx_weight;
	if (y > UPEL_DIST && x > UPEL_DIST && x <= dec->width - UPEL_DIST) {
		for (k = 0; k < NUM_UPELS; k++) {
			ry = y + dyx[k].y;
			rx = x + dyx[k].x;
			u += err[ry][rx] * (*wt_p++);
		}
	} else if (y == 0) {
		if (x == 0) {
			for (k = 0; k < NUM_UPELS; k++) {
				u += ((dec->maxval + 1) >> 2) * (*wt_p++);
			}
		} else {
			ry = 0;
			for (k =0; k < NUM_UPELS; k++) {
				rx = x + dyx[k].x;
				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;
				u += err[ry][rx] * (*wt_p++);
			}
		}
	} else {
		if (x == 0) {
			for (k = 0; k < NUM_UPELS; k++) {
				ry = y + dyx[k].y;
				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;
				rx = x + dyx[k].x;
				if (rx < 0) rx = 0;
				u += err[ry][rx] * (*wt_p++);
			}
		} else {
			for (k = 0; k < NUM_UPELS; k++) {
				ry = y + dyx[k].y;
				if (ry < 0) ry = 0;
				rx = x + dyx[k].x;
				if (rx < 0) rx = 0;
				else if (rx >= dec->width) rx = dec->width - 1;
				u += err[ry][rx] * (*wt_p++);
			}
		}
	}
	u >>= 6;
	if (u > MAX_UPARA) u = MAX_UPARA;
	return (u);
}

int calc_prd(IMAGE *img, DECODER *dec, int cl, int y, int x)
{
	int k, prd, prd_order, rx, ry, *coef_p, *nzc_p, i;

	prd_order = dec->num_nzcoef[cl];
	prd = 0;
	coef_p = dec->predictor[cl];
	nzc_p = dec->nzconv[cl];
	if (y == 0) {
		if (x == 0) {
			for (k = 0; k < prd_order; k++) {
				prd += coef_p[nzc_p[k]];
			}
			prd *= ((img->maxval + 1) >> 1);
		} else {
			ry = 0;
			for (k = 0; k < prd_order; k++) {
				i = nzc_p[k];
				rx = x + dyx[i].x;
				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;
				prd += coef_p[i] * img->val[ry][rx];
			}
		}
	} else {
		if (x == 0) {
			for (k = 0; k < prd_order; k++) {
				i = nzc_p[k];
				ry = y + dyx[i].y;
				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;
				rx = x + dyx[i].x;
				if (rx < 0) rx = 0;
				prd += coef_p[i] * img->val[ry][rx];
			}
		} else {
			for (k = 0; k < prd_order; k++) {
				i = nzc_p[k];
				ry = y + dyx[i].y;
				if (ry < 0) ry = 0;
				rx = x + dyx[i].x;
				if (rx < 0) rx = 0;
				else if (rx >= img->width) rx = img->width - 1;
				prd += coef_p[i] * img->val[ry][rx];
			}
		}
	}
        prd = CLIP(0, dec->maxprd, prd);
	return (prd);
}

*/

struct TM_Member   //データ
{
    int id;    //ID
    int by;   // 位置情報
    int bx;  // 位置情報
    int sum;    // TMの値
    int ave_o;
} ;



void TempleteM(IMAGE *img ,DECODER *dec,int dec_x ,int dec_y , int *array, int *area1 ,int *area_o ){


int a = 0,x = 0, y = 0, bx, by ,g,h, i, j, k,count/*,*area1 ,*area_o, *tm_array*/ ,sum, m  , *roff_p, *org_p;

    int x_size = X_SIZE;
  struct TM_Member tm[Y_SIZE * X_SIZE * 2 + X_SIZE ];
  int tm_array[(Y_SIZE * X_SIZE * 2 + X_SIZE)*4] = {0};

//相関係数用
    int sum1,sum_o;
    double ave1,ave_o;

//nas用
double nas;



///////////////////////////
////////画像の走査/////////
///////////////////////////
if(dec_y == 0 && dec_x == 0){

//return;

}else{
  for (y = 0; y <= dec_y ; y++) {
    if(dec_y == y){
		  for (x = 0; x < dec_x ; x++) {
			  dec->org[y][x] = img->val[y][x];
      }
    }else{
  	  for (x = 0; x < dec->width ; x++) {
		  	dec->org[y][x] = img->val[y][x];
      }
	  }	
	}//fin dec_y
}

bzero(&tm, sizeof(tm));


roff_p = dec->roff[dec_y][dec_x];
org_p = &dec->org[dec_y][dec_x];

for(i=0;i < AREA; i++){//市街地距離AREA個分
  area1[i] = 0; 
	area1[i] = org_p[roff_p[i]];
}
				
///////////////////////////
//テンプレートマッチング///
///////////////////////////

j = 0;
if(dec_y == 0 || dec_y == 1 || dec_y == 2){
  x_size = 50;
}else{
  x_size = X_SIZE;
}
			
	for (by = dec_y - Y_SIZE ; by <= dec_y ; by++) {//参照領域　参考テンプレ．縦に動いた分．走査→○,domain
		if((by < 0) || (by > dec->height))continue;//画像の範囲外の場合の設定
			for (bx = dec_x - x_size ; bx <= dec_x + x_size - 1; bx++) {//横に動いた分
        
        if((by == dec_y) && (bx == dec_x) )break;
				if((bx < 0) || (bx > dec->width))continue;
				
				roff_p = dec->roff[by][bx];//init_ref_offsetが入っている．予測器の範囲指定と番号付け
				org_p = &dec->org[by][bx];

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




      	sum = 0;
        nas = 0;  
				for(m = 0; m < AREA ; m++){//テンプレートマッチングの計算
				  nas += fabs( ((double)area1[m] - ave1) - ((double)area_o[m] - ave_o) );
//					a = ( area1[m] - area_o[m] ) * ( area1[m] - area_o[m]);
//					sum = sum + a;
				}

        if(sum == 0){sum = 1;}
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

dec->temp_num[dec_y][dec_x] = j;

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
  tm_array[k * 4 + count] = tm[k].id;
  count++;
  tm_array[k * 4 + count] = tm[k].by; 
  count++;
  tm_array[k * 4 + count] = tm[k].bx;
  count++;
  tm_array[k * 4 + count] = tm[k].sum;
}

for(k = 0 ; k < MAX_DATA_SAVE ; k++){
array[k] = tm_array[k];//用意した一次元配列に入れよう

}

for(k = 0 ; k < MAX_DATA_SAVE_DOUBLE ; k++){
dec->array[k] = tm[k].ave_o;
}


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



 void design_predictor_for_temp(DECODER *dec, int x,int y,int u,int f_mmse)
{
  double **mat, *weight, w, e, d, pivot;
  int i, j, k, gr, pivpos, *index, *roff_p, *org_p;
  int *predictor;
  int *th_p;
  struct point pos;
  int num,bx,by;

  predictor = (int*)alloc_mem(sizeof(int) * PRED_AREA);

  mat = (double **)alloc_2d_array(PRED_AREA, PRED_AREA + 1, sizeof(double));
  index = (int *)alloc_mem(sizeof(int) * PRED_AREA);
  weight = (double *)alloc_mem(sizeof(double) * dec->num_group);

  for (gr = 0; gr < dec->num_group; gr++) {//gr は分散の番号16種類
    if (f_mmse) {//f_mmseは基本0
      weight[gr] = 1.0;
    } else {
      weight[gr] = 1.0 / (dec->sigma[gr] * dec->sigma[gr]);
    }
  }
    for (i = 0; i < PRED_AREA; i++) {
      for (j = 0; j <= PRED_AREA; j++) {
        mat[i][j] = 0.0;//初期化
      }
    }

        if(f_mmse){
          gr = 5;
        }else{

          th_p = dec->threshold;
          for (gr = 0; gr < dec->num_group; gr++) {
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
          if((by < 0) || (by > dec->height))continue;
          if((bx < 0) || (bx > dec->width))continue;

          roff_p = dec->roff[by][bx];
          org_p = &dec->org[by][bx];
          for (i = 0; i < PRED_AREA; i++) {
            w = weight[gr] * org_p[roff_p[i]];
              for (j = i; j < PRED_AREA; j++) {
                mat[i][j] += w * org_p[roff_p[j]];
              }
            mat[i][PRED_AREA] += w * org_p[0];//当該の答え
          }
        }

        for (i = 0; i < PRED_AREA; i++) {
          index[i] = i;
          for (j = 0; j < i; j++) {
            mat[i][j] = mat[j][i];//対称行列にする
          }
        }
        for (i = 0; i < PRED_AREA; i++) {
          pivpos = i;
          pivot = fabs(mat[index[i]][i]);//浮動小数点の絶対値化
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
        w = (1 << dec->coef_precision);
        e = 0.0;
        for (i = 0; i < PRED_AREA; i++) {
          if (fabs(mat[index[i]][i]) > 1E-10) {
            d = mat[index[i]][PRED_AREA] * w;
          } else {
            d = 0.0;
          }
          k = (int)d;
          if (k > d) k--;
          if (k < -dec->max_coef) {
            d = k = -dec->max_coef;
          } else if (k > dec->max_coef) {
            d = k = dec->max_coef;
          }
          predictor[i] = k;
          d -= k;
          e += d;
          mat[index[i]][PRED_AREA] = d;
        }

  dec->temp_predictor = predictor;


  free(weight);
  free(index);
  free(mat);

}





double continuous_GGF(DECODER *dec, double e){

int lngamma(double);
int cn,num_pmodel = dec->num_pmodel;
double sigma,delta_c,shape,eta,p;
double accuracy = 1 / (double)NAS_ACCURACY;
  cn = WEIGHT_CN; 
  sigma = dec->sigma[dec->w_gr];

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

void composite_MMF(DECODER *dec, PMODEL *pm,int x,int y,int u,
 uint *freq_shift_array,uint *cumfreq_shift_array,int limit){

int i,prd;
int start = dec->coef_start,end = dec->coef_end;
double coef = 0;
uint *freq_array,*array_original;

  if(u > limit){

  if( (u < start) || (end < u) ){
    if( u < start ){return;;
    }else if( end < u ){coef = dec->prd_coef[end - start - 1];
    }else{printf("error - accidental U\n");}
  }else{
  coef = dec->prd_coef[u - start];
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

  prd = dec->prd >> 6;
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



int Shift_freq(IMAGE *img ,PMODEL *pm,  DECODER *dec,int x,int y,
int *a,uint **freq_shift_array_save,uint *freq_shift_array,uint *cumfreq_shift_array, unsigned long long int *sum_freq_shift_array,
int multimodal ,int flg){




  int i,j,by[multimodal],bx[multimodal],e[multimodal];
  int flag2 = 0;
  double er_sum = 0.0,er_coef = 0.0;
  int w[multimodal];
  int w_cn; 

  int flg_min = 0;
  int mc[multimodal];
  double mc_w[multimodal];
//  int multi = multimodal;


int multi = multimodal - 1;//稲村






//T.Mの値を参照する

//fundamental用
int *org_p,*roff_p;
int area1_sum;
double area1_ave;

int w_gr;
int mc_int;
double mc_double;

if(flg == 0){
for(i = 0 ; i < multi ; i++){//いくつの山にするか
    by[i] = a[i*4+1];
    bx[i] = a[i*4+2];  

    e[i] = 0;

 for(j = 0;j < pm->size ; j++){
    freq_shift_array[j+e[i]]  = pm->freq[j];
  }
}

}


else{

if( (x < 5) && (y == 0)){//最初の5画素については山の数は画素分のみ

multi = x;
flag2 = 1;
flg_min = 1;

}
if((flg_min != 1) && ( dec->temp_num[y][x] < multi) ){//それ以降はTMの値の個数に依存

multi = dec->temp_num[y][x];
flag2 = 0;

}

for(j = 0;j < pm->size  ; j++){//初期化
    for(i = 0;i < multi ; i++){
      freq_shift_array_save[i][j] = 1;
    }
    freq_shift_array[j] = 1;
sum_freq_shift_array[j] = 0;
}

	  roff_p = dec->roff[y][x];
		org_p = &dec->org[y][x];
    area1_sum = 0;
    area1_ave = 0;
		for(i=0;i < AREA; i++){//市街地距離AREA個分
			area1_sum += org_p[roff_p[i]];
	  }
    area1_ave = area1_sum / AREA;


for(i = 0 ; i < multi ; i++){//いくつの山にするか
    by[i] = a[i*4+1];
    bx[i] = a[i*4+2];  


    e[i] = (int)(img->val[by[i]][bx[i]] - dec->array[i] + area1_ave);
    if(e[i] < 0 || e[i] > dec->maxval){e[i] = area1_ave;}
    mc[i] = a[i*4+3];

 for(j = 0;j < pm->size ; j++){
    freq_shift_array_save[i][j+e[i]]  = pm->freq[j];
    }
  }//multi fin






//puts("t1");




if (y == 0 || x == 0){


//puts("t2");





  by[multi] = a[multi*4+1];
    bx[multi] = a[multi*4+2];  


    e[multi] = (int)(img->val[by[multi]][bx[multi]] - dec->array[multi] + area1_ave);
    if(e[multi] < 0 || e[multi] > dec->maxval){e[multi] = area1_ave;}
    mc[multi] = a[multi*4+3];

 for(j = 0;j < pm->size ; j++){
    freq_shift_array_save[multi][j+e[multi]]  = pm->freq[j];
    }


//printf("em%d\n",e[multi]);

multi++;

}

else{


//puts("t3");





    by[multi] = a[multi*4+1];
    bx[multi] = a[multi*4+2];  


    e[multi] = (int)((dec->org[y - 1][x] + dec->org[y][x - 1] +  dec->org[y - 1][x - 1]) / 3);//稲村変更
    if(e[multi] < 0 || e[multi] > dec->maxval){e[multi] = area1_ave;}
    mc[multi] = a[multi*4+3];//稲村変更



//printf("%d,%d,%lf\n",img->val[by[multi]][bx[multi]],dec->array[multi],area1_ave);



 for(j = 0;j < pm->size ; j++){
    freq_shift_array_save[multi][j+e[multi]]  = pm->freq[j];
    }






multi++;   //稲村追加

}












///////////////////////////
//重み付け，正規化処理/////
///////////////////////////

if(flag2 == 1){

for(i = 0 ; i < multi ; i++){
  for(j = 0;j < pm->size ; j++){
    sum_freq_shift_array[j] += freq_shift_array_save[i][j];
  }
}

for(j = 0;j < pm->size ; j++){
  sum_freq_shift_array[j] = sum_freq_shift_array[j] / multi;
}

}else{

#ifdef VARIABLE_WEIGHT

//////////////////////////
////ガウス分布を求める////
//////////////////////////



w_gr = dec->w_gr;
w_cn = WEIGHT_CN; //とりあえずmrp.hで固定しておく

for(i = 0 ; i < multi ; i++){

mc_double = 0;
if(mc[i] != 0){
mc_double = (double)mc[i] / NAS_ACCURACY; 
}
mc_int = (int)mc_double; 
mc_w[i] = continuous_GGF(dec, mc_double);

er_sum = er_sum + mc_w[i];
}
if(er_sum == 0){er_sum = 1;}


er_coef = (100) / er_sum;//er_coefが0になるのを未然に防ぐ
for(i = 0 ; i < multi ; i++){
  w[i] = (int)(mc_w[i] * er_coef);
}





#endif


#ifdef RECIPROCAL_WEIGHT


for(i = 0 ; i < multi ; i++){//いくつの山にするか
    mc_w[i] = 1 / (double)mc[i] / NAS_ACCURACY;
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

#endif

for(i = 0 ; i < multi ; i++){
  for(j = 0;j < pm->size ; j++){
     freq_shift_array_save[i][j] = freq_shift_array_save[i][j] * w[i];//全てに重みをかけて保存
  }
  
  for(j = 0;j < pm->size ; j++){
    sum_freq_shift_array[j] += freq_shift_array_save[i][j];//重みをつけたものをsumに入れていく
  }
}//multi fin

for(j = 0;j < pm->size ; j++){
  sum_freq_shift_array[j] = sum_freq_shift_array[j] / 100;
  if(sum_freq_shift_array[j] == 0){
    sum_freq_shift_array[j] = 1;
  }
}




}//else fin
 
  for(j = 0;j < pm->size ; j++){
     pm->freq[j] = sum_freq_shift_array[j] ; //ずらした結果を代入
    
  }

 
}//else(flg == 0) fin


pm->cumfreq[0] = 0;
	for (i = 0; i < pm->size; i++) {
		pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
	}

return(0);
}







void convert_roff_value(DECODER *dec,int r ,int *roff_x,int *roff_y){
//roffの値から実際の座標を割り出します．
int roff;
int quotient,scs,k;
int width = dec->width;

roff = r;
quotient = roff / width;
scs = roff % width;
k = scs + width;
if((width / 2 < k) && (width >= k)){k = k - width;}
*roff_y = quotient;
*roff_x = k;

return;
}




void decode_threshold(FILE *fp, DECODER *dec)
{
	int gr, m, k,l;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = 16;
	pm->cumfreq[0] = 0;
	for (k = 0; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	m = rc_decode(fp, dec->rc, pm, 0, pm->size);
	set_spmodel(pm, MAX_UPARA + 2, 0);
	k = 0;
  dec->threshold[0] = 0;
		for (gr = 1; gr < m  ; gr++) {
			if (k <= MAX_UPARA) {
				k = rc_decode(fp, dec->rc, pm, 0, pm->size);
			}
			dec->threshold[gr] = k;
	}
  for(l = m ;l < dec->num_group ;l++ ){
    dec->threshold[l] = MAX_UPARA + 1;
  }
   for(gr = 0 ;gr < dec->num_group ;gr++ ){
    printf("threshold[%d] = %d\n",gr,dec->threshold[gr] );
  }
 

	return;
}

void decode_prd_threshold(FILE *fp, DECODER *dec)
{
	int gr, m, k,l;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = 16;
	pm->cumfreq[0] = 0;
	for (k = 0; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	m = rc_decode(fp, dec->rc, pm, 0, pm->size);
	set_spmodel(pm, MAX_UPARA + 2, 0);
	k = 0;
  dec->prd_threshold[0] = 0;
		for (gr = 1; gr < m  ; gr++) {
			if (k <= MAX_UPARA) {
				k = rc_decode(fp, dec->rc, pm, 0, pm->size);
			}
			dec->prd_threshold[gr] = k;
	}
  for(l = m ;l < dec->num_group ;l++ ){
    dec->prd_threshold[l] = MAX_UPARA + 1;
  }
   for(gr = 0 ;gr < dec->num_group ;gr++ ){
    printf("threshold[%d] = %d\n",gr,dec->prd_threshold[gr] );
  }
 

	return;
}



void set_up_prd(DECODER *dec)
{
  int i,k;
  int shift;

	dec->bconv = (img_t *)alloc_mem((dec->maxprd + 1) * sizeof(img_t));
	shift = dec->coef_precision - dec->pm_accuracy;
	for (k = 0; k <= dec->maxprd; k++) {
		i = (dec->maxprd - k + (1 << shift) / 2) >> shift;
		dec->bconv[k] = (i >> dec->pm_accuracy);
	}

}


IMAGE *decode_image(FILE *fp, DECODER *dec)
{
  int x, y, gr , cn, prd = 0, u = 0, j, p, k, base,flg,multi;
  int *th_p, *array, *area1, *area_o;
  uint *freq_array, *cumfreq_array , **freq_array_save;
  IMAGE *img;
  PMODEL *pm,*pm_save,*pm_s;

  int cumbase,c;
  double a,cost,subcost;
  unsigned long long int *sum_freq_array;


cn = dec->cn;
multi = dec->multi = MAX_MULTIMODAL;

//set_up_prd(dec);


printf("Strating Decode Image\n");
setbuf(stderr,NULL);


/////////////////////////
///////メモリ確保////////
/////////////////////////
  freq_array_save = (uint **)alloc_2d_array(MAX_MULTIMODAL ,(2 * 511 + 1) , sizeof(uint));//一次元配列
  freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
  cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列 
  sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));

	area1 = (int *)alloc_mem(AREA * sizeof(int  ));//一次元配列メモリ確保.初期化しようとしたらエラーが出た!
	area_o = (int *)alloc_mem(AREA * sizeof(int ));//一次元配列メモリ確保;


int **cost_save;
  cost_save = (int **)alloc_2d_array(  dec->height+1 , dec->width ,sizeof(int));
	cost_save[dec->height][0] = 0;

  array = (int *)alloc_mem(sizeof(int) * MAX_DATA_SAVE);
  img = alloc_image(dec->width, dec->height, dec->maxval);



	a = 1.0 / log(2.0);

double *wt_p;
int *org_p;
int *roff_p;
int area_cost_qt,count;
double a_cost,area_cost;
	wt_p = dec->ctx_weight_double;
int prd_buf,i;

dec->prd = 0;
prd_buf = 0;
int prd_max = dec->maxval << 6;




//puts("test5");






 for (y = 0; y < dec->height; y++) {
    for (x = 0; x < dec->width; x++) {


      TempleteM(img , dec, x , y , array,area1,area_o);

//特徴量Uの算出
      roff_p = dec->roff[y][x];
    	org_p = &dec->org[y][x];
      prd_buf = dec->prd;
      prd = 0;
      design_predictor_for_temp(dec, x, y, u, 1);
      for (i = 0; i < PRED_AREA; i++) {
         prd += org_p[roff_p[i]] * dec->temp_predictor[i];
      }
      if(prd < 0 || prd_max < prd ){prd = prd_buf;}
      dec->prd = prd;
 
      if(y == 0 && x == 0){
      u = 0;
      }else{
		  roff_p = dec->roff[y][x];
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
      th_p = dec->threshold;
      for (gr = 0; gr < dec->num_group; gr++) {
        if (u < th_p[gr]) break;
      }

      cn = dec->TM_cn[gr];

      pm = dec->pmodels[gr][cn];// + (base & mask);
      pm_save = dec->pmodels_save[gr][cn];// + (base & mask);
   
  
        if(x == 0 && y == 0){


//puts("test6");






      flg = 0;
      Shift_freq(img ,pm  ,dec,x,y,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,flg);}
        else{
      flg = 1;
      Shift_freq(img ,pm  ,dec,x,y,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,flg);
        }
 
//puts("test7");




//画素適応予測との併用
/*
      th_p = dec->prd_threshold;
      for (gr = 0; gr < dec->num_group; gr++) {
        if (u < th_p[gr]) break;
      }
      cn = dec->PRD_cn;
			pm_s = dec->pmodels_s[gr][cn];

      composite_MMF(dec,pm_s, x, y,u,pm->freq,pm->cumfreq,dec->limit_PRD);
*/
//復号		



      base = 255 ;
      p = rc_decode(fp, dec->rc, pm, base, base+dec->maxval+1)
      - base  ;


//どこまでデコードできたか表示する
//printf("%d,%d,%d\n",x,y,p);


 
//puts("test7-1");


//画素毎の見積もり符号量を算出
			cumbase = pm->cumfreq[base];//baseまでの高さの和
      cost = (float)(-a * log(pm->freq[base + p]));
      c = pm->cumfreq[base + dec->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));




 
//puts("test7-2");




      
      cost_save[y][x] = (int)((cost + subcost) * COST_ACCURACY);

      img->val[y][x] = p;



 
//puts("test7-3");



      for(j = 0;j < pm->size + 1; j++){
        pm->freq[j] = pm_save->freq[j] ; //元に戻しておく
        pm->cumfreq[j] = pm_save->cumfreq[j];

      } 





 
//puts("test7-4");







//puts("test8");



    }


//puts("test9");


  }


puts("test10");

return (img);
}





void write_pgm(IMAGE *img, char *filename)
{
	int i, j;
	FILE *fp;
	fp = fileopen(filename, "wb");
	fprintf(fp, "P5\n%d %d\n%d\n", img->width, img->height, img->maxval);
	for (i = 0; i < img->height; i++) {
		for (j = 0; j < img->width; j++) {
			putc(img->val[i][j], fp);
		}
	}
	fclose(fp);
	return;
}

union floatUnion
{
  float f;
  char bytes[4];
};

float getfloat(FILE *fp){
int i;
float n = 0;
union floatUnion fu;

  for(i = 0; i < 4; i++){
    fu.bytes[i] = (char)getbits(fp,8);
    printf("%d\n",fu.bytes[i]);
  }
  n = fu.f;

return(n);
}

void decode_prd_factor(FILE *fp,DECODER *dec){

int k;
double *prd_coef;

  dec->coef_start= getbits(fp,9);
  dec->coef_end = getbits(fp,9);
  dec->limit_PRD = getbits(fp,9);
  dec->trendline_a = getfloat(fp);
  dec->trendline_b = getfloat(fp);

int num_data = dec->coef_end - dec->coef_start;
  prd_coef = (double *)alloc_mem( (num_data + 1 ) * sizeof(double )); 

  for(k = 0; k < num_data; k++){
    prd_coef[k] = dec->trendline_a * (double)k + dec->trendline_b;
  }
  dec->prd_coef = prd_coef;

return;
}




void decode_cn(FILE *fp, DECODER *dec){

int i;
  for(i = 0; i < dec->num_group ; i++){
    dec->TM_cn[i] = getbits(fp,4);
    printf("TM_cn[%d]= %d\n",i,dec->TM_cn[i]);
  }
return;
}
/*
void decode_cn_prd(FILE *fp, DECODER *dec){

int i;
  for(i = 0; i < dec->num_group ; i++){
    dec->PRD_cn[i] = getbits(fp,4);
    printf("PRD_cn[%d]= %d\n",i,dec->PRD_cn[i]);
  }
return;
}
*/


int main(int argc, char **argv)
{
	int i;
	IMAGE *img;
	DECODER *dec;
	char *infile, *outfile;
	FILE *fp;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;
	for (i = 1; i < argc; i++) {
		if (infile == NULL) {
			infile = argv[i];
		} else {
			outfile = argv[i];
		}
	}
	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", 0.01 * VERSION);
		printf("usage: decmrp infile outfile\n");
		printf("infile:     Input file\n");
		printf("outfile:    Output file\n");
		exit(0);
	}
	fp = fileopen(infile, "rb");
	dec = init_decoder(fp);


puts("test1");



	dec->rc = rc_init();


//	rc_startdec(fp, dec->rc);
/*
	decode_class(fp, dec);
	decode_predictor(fp, dec);
	decode_threshold(fp, dec);
*/
	dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel,
		dec->pm_accuracy, NULL/*dec->pm_idx*/, dec->sigma,
		dec->maxval + 1);

	dec->pmodels_save = init_pmodels(dec->num_group, dec->num_pmodel,
		dec->pm_accuracy, NULL/*dec->pm_idx*/, dec->sigma,
		dec->maxval + 1);

	dec->pmodels_s = init_pmodels(dec->num_group, dec->num_pmodel,
		dec->pm_accuracy, NULL/*dec->pm_idx*/, dec->sigma,
		dec->maxval + 1);


// decode_prd_factor(fp, dec);
 decode_cn(fp, dec);


puts("test2");







// dec->PRD_cn = getbits(fp,4);
 dec->w_gr = getbits(fp,4); 




puts("test3");





printf("w_gr %d\n",dec->w_gr); 
printf("start %d\n",dec->coef_start);
printf("end %d\n",dec->coef_end);
printf("limit_PRD %d\n",dec->limit_PRD);
printf("trendline_a %f\n",dec->trendline_a);
printf("trendline_b %f\n",dec->trendline_b);
 
	rc_startdec(fp, dec->rc);







	decode_threshold(fp, dec);
//	decode_prd_threshold(fp, dec);





//puts("test4");




	img = decode_image(fp, dec);


//puts("test11");




	fclose(fp);



	write_pgm(img, outfile);
	printf("cpu time: %.2f sec.\n", cpu_time());


#if LOG_PUT_OUT_DEC
#if defined(_WIN32)
	if (set_directory()) {
		fprintf(stderr, "check \"DIR\"!\n");
		exit(1);
	}
#endif
	print_predictor(dec->predictor, dec->max_prd_order, dec->num_class, dec->max_coef, outfile);
	print_threshold(dec->th, dec->num_group, dec->num_class, NULL, dec->pm_idx, outfile);
//	print_class(dec->class, dec->num_class, dec->height, dec->width, outfile);
	print_class_color(dec->class, dec->num_class, dec->height, dec->width, outfile);
	print_class_and_block(dec->class, dec->num_class, dec->qtmap, dec->quadtree_depth, dec->height, dec->width, outfile);
	print_amp_chara(dec->predictor, dec->max_prd_order, dec->num_class, dec->height, dec->width, outfile);
#endif
	return(0);
}
