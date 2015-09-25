#if defined(_WIN32)
# include <windows.h>	// For "make_direcotry"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mrp.h"

extern CPOINT dyx[];

char pallet[256][3] = {
{255, 128, 128}, {255, 255, 128}, {128, 255, 128}, {0, 255, 128},   {128, 255, 255},
{0, 128, 255},   {255, 128, 192}, {255, 128, 255}, {255, 0, 0},     {255, 255, 0},
{128, 255, 0},   {0, 255, 64},    {0, 255, 255},   {0, 128, 192},   {192, 128, 0},
{255, 0, 255},   {128, 64, 64},   {255, 128, 64},  {0, 255, 0},     {0, 128, 128},
{0, 64, 128},    {128, 128, 255}, {128, 0, 64},    {255, 0, 128},   {128, 0, 0},
{255, 128, 0},   {0, 128, 0},     {0, 128, 64},    {0, 0, 255},     {0, 0, 160},
{128, 0, 128},   {128, 0, 255},   {32, 32, 32},    {33, 33, 33},    {34, 34, 34},
{35, 35, 35},    {36, 36, 36},    {37, 37, 37},    {38, 38, 38},    {39, 39, 39},
{40, 40, 40},    {41, 41, 41},    {42, 42, 42},    {43, 43, 43},    {44, 44, 44},
{45, 45, 45},    {46, 46, 46},    {47, 47, 47},    {48, 48, 48},    {49, 49, 49},
{50, 50, 50},    {51, 51, 51},    {52, 52, 52},    {53, 53, 53},    {54, 54, 54},
{55, 55, 55},    {56, 56, 56},    {57, 57, 57},    {58, 58, 58},    {59, 59, 59},
{60, 60, 60},    {61, 61, 61},    {62, 62, 62},    {63, 63, 63},    {64, 64, 64},
{65, 65, 65},    {66, 66, 66},    {67, 67, 67},    {68, 68, 68},    {69, 69, 69},
{70, 70, 70},    {71, 71, 71},    {72, 72, 72},    {73, 73, 73},    {74, 74, 74},
{75, 75, 75},    {76, 76, 76},    {77, 77, 77},    {78, 78, 78},    {79, 79, 79},
{80, 80, 80},    {81, 81, 81},    {82, 82, 82},    {83, 83, 83},    {84, 84, 84},
{85, 85, 85},    {86, 86, 86},    {87, 87, 87},    {88, 88, 88},    {89, 89, 89},
{90, 90, 90},    {91, 91, 91},    {92, 92, 92},    {93, 93, 93},    {94, 94, 94},
{95, 95, 95},    {96, 96, 96},    {97, 97, 97},    {98, 98, 98},    {99, 99, 99},
{100, 100, 100}, {101, 101, 101}, {102, 102, 102}, {103, 103, 103}, {104, 104, 104},
{105, 105, 105}, {106, 106, 106}, {107, 107, 107}, {108, 108, 108}, {109, 109, 109},
{110, 110, 110}, {111, 111, 111}, {112, 112, 112}, {113, 113, 113}, {114, 114, 114},
{115, 115, 115}, {116, 116, 116}, {117, 117, 117}, {118, 118, 118}, {119, 119, 119},
{120, 120, 120}, {121, 121, 121}, {122, 122, 122}, {123, 123, 123}, {124, 124, 124},
{125, 125, 125}, {126, 126, 126}, {127, 127, 127}, {128, 128, 128}, {129, 129, 129},
{130, 130, 130}, {131, 131, 131}, {132, 132, 132}, {133, 133, 133}, {134, 134, 134},
{135, 135, 135}, {136, 136, 136}, {137, 137, 137}, {138, 138, 138}, {139, 139, 139},
{140, 140, 140}, {141, 141, 141}, {142, 142, 142}, {143, 143, 143}, {144, 144, 144},
{145, 145, 145}, {146, 146, 146}, {147, 147, 147}, {148, 148, 148}, {149, 149, 149},
{150, 150, 150}, {151, 151, 151}, {152, 152, 152}, {153, 153, 153}, {154, 154, 154},
{155, 155, 155}, {156, 156, 156}, {157, 157, 157}, {158, 158, 158}, {159, 159, 159},
{160, 160, 160}, {161, 161, 161}, {162, 162, 162}, {163, 163, 163}, {164, 164, 164},
{165, 165, 165}, {166, 166, 166}, {167, 167, 167}, {168, 168, 168}, {169, 169, 169},
{170, 170, 170}, {171, 171, 171}, {172, 172, 172}, {173, 173, 173}, {174, 174, 174},
{175, 175, 175}, {176, 176, 176}, {177, 177, 177}, {178, 178, 178}, {179, 179, 179},
{180, 180, 180}, {181, 181, 181}, {182, 182, 182}, {183, 183, 183}, {184, 184, 184},
{185, 185, 185}, {186, 186, 186}, {187, 187, 187}, {188, 188, 188}, {189, 189, 189},
{190, 190, 190}, {191, 191, 191}, {192, 192, 192}, {193, 193, 193}, {194, 194, 194},
{195, 195, 195}, {196, 196, 196}, {197, 197, 197}, {198, 198, 198}, {199, 199, 199},
{200, 200, 200}, {201, 201, 201}, {202, 202, 202}, {203, 203, 203}, {204, 204, 204},
{205, 205, 205}, {206, 206, 206}, {207, 207, 207}, {208, 208, 208}, {209, 209, 209},
{210, 210, 210}, {211, 211, 211}, {212, 212, 212}, {213, 213, 213}, {214, 214, 214},
{215, 215, 215}, {216, 216, 216}, {217, 217, 217}, {218, 218, 218}, {219, 219, 219},
{220, 220, 220}, {221, 221, 221}, {222, 222, 222}, {223, 223, 223}, {224, 224, 224},
{225, 225, 225}, {226, 226, 226}, {227, 227, 227}, {228, 228, 228}, {229, 229, 229},
{230, 230, 230}, {231, 231, 231}, {232, 232, 232}, {233, 233, 233}, {234, 234, 234},
{235, 235, 235}, {236, 236, 236}, {237, 237, 237}, {238, 238, 238}, {239, 239, 239},
{240, 240, 240}, {241, 241, 241}, {242, 242, 242}, {243, 243, 243}, {244, 244, 244},
{245, 245, 245}, {246, 246, 246}, {247, 247, 247}, {248, 248, 248}, {249, 249, 249},
{250, 250, 250}, {251, 251, 251}, {252, 252, 252}, {253, 253, 253}, {254, 254, 254},
{0, 0, 0}
};

/* Write Predictive Coef */
void print_predictor(int **predictor, int prd_order, int num_class, int max_coef, char *outfile)
{
	int i, j, cl, k, **p_prd;
	int max_mhd = 10;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_PRED_DIR"%s_prd.csv", name);
	fp = fileopen(file, "wb");
	p_prd = (int **)alloc_2d_array(max_mhd + 1, max_mhd << 1, sizeof(int));
	for (cl = 0; cl < num_class; cl++) {
		for (i = 0; i < max_mhd + 1; i++) {
			for (j = 0; j < max_mhd << 1; j++) {
				p_prd[i][j] = max_coef + 1;
			}
		}
		for (k = 0; k < prd_order; k++) {
			p_prd[max_mhd + dyx[k].y][max_mhd + dyx[k].x] = predictor[cl][k];
		}
		for (i = 0; i < max_mhd + 1; i++) {
			for (j = 0; j < max_mhd << 1; j++) {
				if (p_prd[i][j] == max_coef + 1) fprintf(fp, ",");
				else fprintf(fp, "%d,", p_prd[i][j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	free(p_prd);
	return;
}

/* Write Threshold */
void print_threshold(int **th, int num_group, int num_class, 
					 PMODEL **pmlist, int *pm_idx, char *outfile)
{
	int cl, gr;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_TH_DIR"%s_th.csv", name);
	fp = fileopen(file, "wb");
	for (cl = 0; cl < num_class; cl++) {
		for (gr = 0; gr < num_group - 1; gr++) {
			fprintf(fp, "%d,", th[cl][gr]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	if (pmlist != NULL) {
		for (gr = 0; gr < num_group; gr++) {
			fprintf(fp, "%d,", pmlist[gr]->id);
		}
		fprintf(fp, "\n");
	}
	if (pm_idx != NULL) {
		for (gr = 0; gr < num_group; gr++) {
			fprintf(fp, "%d,", pm_idx[gr]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

/* Draw Class Map */
void print_class(char **class, int num_class, int height, int width, char *outfile)
{
	int i, j, step;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_CL_DIR"%s_class.pgm", name);
	fp = fileopen(file, "wb");
	step = 255 / num_class;
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			putc(class[i][j] * step, fp);
		}
	}
	fclose(fp);
	return;
}

void print_class_color(char **class, int num_class, int height, int width, char *outfile)
{
	int y, x;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_CL_DIR"%s_class.ppm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", width, height);
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			putc(pallet[(int)(class[y][x])][0],fp);
			putc(pallet[(int)(class[y][x])][1],fp);
			putc(pallet[(int)(class[y][x])][2],fp);
		}
	}
	fclose(fp);
	return;
}

/* Draw Block Size */
void print_block_size(int **org, char ***qtmap, int quadtree_depth, 
					  int height, int width, char *outfile)
{
	int y, x, i, j, depth, bs;
	int **level;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr(outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_VBS_DIR"%s_blocksize.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	level = (int **)alloc_2d_array(height, width, sizeof(int));
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			level[y][x] = 0;
		}
	}
	if (quadtree_depth > 0) {
		for (depth = 0; depth < quadtree_depth; depth++) {
			for (y = 0; y < height; y++) {
				for (x = 0; x < width; x++) {
					j = y / (MAX_BSIZE >> depth);
					i = x / (MAX_BSIZE >> depth);
					if (qtmap[quadtree_depth - depth - 1][j][i] == 1) {
						level[y][x]++;
					}
				}
			}
		}
	}else {
		depth = MAX_BSIZE / BASE_BSIZE / 2;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				level[y][x] = depth;
			}
		}
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			depth = level[y][x];
			bs = MAX_BSIZE >> depth;
			if (y % bs == 0 || x %  bs == 0 || y == height - 1 || x == width - 1) {
				putc(255, fp);
			}else {
				putc(org[y][x], fp);
			}
		}
	}
	fclose(fp);
	free(level);
	return;
}

void print_class_and_block(char **class, int num_class, char ***qtmap, int quadtree_depth,
						   int height, int width, char *outfile)
{
	int y, x, i, j, depth, bs;
	int **level;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_CL_DIR"%s_class_block.ppm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", width, height);
	level = (int **)alloc_2d_array(height, width, sizeof(int));
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			level[y][x] = 0;
		}
	}
	if (quadtree_depth > 0) {
		for (depth = 0; depth < quadtree_depth; depth++) {
			for (y = 0; y < height; y++) {
				for (x = 0; x < width; x++) {
					j = y / (MAX_BSIZE >> depth);
					i = x / (MAX_BSIZE >> depth);
					if (qtmap[quadtree_depth - depth - 1][j][i] == 1) {
						level[y][x]++;
					}
				}
			}
		}
	}else {
		depth = MAX_BSIZE / BASE_BSIZE / 2;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				level[y][x] = depth;
			}
		}
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			depth = level[y][x];
			bs = MAX_BSIZE >> depth;
			if (y % bs == 0 || x %  bs == 0 || y == height - 1 || x == width - 1) {
				putc(255, fp);
				putc(255, fp);
				putc(255, fp);
			}else {
				putc(pallet[(int)(class[y][x])][0],fp);
				putc(pallet[(int)(class[y][x])][1],fp);
				putc(pallet[(int)(class[y][x])][2],fp);
			}
		}
	}
	fclose(fp);
	free(level);
	return;
}

/* Write Predictor Characteristic (Amplitude) */
void print_amp_chara(int **predictor, int prd_order, int num_class,
					 int height, int width, char *outfile)
{
    int fx, fy, j, cl, delta_y, delta_x;
    int *coef;
    double Re, Im, H;
    char *name;
    char file[256];
    FILE *fp;

    name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
    sprintf(file, LOG_AMP_CH_DIR"%s_amp_chara.csv", name);
    fp = fileopen(file, "wb");

    for(cl = 0; cl < num_class; cl++) {
        coef = &predictor[cl][0];
        fprintf(fp, "Predictor,%d\n\n", cl);
        fprintf(fp, " , , ");
        for(fx = 0; fx < (width / 2); fx += 5) {
            fprintf(fp, "%d,", fx);
        }
        fprintf(fp, "\n");
        for(fy = 0; fy < (height / 2); fy += 5) {
            fprintf(fp, " ,%d,", fy);
            for(fx = 0; fx < (width / 2); fx += 5) {
                Re = Im = 0.0;
                for(j = 0; j < prd_order; j++) {
                    delta_y = dyx[j].y;
                    delta_x = dyx[j].x;
                    Re += coef[j]
                      * cos(2 * M_PI * (((double)fx / width) * delta_x
                                        + ((double)fy / height) * delta_y));
                    Im += coef[j]
                      * sin(2 * M_PI * (((double)fx / width) *delta_x
                                        + ((double)fy / height) * delta_y));
                }
                H = sqrt(pow(Re, 2) + pow(Im, 2));
                fprintf(fp, "%f,", H);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n\n");
    }
    fclose(fp);
    return;
}






double sigma_b[] = {0.15, 0.26, 0.38, 0.57, 0.83, 1.18, 1.65, 2.31,
3.22, 4.47, 6.19, 8.55, 11.80, 16.27, 22.42, 30.89};




double continuous_GGF_map(ENCODER *enc, double e,int w_gr){

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

void Shift_freq_map(PMODEL *pm, int x , int y, ENCODER *enc,int ***a
,uint **freq_shift_array_save,uint *freq_shift_array,uint *cumfreq_shift_array,
unsigned long long int *sum_freq_shift_array,int multimodal,int w_gr,int w_func,int flg){

  int i,j,by[multimodal],bx[multimodal],e[multimodal];
  int flag2 = 0;
  double er_sum = 0.0,er_coef = 0.0;
  int w[multimodal];
  int multi = multimodal;
  
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
mc_w[i] = continuous_GGF_map(enc, mc_double,w_gr);
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


}








// Draw Predictive error rate 
void print_rate_map(ENCODER *enc,int ***array, char *outfile)
{
	int y, x, step,flg,cn,multi,c,count;
	int gr, base ;
	int max_rate = 20;
	double cost,a,subcost;
	char *name;
	char file[256];
	PMODEL *pm;
	FILE *fp;
  uint *freq_array,*cumfreq_array,**freq_array_save;

  freq_array_save = (uint **)alloc_2d_array(enc->multi ,(2 * 511 + 1) , sizeof(uint));//一次元配列
  freq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列
  cumfreq_array =  (uint *)alloc_mem((2 * 511 + 1) * sizeof(uint));//一次元配列

unsigned long long int *sum_freq_array;
sum_freq_array =  (unsigned long long int *)alloc_mem((512 * 2 + 1) * sizeof(unsigned long long int));


  int *th_p;
  int u,e,cumbase,w_gr = enc->w_gr,w_func = 1;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_RATE_DIR"%s_rate.pgm", name);
	fp = fileopen(file, "wb");
	step = 255 / max_rate;
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);

  gr = enc->gr;
  cn = enc->cn;
  multi = enc->multi;
  base = 255;
  pm = enc->pmodels[gr][cn];

	a = 1.0 / log(2.0);
  count = 0;


	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
      th_p = enc->threshold;


      u = enc->TM_U[y][x];
      gr = enc->TM_gr[y][x];
      cn = enc->TM_cn[gr];
			pm = enc->pmodels[gr][cn];

        if(x == 0 && y == 0){
      flg = 0;
      Shift_freq_map(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg);
      }else{
      flg = 1;
      Shift_freq_map(pm , x , y ,enc,array,freq_array_save,freq_array,cumfreq_array,sum_freq_array,multi,w_gr,w_func,flg);
        }

      e = enc->encval[y][x];//輝度値代入

      base = 255;
			cumbase = cumfreq_array[base];//baseまでの高さの和
      cost = (float)(-a * log(freq_array[base + e]));
      c = cumfreq_array[base + enc->maxval + 1] - cumbase;
      subcost = (float)(a * log(c));

			cost = (cost + subcost) * step;

			if(cost < 0) cost = 0;
			else if(cost > 255) cost = 255;
			putc((int)cost, fp);

		}//x fin
	}//y fin

	fclose(fp);
	return;
}


void print_temp_map(ENCODER *enc, int ***a, char *outfile)
{
	int i, j=0,g,h,x;
	char *name;
	char file[256];
	FILE *fp;
  int *array_save,temp = 0;
  double step;
  int count=0; 
 
  array_save =  (int *)alloc_mem((enc->height * enc->width * 2) * sizeof(int));//一次元配列
  for (i = 0; i < enc->height; i++) {
    for (j = 0; j < enc->width; j++) {
      array_save[i * enc->height + j] = a[i][j][3];
      count++;
    }
  }

  j = enc->height * enc->width;

  for (g = 0; g < j - 1; g++) {
        for (h = j - 1; h > g; h--) {
            if (array_save[h - 1] < array_save[h]) {  // 前の要素の方が小さかったら 
                temp = array_save[h];        // 交換する 
                array_save[h] = array_save[h - 1];
                array_save[h - 1] = temp;
            }
        } 
    }

  
  printf("%d\n",array_save[0]);  

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_TEMP_DIR"%s_temp.pgm", name);

	fp = fileopen(file, "wb");

	//step = (255 / (double)array_save[0]);

	step = (255 / (double)5000.0);


	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);
	for (i = 0; i < enc->height; i++) {
		for (j = 0; j < enc->width; j++) {
      x = a[i][j][3] * step;
			putc(x, fp);
		}
	}
	fclose(fp);
	return;
}




/* Write upara - variance of predictive error */
void calc_var_upara( ENCODER *enc, char *outfile)
{
    int y, x, u, prd, base, pels, upara;
    int maximum_upara = 0;
    double ave, var;
	char *name;
	char file[256];
    FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
    sprintf( file, LOG_VARI_DIR"%s_vari_upara.csv", name);
    fp = fileopen( file, "wb");

    for( y = 0; y < enc->height; y++) {
        for( x = 0; x < enc->width; x++) {
            if( enc->upara[y][x] > maximum_upara)
              maximum_upara = enc->upara[y][x];
        }
    }

	fprintf(fp, "Upara,Num_Pels,Variance,Deviation\n");
    for( u = 0; u <= maximum_upara; u++) {
        ave = 0;
        pels = 0;
        var = 0.0;
        for( y = 0; y < enc->height; y++) {
            for( x = 0; x < enc->width; x++) {
                upara = enc->upara[y][x];

                if( u != upara)
                  continue;

                prd = enc->prd[y][x];
                if( prd < 0) prd = 0;
                else if( prd > enc->maxprd) prd = enc->maxprd;
                base = enc->bconv[prd];

                ave += base + enc->org[y][x];
                pels++;
            }
        }
        if( pels <= 0) {
            ave = 0;
        }else {
            ave /= pels;
        }
        for( y = 0; y < enc->height; y++) {
            for( x = 0; x < enc->width; x++) {
                upara = enc->upara[y][x];

                if( u != enc->upara[y][x])
                  continue;

                prd = enc->prd[y][x];
                if( prd < 0) prd = 0;
                else if( prd > enc->maxprd) prd = enc->maxprd;
                base = enc->bconv[prd];

                var += pow( base + enc->org[y][x] - ave, 2.0);
            }
        }
        if( pels <= 0) {
            var = 0;
            fprintf( fp, "%d,%d, , \n", u, pels);
        }else {
            var /= pels;
            fprintf( fp, "%d,%d,%.4f,%.4f\n", u, pels, var, sqrt(var));
        }
    }
    fprintf( fp,"\n");
    return;
}

void init_log_sheet(ENCODER *enc, char *outfile)
{
	char *name;
	char date[64];
	time_t ct = time(NULL);
	struct tm *lst = localtime(&ct);
	FILE *fp;

	if ((name = strrchr(outfile, BS)) == NULL) {
		name = outfile;
	}else {
		name++;
	}
	
	sprintf(date, "%4d/%2d/%2d_%2d:%2d",
		1900 + lst->tm_year, 1 + lst->tm_mon, lst->tm_mday, lst->tm_hour, lst->tm_min);
    
    if ((fp = fopen(LOG_LIST, "rb")) == NULL) {
	fp = fileopen(LOG_LIST, "wb");
	fprintf(fp, "Date,Input,Height,Width,Init_Class,Use_Class,Prd_Order,\
Header[bits],Class[bits],Predictor[bits],Threshold[bits],\
Pred. Errors[bits],Total info.[bits],Coding Rate[bits/pel], ,");
	
	fprintf(fp, "Auto_Del_CL,\
Auto_Set_Coef,\
BASE_BSIZE,\
QUADTREE_DEPTH,MIN_BSIZE,MAX_BSIZE,\
COEF_PRECISION,PM_ACCURACY,NUM_GROUP,UPEL_DIST,\
CTX_WEIGHT,");
	
	fprintf(fp, "\n");
	fclose(fp);
    }
    
    fp = fileopen(LOG_LIST, "ab");
    fprintf(fp, "%s,%s,%d,%d,%d,", date, name, enc->height, enc->width, enc->num_class);
	fclose(fp);

	return;
}

void finish_log_sheet(ENCODER *enc, int header_info, int class_info, int pred_info,
							  int th_info, int err_info, int total_info, double rate)
{
	FILE *fp;

	fp = fileopen(LOG_LIST, "ab");
	fprintf(fp, "%d,%d,", enc->num_class, enc->max_prd_order);
	fprintf(fp, "%d,%d,%d,%d,%d,%d,", header_info, class_info, pred_info, th_info, err_info, total_info);
	fprintf(fp, "%f, ,", rate);

	if(AUTO_DEL_CL) {
		fprintf(fp, "ON,");
	}else {
		fprintf(fp, "OFF,");
	}

	if(AUTO_PRD_ORDER) {
		fprintf(fp, "ON,");
	}else {
		fprintf(fp, "OFF,");
	}
	
	fprintf(fp, "%d,", BASE_BSIZE);
	fprintf(fp, "%d,%d,%d,", enc->quadtree_depth, MIN_BSIZE, MAX_BSIZE);
	fprintf(fp, "%d,%d,%d,%d,", enc->coef_precision, enc->pm_accuracy,enc->num_group, UPEL_DIST);
	
	if(CTX_WEIGHT) {
		if(MHD_WEIGHT) {
			fprintf(fp, "mhd,");
		}else {
			fprintf(fp, "euclid,");
		}
	}else {
		fprintf(fp, "Non,");
	}

	fprintf(fp, "\n");
	fclose(fp);

	return;
}

#if defined(_WIN32)
int check_pathname(LPTSTR lpPathName)
{
	if (lpPathName[1] != ':' || lpPathName[2] != '\\') return 0;
	if (strpbrk(lpPathName + 2, TEXT("/,:;*\"<>|"))) return 0;
	if (strstr(lpPathName, TEXT("\\\\"))) return 0;
	if (strstr(lpPathName, TEXT(" \\"))) return 0;

	return 1;
}

int make_directory_sub(LPTSTR lpPathName)
{
	char lpPrePath[MAX_PATH], lpCorrPath[MAX_PATH];
	char *tmp;
	int flag;

	strcpy(lpCorrPath, lpPathName);
	tmp = strrchr(lpCorrPath, '\\');
	if (!tmp) return 0;
	tmp[0] = '\0';

	flag = CreateDirectory(lpCorrPath, NULL);
	if (!flag) {
		strcpy(lpPrePath, lpCorrPath);
		make_directory_sub(lpPrePath);
		flag = CreateDirectory(lpCorrPath, NULL);
	}
	return (flag);
}

int make_directory(LPTSTR lpPathName)
{
	int flag = 0;

	if (check_pathname(lpPathName)) {
		CreateDirectory(lpPathName, NULL);
		if (ERROR_ALREADY_EXISTS == GetLastError()){
			flag = 1;
		}else {
			flag = make_directory_sub(lpPathName);
		}
	}
	return (flag);
}

int set_directory(void)
{
	int flag = 1;

	flag *= make_directory(LOG_RATE_DIR);
	flag *= make_directory(LOG_VARI_DIR);
	flag *= make_directory(LOG_CL_DIR);
	flag *= make_directory(LOG_VBS_DIR);
	flag *= make_directory(LOG_TH_DIR);
	flag *= make_directory(LOG_PRED_DIR);
	flag *= make_directory(LOG_AMP_CH_DIR);
 
	if (flag == 0) {
		fprintf(stderr, "Can't make directory!\n");
		return (1);
	}
	return (0);
}
#endif
