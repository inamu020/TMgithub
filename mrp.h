/****** LOG ********************************/
#define LOG_LIST_MODE		1
#define LOG_PUT_OUT_ENC		1
#define LOG_PUT_OUT_DEC		1

#if defined(_WIN32)		// Windows OS ?
# define BOUNDARY		"\\"
# define BS				'\\'
# define DIR			"C:\\test\\"	// Directory for output (Set this value according to your environment.)
															// [ATTENTION] This program does not support relative path!
#else
# define BOUNDARY		"/"
# define BS				'/'
# define DIR			"./Info/" // Directory for output (Must not change!)
#endif

// For Encoder only
#define LOG_LIST		DIR"Log_List.csv"
#define LOG_RATE_DIR	DIR"Rate_Map"BOUNDARY
#define LOG_VARI_DIR	DIR"Var_Upara"BOUNDARY

// For Encoder and Decoder
#define LOG_CL_DIR		DIR"Class_Map"BOUNDARY
#define LOG_VBS_DIR		DIR"Block_Size"BOUNDARY
#define LOG_TH_DIR		DIR"Threshold"BOUNDARY
#define LOG_PRED_DIR	DIR"Predictor"BOUNDARY
#define LOG_AMP_CH_DIR	DIR"Amp_Chara"BOUNDARY
#define LOG_TEMP_DIR  DIR"Temp_Map"BOUNDARY
#define LOG_MULTI_DIR  DIR"Multi_Jen"BOUNDARY
#define LOG_ETC_DIR  DIR"etc"BOUNDARY

/****** MRP-VERSION ************************/
#define MAGIC_NUMBER    ('M' << 8) + 'R'
#define BANNER          "ENCMRP/DECMRP version %.2f (Mar. 2010)"
#define VERSION         610

/****** OPTIMIZE ***************************/
#define OPT_SIDEINFO	1	// 1 : side-info into consideration (standard), 0 : neglect side-info
#define MAX_ITERATION   100
#define EXTRA_ITERATION 10
#define AUTO_DEL_CL		0
#define AUTO_PRD_ORDER	0

/***** BLOCK_SIZE **************************/
#define BASE_BSIZE      8
#define QUADTREE_DEPTH	4
#define MAX_BSIZE       32
#define MIN_BSIZE       (MAX_BSIZE >> QUADTREE_DEPTH)

/***** CLASS *******************************/
#define NUM_CLASS       -1
#define MAX_CLASS		63

/***** PREDICTOR ***************************/
#define COEF_PRECISION  6
#define PRD_ORDER       -1
#define BASE_PRD_ORDER  20
#define MAX_PRD_ORDER   110

/***** GROUP *******************************/
#define NUM_GROUP       16

/***** UPARA *******************************/
#define MAX_UPARA       512
#define UPEL_DIST       3
#define NUM_UPELS       (UPEL_DIST * (UPEL_DIST + 1))
#define CTX_WEIGHT		1	// 1 : weight on (standard) , 0 : weight off
#if CTX_WEIGHT
#	define MHD_WEIGHT	1	// 1 : Manhattan Distance (standard), 0 : Euclid Distance
#endif

/***** PMODEL ******************************/
#define PM_ACCURACY     3
#define NUM_PMODEL      16
#define MIN_FREQ        1
#define PMCLASS_MAX		16
#define PMCLASS_LEVEL	32
#define NUM_ZMODEL      49
#define TOT_ZEROFR      (1 << 10)
#define MAX_SYMBOL		1024	// must be >> MAX_UPARA

/***** TEMPLATE MATCHING *******************/
#define AREA						6
#define Y_SIZE          20  //54 * 27 when 512*512
#define X_SIZE          20
#define MAX_DATA_SAVE   200//64 * 4 //struct_num = 4. so data_num is 25. 
#define MAX_MULTIMODAL  45
#define BASE_MULTIMODAL 11
//#define MULTIMADAL_ON     //コメントを外した場合山の自動選択を行う.時間がかかるので注意されたし
#define REPETITION      1  //gr,cnの修正反復を何回やるか
#define VAR_ADJUSTMENT  1   //重みに正規分布を使用した時の調整.1で調整なし
#define WEIGHT_CN       5   //重み付けのための確立モデルのcn
#define TH_SET          64  //thを全選択する場合何個から選ぶか
#define TH_BITS         6   //thの情報を何ビット分送るか.
#define VARIABLE_WEIGHT
//#define RECIPROCAL_WEIGHT
#define W_SIZE          256
#define U_AREA          12
#define COST_ACCURACY   100
#define MAX_DATA_SAVE_DOUBLE 50
#define NAS_ACCURACY    100
#define PRED_AREA       6
#define TRAIN_AREA      111
#define PRED_ON         0

/***** RangeCoder **************************/
#define HAVE_64BIT_INTEGER	1
#if HAVE_64BIT_INTEGER
#  define RANGE_SIZE 64
#  if defined(__INTEL_COMPILER) || defined(_MSC_VER) || defined(__BORLANDC__)
#    define range_t unsigned __int64
#  else
#    define range_t unsigned long long
#  endif
#  define MAX_TOTFREQ (1 << 20)	/* must be < RANGE_BOT */
#else
#  define RANGE_SIZE 32
#  define range_t unsigned int
#  define MAX_TOTFREQ (1 << 14)	/* must be < RANGE_BOT */
#endif
#define RANGE_TOP  ((range_t)1 << (RANGE_SIZE - 8))
#define RANGE_BOT  ((range_t)1 << (RANGE_SIZE - 16))

/***** TYPE DEFINE *************************/
#define uint            unsigned int
#define img_t           unsigned char
#define cost_t          double

/***** PI **********************************/
#ifndef M_PI
#	define M_PI 3.14159265358979323846
#endif

/***** MACRO DEFINE ************************/
#define CLIP(min, max, i)       ((i < min) ? min : ((i > max) ? max : i))

/***** TIME ********************************/
#define HAVE_CLOCK

/***** STRUCTURE ***************************/
typedef struct {
	int height;
	int width;
	int maxval;
	img_t **val;
} IMAGE;

typedef struct {
	int size;
	int id;
	uint *freq;
	uint *cumfreq;
	float *cost;
	float *subcost;
	double norm;
} PMODEL;

typedef struct {
	range_t low;
	range_t code;
	range_t range;
} RANGECODER;

typedef struct {
	int y, x;
} CPOINT;

typedef struct {

int x;
int y;

}POINT;

typedef struct {

int bx;
int by;
double cor;
double var;

}TempleteM_S;



typedef struct {
	int height;
	int width;
	int maxval;
	int num_class;
	int num_group;
	int prd_order;
	int max_prd_order;
	int coef_precision;
	int num_pmodel;
	int pm_accuracy;
	int maxprd;
	int max_coef;
	int quadtree_depth;
	int optimize_loop;
	int **predictor;
	int **nzconv;
	int *num_nzcoef;
	int **th;
	int **upara;
	int **prd;
	int **encval;
	int **err;
	int **org;
	int *ctx_weight;
	int ***roff;
	int qtctx[QUADTREE_DEPTH << 3];
	char **qtmap[QUADTREE_DEPTH];
	char **class;
	char **group;
	char **uquant;
	int **econv;
	img_t *bconv;
	img_t *fconv;
	PMODEL ***pmodels;
	PMODEL **pmlist;
	PMODEL spm;
	RANGECODER *rc;
	double *sigma;
	int *mtfbuf;
	int *cl_hist;
	uint **cl_pos;
	int *coef_m;
#if AUTO_PRD_ORDER
	int prd_mhd;
	int *ord2mhd;
	int *num_search;
	int *zero_m;
	int *zero_fr;
	cost_t ***coef_cost;
#else
	cost_t **coef_cost;
#endif
	cost_t *th_cost;
	cost_t *class_cost;
	cost_t qtflag_cost[QUADTREE_DEPTH << 3];
#if AUTO_DEL_CL
	cost_t **err_cost;
#endif
   //using for TM
  int gr;
  int cn;
  int th_gr;
  int th_cn;
  int multi;
  int bit_count;
  int tm_th;
  int tm_sum;
  int tm_average;
  double tm_cost;
	PMODEL ***w_pmodels;
  POINT *point_th;
  int th_number;
  TempleteM_S ***tm;
  double **cost_save;
  int **temp_num;
  double *ctx_weight_double;
  int **TM_U;
  int *threshold;
  int *prd_threshold;
  int ***array;
  int w_gr;
  int ***w_pm_accuracy;
  int ***temp_predictor;
  int picsel_th;
  int *TM_cn;
  int **TM_gr;
  int **PRD_gr;
  int PRD_cn;
  double *prd_coef;
  uint coef_start;
  uint coef_end;
  uint limit_PRD;
  float trendline_a;
  float trendline_b;
  int **prd_b;
  int **comp;
} ENCODER;

typedef struct {
	int version;
	int height;
	int width;
	int maxval;
	int num_class;
	int num_group;
	int max_prd_order;
	int num_pmodel;
	int pm_accuracy;
	int maxprd;
	int max_coef;
	int coef_precision;
	int quadtree_depth;
	int **predictor;
	int **nzconv;
	int *num_nzcoef;
	int **th;
	int **err;
	int *ctx_weight;
	char **qtmap[QUADTREE_DEPTH];
	char **class;
	int *pm_idx;
	PMODEL ***pmodels;
	PMODEL spm;
	RANGECODER *rc;
	double *sigma;
	int *mtfbuf;
	int ***roff;
	int **org;
#if AUTO_PRD_ORDER
	int prd_mhd;
	int *zero_fr;
	int *ord2mhd;
#endif
  int multi;
  int cn;
  int th_gr;
  int th_cn;
	PMODEL ***w_pmodels;
  TempleteM_S ***tm;
  double **cost_save;
  int **temp_num;
	PMODEL ***pmodels_save;
	PMODEL ***pmodels_s;
  int *threshold;
  double *ctx_weight_double;
  int *array;
  int w_gr;
  int *w_pm_accuracy;
  int picsel_th;
	img_t *bconv;
  int *temp_predictor;
  int *TM_cn;
  int *prd_threshold;
  int PRD_cn;
  double *prd_coef;
  uint coef_start;
  uint coef_end;
  uint limit_PRD;
  float trendline_a;
  float trendline_b;
  int prd;
} DECODER;

/***** FUNC - common.c ****************************/
FILE *fileopen(char *, char *);
void *alloc_mem(size_t);
void **alloc_2d_array(int, int, int);
void  ***alloc_3d_array(int, int, int,int);
IMAGE *alloc_image(int, int, int);
PMODEL ***init_pmodels(int, int, int, int *, double *, int);
void set_spmodel(PMODEL *, int, int);
int *init_ctx_weight(void);
double *init_ctx_weight_double(void);
void mtf_classlabel(char **, int *, int, int, int, int, int);
double cpu_time(void);

/***** FUNC - rc.c ********************************/
RANGECODER *rc_init(void);
void rc_encode(FILE *, RANGECODER *, uint, uint, uint,int);
void rc_finishenc(FILE *, RANGECODER *,int);
int rc_decode(FILE *, RANGECODER *, PMODEL *, int, int);
void rc_startdec(FILE *, RANGECODER *);

/***** FUNC - log.c *******************************/
void print_predictor(int **, int, int, int, char *);
void print_threshold(int **, int, int, PMODEL **, int *, char *);
void print_class(char **, int, int, int, char *);
void print_class_color(char **, int, int, int, char *);
void print_class_and_block(char **, int, char ***, int, int, int, char *);
void print_amp_chara(int **, int, int, int, int, char *);
void print_etc(ENCODER *,int ***, char *);
//Lower funcs are can use in encoder only.
int elect_nearest_var_log(int);
void print_multimodal(ENCODER *,int ***, char *);
void Shift_freq_map(PMODEL *, int, int, ENCODER *,int ***,uint **,uint *,uint *,unsigned long long int *,int,int,int,int);
void print_rate_map(ENCODER *,int ***, char * );
void print_temp_map(ENCODER *,int ***, char * );
void print_block_size(int **, char ***, int, int, int, char *);
void calc_var_upara( ENCODER *, char *);
void init_log_sheet(ENCODER *, char *);
void finish_log_sheet(ENCODER *, int, int, int, int, int, int, double);


#if defined(_WIN32)
	int set_directory(void);
#endif

