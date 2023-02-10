#define __STDC_FORMAT_MACROS 1
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <limits.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include "izlib.h"
#include "ksw.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAXBQ 64
#define VERSION "r16"
#define abs(x) ((x)>0?(x):-(x))
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)>(y)?(y):(x))
#define basename(str) (strrchr(str, '/') ? strrchr(str, '/') + 1 : str)
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

/***************
 * CMD options *
 ***************/

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

typedef struct {
	int type, len;
	uint8_t *seq;
	uint64_t cnt;
} ta_adap_t;

typedef struct {
	uint64_t reads, *bases, *quals, q20, q30;
} qc_sta_t;

typedef struct {
	int sa, sb, go, ge; // scoring; not exposed to the command line, for now
	int min_sc, min_len, min_trim_len, min_rlen;
	unsigned kn, sl; // selected length to output
	char mskr, qc;
	int n_threads;
	int chunk_size;
	double max_diff;
	int n_adaps, m_adaps;
	ta_adap_t *adaps;
	qc_sta_t *qcstat;
	int8_t mat[25];
	kseq_t *ks;
} ta_opt_t;

void ta_opt_set_mat(int sa, int sb, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
}

void ta_opt_init(ta_opt_t *opt)
{
	memset(opt, 0, sizeof(ta_opt_t));
	opt->sa = 1, opt->sb = 2, opt->go = 1, opt->ge = 3;
	opt->min_sc = 15, opt->min_len = 8, opt->min_trim_len = 1000000, opt->min_rlen = 35, opt->max_diff = .15f;
	opt->sl = UINT_MAX;
	opt->mskr = 'X';
	opt->qc = 0;
	opt->qcstat = calloc(1, sizeof(qc_sta_t));
	opt->qcstat->reads = opt->qcstat->q20 = opt->qcstat->q30 = 0;
	opt->qcstat->bases = calloc(5, sizeof(uint64_t));
	opt->qcstat->quals = calloc(MAXBQ, sizeof(uint64_t));
	opt->n_threads = 1, opt->chunk_size = 10000000;
	ta_opt_set_mat(opt->sa, opt->sb, opt->mat);
}

void ta_opt_add_adap(ta_opt_t *opt, int type, const char *adap)
{
	ta_adap_t *p;
	uint8_t *q;
	if (opt->m_adaps == opt->n_adaps) {
		opt->m_adaps = opt->m_adaps? opt->m_adaps<<1 : 4;
		opt->adaps = realloc(opt->adaps, opt->m_adaps * sizeof(ta_adap_t));
	}
	p = &opt->adaps[opt->n_adaps++];
	p->len = strlen(adap);
	assert(p->len * opt->sa < 256);
	p->seq = (uint8_t*)strdup(adap);
	p->type = type;
	p->cnt = 0;
	for (q = p->seq; *q; ++q)
		*q = seq_nt4_table[*q];
}

void ta_opt_default_adaps(ta_opt_t *opt)
{
	if (opt->n_adaps) return;
	// Illumina
	ta_opt_add_adap(opt, 5, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
	ta_opt_add_adap(opt, 3, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
	ta_opt_add_adap(opt, 3, "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT");
	ta_opt_add_adap(opt, 3, "ATCTCGTATGCCGTCTTCTGCTTG");
	// DNBSEQ
	ta_opt_add_adap(opt, 5, "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG");
	ta_opt_add_adap(opt, 3, "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA");
	ta_opt_add_adap(opt, 3, "GAACGACATGGCTACGATCCGACTT");
}

void ta_opt_open(ta_opt_t *opt, const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	opt->ks = kseq_init(fp);
}

void ta_opt_free(ta_opt_t *opt)
{
	int i;
	gzFile fp = opt->ks->f->f;
	for (i = 0; i < opt->n_adaps; ++i)
		free(opt->adaps[i].seq);
	kseq_destroy(opt->ks);
	gzclose(fp);
	free(opt->adaps);
	free(opt->qcstat->bases);
	free(opt->qcstat->quals);
	free(opt->qcstat);
}

/*************
 * FASTQ I/O *
 *************/

typedef struct {
	int l_seq;
	char *name, *seq, *qual, *comment;
} bseq1_t;

bseq1_t *bseq_read(kseq_t *ks, int chunk_size, int *n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->comment = ks->comment.l? strdup(ks->comment.s) : 0;
		s->l_seq = ks->seq.l;
		size += seqs[n++].l_seq;
		if (size >= chunk_size) break;
	}
	*n_ = n;
	return seqs;
}

/****************
 * Trim one seq *
 ****************/

void ta_trim1(ta_opt_t *opt, char *seq)
{
	int i, j, k;
	kstring_t _str = {0,0,0}, *str = &_str;
	str->l = strlen(seq);
	str->s = malloc(str->l);
	for (i = 0; i < str->l; ++i)
		str->s[i] = seq_nt4_table[(uint8_t)seq[i]];
	for (j = 0; j < opt->n_adaps; ++j) {
		kswr_t r;
		double diff;
		int type;
		ta_adap_t *p = &opt->adaps[j];
		r = ksw_align(p->len, p->seq, str->l, (uint8_t*)str->s, 5, opt->mat, opt->go, opt->ge, KSW_XBYTE|KSW_XSTART|(opt->min_len * opt->sa), 0);
		++r.te; ++r.qe; // change to 0-based
		k = r.qe - r.qb < r.te - r.tb? r.qe - r.qb : r.te - r.tb;
		diff = (double)(k * opt->sa - r.score) / opt->sb / k;
		//printf("%d:%.3f [%d,%d):%d <=> [%d,%d):%d\n", r.score, diff, r.qb, r.qe, p->len, r.tb, r.te, (int)str.l);
		if (r.qb <= r.tb && p->len - r.qe <= str->l - r.te) { // contained
			if (r.qb * opt->sa > opt->sa + opt->sb) continue;
			if ((p->len - r.qe) * opt->sa > opt->sa + opt->sb) continue;
			type = 1;
		} else if (r.qb <= r.tb) { // 3'-end overlap
			if (r.qb * opt->sa > opt->sa + opt->sb) continue;
			if ((str->l - r.te) * opt->sa > opt->sa + opt->sb) continue;
			type = 2;
		} else { // 5'-end overlap
			if ((p->len - r.qe) * opt->sa > opt->sa + opt->sb) continue;
			if (r.tb * opt->sa > opt->sa + opt->sb) continue;
			type = 3;
		}
		if (p->type == 5) {
			if (r.tb == 0 && r.qe == p->len && (r.te - r.tb) * opt->sa == r.score)
				type = 4;
		} else if (p->type == 3) {
			if (r.qb == 0 && r.te == str->l && (r.te - r.tb) * opt->sa == r.score)
				type = 4;
		}
		if (type == 4) { // exact match
			if (r.te - r.tb < opt->min_len) continue;
		} else { // inexact match
			if (r.score < opt->min_sc || diff > opt->max_diff) continue;
		}
		__sync_fetch_and_add(&p->cnt, 1);
		if (p->type == 5) {
			k = r.te + (p->len - r.qe);
			k = k < str->l? k : str->l;
			for (i = 0; i < k; ++i) seq[i] = opt->mskr;
		} else if (p->type == 3) {
			k = r.tb > r.qb? r.tb - r.qb : 0;
			for (i = k; i < str->l; ++i) seq[i] = opt->mskr;
		}
	}
	free(str->s);
}

static void apply_trim(int min_trim, int l_seq, char *seq, char *qual, char mskr)
{
	int i;
	while (l_seq > min_trim && seq[l_seq-1] == mskr) --l_seq;
	for (i = 0; l_seq - i > min_trim && seq[i] == mskr;) ++i;
	if (i > 0) {
		memmove(seq, seq + i, l_seq - i + 1);
		if (qual) memmove(qual, qual + i, l_seq - i + 1);
		l_seq -= i;
	}
	seq[l_seq] = 0;
	if (qual) qual[l_seq] = 0;
}

static int trim_len(const int l_seq, const char *seq, const char mskr)
{
	int i, n = 0;
	for (i = 0; i < l_seq; ++i)
		n += (seq[i] == mskr);
	return n;
}

static void do_qc(const bseq1_t *s, qc_sta_t *q)
{
	int i;
	++q->reads;
	for (i = 0; i < s->l_seq; ++i)
	{
		++q->bases[seq_nt16_int[seq_nt16_table[s->seq[i]]]];
		if (s->qual)
			++q->quals[s->qual[i] - 33];
	}
}

static void dump_qc(qc_sta_t *q)
{
	int i;
	for (i = 29; i < MAXBQ; ++i) q->q30 += q->quals[i];
	for (i = 19; i < 29; ++i) q->q20 += q->quals[i];
	q->q20 += q->q30;
	fputs("Reads\tBases\tQ20_bases\tQ30_bases\tGC_bases\n", stderr);
	fprintf(stderr, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",
		q->reads, q->bases[0] + q->bases[1] + q->bases[2] + q->bases[3] + q->bases[4],
		q->q20, q->q30, q->bases[1] + q->bases[2]);
}

/**********************
 * Callback functions *
 **********************/

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	ta_opt_t *opt;
} data_for_t;

static void worker_for(void *_data, long i, int tid)
{
	data_for_t *data = (data_for_t*)_data;
	ta_trim1(data->opt, data->seqs[i].seq);
}

static unsigned strn(const char *a)
{
	unsigned n = 0;
	const char *p = a;
	while (*p)
		n += (bool)(*p++ == 'N');
	return n;
}

static void *worker_pipeline(void *shared, int step, void *_data)
{
	int i;
	ta_opt_t *opt = (ta_opt_t*)shared;
	if (step == 0) {
		data_for_t *ret;
		ret = calloc(1, sizeof(data_for_t));
		ret->seqs = bseq_read(opt->ks, opt->chunk_size, &ret->n_seqs);
		ret->opt = opt;
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		data_for_t *data = (data_for_t*)_data;
		kt_for(opt->n_threads, worker_for, data, data->n_seqs);
		return data;
	} else if (step == 2) {
		data_for_t *data = (data_for_t*)_data;
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			// filter out reads shorter than threshold specified
			if (s->l_seq - trim_len(s->l_seq, s->seq, opt->mskr) < opt->min_rlen)
				goto next;
			if (opt->min_trim_len < s->l_seq)
				apply_trim(opt->min_trim_len, s->l_seq, s->seq, s->qual, opt->mskr);
			if (opt->kn && strn(s->seq) >= opt->kn)
				goto next;
			putchar(s->qual? '@' : '>'); fputs(s->name, stdout);
			if (s->comment) {
				putchar(' ');
				puts(s->comment);
			} else putchar('\n');
			if (s->l_seq > opt->sl) // select first n bases
			{
				s->seq[opt->sl] = '\0';
				if (s->qual)
					s->qual[opt->sl] = '\0';
				s->l_seq = opt->sl;
			}
			puts(s->seq);
			if (s->qual) {
				puts("+"); puts(s->qual);
			}
			// basic qc
			if (opt->qc) do_qc(s, opt->qcstat);
next:
			free(s->seq); free(s->qual); free(s->name); free(s->comment);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

/*****************
 * Main function *
 *****************/

int main(int argc, char *argv[])
{
	int c, i, j;
	ta_opt_t opt;
	ta_opt_init(&opt);
	while ((c = getopt(argc, argv, "5:3:s:p:l:t:r:m:n:k:qvh")) >= 0) {
		if (c == 'h') goto usage;
		else if (c == '5' || c == '3') ta_opt_add_adap(&opt, c - '0', optarg);
		else if (c == 's') opt.min_sc = atoi(optarg);
		else if (c == 'd') opt.max_diff = atof(optarg);
		else if (c == 'l') opt.min_len = atoi(optarg);
		else if (c == 'p') opt.n_threads = atoi(optarg);
		else if (c == 't') opt.min_trim_len = atoi(optarg);
		else if (c == 'r') opt.min_rlen = atoi(optarg);
		else if (c == 'm') opt.mskr = *optarg;
		else if (c == 'n') opt.sl = atoi(optarg);
		else if (c == 'k') opt.kn = atoi(optarg);
		else if (c == 'q') opt.qc = 1;
		else if (c == 'v') {
			puts(VERSION);
			return 0;
		}
		else
		{
			puts("Invalid option");
			return 0;
		}
	}
	if (!opt.sl)
	{
		fprintf(stderr, "Error: Invalid sequence length via -n: [%d]\n", opt.sl);
		exit(1);
	}
	if (opt.mskr != 'X' && opt.mskr != 'N')
	{
		fprintf(stderr, "Error: Invalid masker character: [%c]\n", opt.mskr);
		exit(1);
	}
	if (opt.n_adaps == 0) ta_opt_default_adaps(&opt);
	if (optind == argc && isatty(fileno(stdin))) {
usage:
		fprintf(stderr, "Usage: \033[1;31m%s\033[0;0m [options] <in.fq>\n", basename(argv[0]));
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -5 STR     5'-end adapter\n");
		fprintf(stderr, "  -3 STR     3'-end adapter\n");
		fprintf(stderr, "  -l INT     min length \033[2m[%d]\033[0m\n", opt.min_len);
		fprintf(stderr, "  -s INT     min score \033[2m[%d]\033[0m\n", opt.min_sc);
		fprintf(stderr, "  -t INT     trim down masked part (Xs) \033[2m[don't trim]\033[0m\n");
		fprintf(stderr, "  -d FLOAT   max difference \033[2m[%.3f]\033[0m\n", opt.max_diff);
		fprintf(stderr, "  -r INT     min read length (w/o trimmed) to output \033[2m[%d]\033[0m\n", opt.min_rlen);
		fprintf(stderr, "  -p INT     number of trimmer threads \033[2m[%d]\033[0m\n", opt.n_threads);
		fprintf(stderr, "  -m CHAR    masker character (X or N) \033[2m[%c]\033[0m\n", opt.mskr);
		fprintf(stderr, "  -n INT     first n bases to keep \033[2m[all]\033[0m\n");
		fprintf(stderr, "  -q         perform basic QC of trimmed reads \033[2m[false]\033[0m\n");
		fprintf(stderr, "             QC info goes to stderr\n");
		fprintf(stderr, "  -k INT     keep reads with <\033[3mk\033[0m non-adapter-Ns in output \033[2m[unlimited]\033[0m\n\n");
		fprintf(stderr, "  -h         print help message\n");
		fprintf(stderr, "  -v         print version number\n");
		fprintf(stderr, "NB: please redirect output fq (stdout) to file\n");
		return 1;
	}
	ta_opt_open(&opt, optind < argc ? argv[optind] : 0);
	kt_pipeline(2, worker_pipeline, &opt, 3);
	for (j = 0; j < opt.n_adaps; ++j) {
		ta_adap_t *p = &opt.adaps[j];
		fprintf(stderr, "%-15ld ", (long)p->cnt);
		for (i = 0; i < p->len; ++i) fputc("ACGTN"[(int)p->seq[i]], stderr);
		fputc('\n', stderr);
	}
	if (opt.qc) dump_qc(opt.qcstat);
	ta_opt_free(&opt);
	return 0;
}
