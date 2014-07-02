#include <sys/types.h>

#include <assert.h>
#include <stdio.h>

#include "ruby.h"
#include "nwam.h"

#define nelems(e)	(sizeof (e) / sizeof *(e))

#define	sizeof_dxr_vec(n)	(offsetof(struct xr_vec, elems[(n)]))

struct wamdata {
	WAM *w;
};

struct syminfo_conv_cookie {
	WAM *w;
	int dir;
};

/* UTILITIES (1) */
static void closed_wam(void);
static void free_wam(struct wamdata *);

/* ENTRY */
static VALUE fwam_close(VALUE);
static VALUE fwam_closed(VALUE);
static VALUE fwam_alloc(VALUE);
static VALUE fwam_initialize(int, VALUE *, VALUE);
static VALUE fwam_s_open(int, VALUE *, VALUE);
static VALUE fwam_size(VALUE, VALUE);
static VALUE fwam_elem_num(VALUE, VALUE, VALUE);
static VALUE fwam_freq_sum(VALUE, VALUE, VALUE);
static VALUE fwam_max_elem_num(VALUE, VALUE);
static VALUE fwam_max_freq_sum(VALUE, VALUE);
static VALUE fwam_total_elem_num(VALUE);
static VALUE fwam_total_freq_sum(VALUE);
static VALUE fwam_id2name(VALUE, VALUE, VALUE);
static VALUE fwam_name2id(VALUE, VALUE, VALUE);
static VALUE fwam_get_vec(VALUE, VALUE, VALUE);
static VALUE fwam_prop_nentries(VALUE, VALUE, VALUE);
static VALUE fwam_prop_gets(VALUE, VALUE, VALUE, VALUE);
static VALUE fwam_fss_gets(VALUE, VALUE, VALUE);
/* UTILITY */ static VALUE wsh0(VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE);
/* UTILITY */ static VALUE ncsb0(VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE);
static VALUE fwam_wsh(VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE);
static VALUE fwam_bex_wsh(VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE);
static VALUE fwam_fss_wsh(VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE);
static VALUE fwam_ncsb(VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE, VALUE);
static VALUE fwam_text2vec(VALUE, VALUE, VALUE);
static VALUE fwam_get_last_nd(VALUE);

/* UTILITIES (2) */
static VALUE vec2value(struct xr_vec const *, WAM *, int);
static int elem2value(void *, void *, VALUE *);
static int cc2value(void *, void *, VALUE *);
static int charpp2value(void *, void *, VALUE *);
static int syminfo2value(void *, void *, VALUE *);
static int value2syminfo_e_r(VALUE, void *, void *);
static struct bxue_t *value2bxue_t(VALUE, WAM *, int, ssize_t *);
static int value2bxue_t_r(VALUE, void *, void *);
static struct fss_query *value2fss_query(VALUE);
static int value2fss_con_query_r(VALUE, void *, void *);
static int value2fss_simple_query_r(VALUE, void *, void *);
static struct xr_vec *value2vec(VALUE, struct syminfo_conv_cookie *);

static VALUE array2value(void *, size_t, size_t, int (*)(void *, void *, VALUE *), void *);
static void *value2array(VALUE, size_t *, size_t, int (*)(VALUE, void *, void *), void (*)(void *, void *), void *);
static void free_vecs(struct syminfo *, size_t);
static void free_fq(struct fss_query *);
static size_t iwl(struct syminfo *, size_t, size_t, WAM *, int);

void Init_rnwam(void);

static VALUE cWAM;
static VALUE eWAMError;
static VALUE cWAM_XR_VEC;
static VALUE cWAM_XR_ELEM;
static VALUE cWAM_SYMINFO_ELEM;
static VALUE cWAM_BXU_ELEM;
static VALUE cWAM_SIMPLE_QUERY;

static ssize_t last_nd = 0;

#define	sym_compar	n_sym_compar

/**************************************************************************
 * UTILITIES (1)
 **************************************************************************/

static void
closed_wam()
{
	rb_raise(eWAMError, "close WAM file");
}

#define GetWAM(obj, wamp) do { \
	Data_Get_Struct(obj, struct wamdata, wamp); \
	if (wamp == 0 || wamp->w == 0) closed_wam(); \
} while (0)

#define	GetWAM2(obj, data, wam) do { \
	GetWAM(obj, data); \
	(wam) = wamp->w; \
} while (0)

static void
free_wam(wamp)
struct wamdata *wamp;
{
	if (wamp) {
		if (wamp->w) {
			wam_close(wamp->w);
		}
		free(wamp);
	}
}

/**************************************************************************
 * ENTRY
 **************************************************************************/

static VALUE
fwam_close(obj)
VALUE obj;
{
	struct wamdata *wamp;

	GetWAM(obj, wamp);
	wam_close(wamp->w);
	wamp->w = 0;

	return Qnil;
}

static VALUE
fwam_closed(obj)
VALUE obj;
{
	struct wamdata *wamp;

	Data_Get_Struct(obj, struct wamdata, wamp);
	if (wamp == 0 || wamp->w == 0) {
		return Qtrue;
	}
	return Qfalse;
}

static VALUE
fwam_alloc(klass)
VALUE klass;
{
	return Data_Wrap_Struct(klass, 0, free_wam, 0);
}

static VALUE
fwam_initialize(argc, argv, obj)
VALUE *argv, obj;
{
	VALUE file;
	struct wamdata *wamp;
	WAM *w;

	rb_scan_args(argc, argv, "1", &file);
	SafeStringValue(file);
	w = wam_open(RSTRING_PTR(file), NULL);
	if (!w) {
		return Qnil;
	}

	wamp = ALLOC(struct wamdata);
	DATA_PTR(obj) = wamp;
	wamp->w = w;
	return obj;
}

static VALUE
fwam_s_open(argc, argv, klass)
VALUE *argv, klass;
{
	VALUE obj = Data_Wrap_Struct(klass, 0, free_wam, 0);

	if (NIL_P(fwam_initialize(argc, argv, obj))) {
		return Qnil;
	}

	if (rb_block_given_p()) {
		return rb_ensure(rb_yield, obj, fwam_close, obj);
	}

	return obj;
}

static VALUE
fwam_size(obj, o_dir)
VALUE obj, o_dir;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	df_t s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);

	s = wam_size(w, dir);	/* never fails */
	return INT2FIX(s);
}

static VALUE
fwam_elem_num(obj, o_dir, o_id)
VALUE obj, o_dir, o_id;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	idx_t id;
	df_t s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	id = NUM2INT(o_id);

	if ((s = wam_elem_num(w, dir, id)) == -1) {
		return Qnil;
	}
	return INT2FIX(s);
}

static VALUE
fwam_freq_sum(obj, o_dir, o_id)
VALUE obj, o_dir, o_id;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	idx_t id;
	freq_t s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	id = NUM2INT(o_id);

	if ((s = wam_freq_sum(w, dir, id)) == -1) {
		return Qnil;
	}
	return INT2NUM(s);
}

static VALUE
fwam_max_elem_num(obj, o_dir)
VALUE obj, o_dir;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	df_t s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);

	s = wam_max_elem_num(w, dir);	/* never fails */
	return INT2FIX(s);
}

static VALUE
fwam_max_freq_sum(obj, o_dir)
VALUE obj, o_dir;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	freq_t s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);

	s = wam_max_freq_sum(w, dir);	/* never fails */
	return INT2NUM(s);
}

static VALUE
fwam_total_elem_num(obj)
VALUE obj;
{
	struct wamdata *wamp;
	WAM *w;
	freq_t s;

	GetWAM2(obj, wamp, w);
	s = wam_total_elem_num(w);	/* never fails */
	return INT2NUM(s);
}

static VALUE
fwam_total_freq_sum(obj)
VALUE obj;
{
	struct wamdata *wamp;
	WAM *w;
	freq_t s;

	GetWAM2(obj, wamp, w);
	s = wam_total_freq_sum(w);	/* never fails */
	return INT2NUM(s);
}

static VALUE
fwam_id2name(obj, o_dir, o_id)
VALUE obj, o_dir, o_id;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	idx_t id;
	char const *name;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	id = NUM2INT(o_id);

	if (!(name = wam_id2name(w, dir, id))) {
		return Qnil;
	}
	return rb_str_new2(name);
}

static VALUE
fwam_name2id(obj, o_dir, o_name)
VALUE obj, o_dir, o_name;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	char const *name;
	idx_t id;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	SafeStringValue(o_name);
	name = RSTRING_PTR(o_name);

	if (!(id = wam_name2id(w, dir, name))) {
		return Qnil;
	}
	return INT2FIX(id);
}

static VALUE
fwam_get_vec(obj, o_dir, o_id)
VALUE obj, o_dir, o_id;
{
	struct wamdata *wamp;
	WAM *w;
	int dir;
	idx_t id;
	struct xr_vec const *v;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	id = NUM2INT(o_id);

	if (wam_get_vec(w, dir, id, &v) == -1) {
		return Qnil;
	}
	return vec2value(v, w, dir);
}

static VALUE
fwam_prop_nentries(obj, o_dir, o_file)
VALUE obj, o_dir, o_file;
{
	struct wamdata *wamp;
	WAM *w, *p;
	int dir;
	char const *file;
	df_t s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	SafeStringValue(o_file);
	file = RSTRING_PTR(o_file);

	if (!(p = wam_prop_open(w, dir, file))) {
		return Qnil;
	}
	s = wam_prop_nentries(p);	/* never fails */
	/* do not close */
	return INT2FIX(s);
}

static VALUE
fwam_prop_gets(obj, o_dir, o_file, o_id)
VALUE obj, o_dir, o_file, o_id;
{
	struct wamdata *wamp;
	WAM *w, *p;
	int dir;
	char const *file;
	idx_t id;
	char const *s;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	SafeStringValue(o_file);
	file = RSTRING_PTR(o_file);
	id = NUM2INT(o_id);

	if (!(p = wam_prop_open(w, dir, file))) {
		return Qnil;
	}
	if (wam_prop_gets(p, id, &s) <= 0) {
		return Qnil;
	}
	/* do not close */
	return rb_str_new2(s);
}

static VALUE
fwam_fss_gets(obj, o_id, o_segid)
VALUE obj, o_id, o_segid;
{
	struct wamdata *wamp;
	WAM *w, *f;
	idx_t id;
	unsigned int segid;
	char const *s;

	GetWAM2(obj, wamp, w);
	id = NUM2INT(o_id);
	segid = NUM2INT(o_segid);

	if (!(f = wam_fss_open(w))) {
		return Qnil;
	}
	if (wam_fss_gets(f, id, segid, &s) <= 0) {
		return Qnil;
	}
	/* do not close */
	return rb_str_new2(s);
}

static VALUE
wsh0(obj, s, e, f, o_dir, o_weight_type, o_maxnd, o_tak_iw_limit, obj2)
VALUE obj, s, e, f, o_dir, o_weight_type, o_maxnd, o_tak_iw_limit, obj2;
{
	struct wamdata *wamp, *wamp2;
	WAM *w, *wt;
	struct syminfo *q = NULL;
	struct bxue_t *qe = NULL;
	struct fss_query *fq = NULL;
	df_t len, nd, total;
	int dir, weight_type;
	struct syminfo_conv_cookie c;
	struct syminfo *dp = NULL;
	VALUE r;
	size_t lp;

	GetWAM2(obj, wamp, w);
	dir = NUM2INT(o_dir);
	weight_type = NUM2INT(o_weight_type);
	nd = NUM2INT(o_maxnd);

	if (obj2 != Qnil) GetWAM2(obj2, wamp2, wt); else wt = NULL;

	c.w = w;
	c.dir = dir;
	if (!(q = value2array(s, &lp, sizeof *q, value2syminfo_e_r, NULL, &c))) {
		return Qnil;
	}
	len = lp;

	if (e != Qnil) {
		df_t elen;
		if (!(qe = value2bxue_t(e, w, dir, &lp))) {
			free_vecs(q, len);
			free(q);
			return Qnil;
		}
		elen = lp;
		dp = bex_wsh(q, len, w, dir, weight_type, &nd, &total, qe, elen, NULL);
		free(qe);
	}
	else if (f != Qnil) {
		if (!(fq = value2fss_query(f))) {
			free_vecs(q, len);
			free(q);
			return Qnil;
		}
		dp = fss_wsh(q, len, w, dir, weight_type, &nd, &total, NULL, fq);
		free_fq(fq);
	}
	else {
		dp = wsh(q, len, w, dir, weight_type, &nd, &total, NULL);
	}
	free_vecs(q, len);
	free(q);

	if (!dp) {
		return Qnil;
	}

	qsort(dp, nd, sizeof *dp, sym_compar);

	if (o_tak_iw_limit != Qnil && wt) {
		nd = iwl(dp, nd, NUM2INT(o_tak_iw_limit), wt, WAM_REVERT_DIRECTION(dir));
	}

	last_nd = total;

	c.w = w;
	c.dir = WAM_REVERT_DIRECTION(dir);
	if ((r = array2value(dp, nd, sizeof *dp, syminfo2value, &c)) == Qnil) {
		return Qnil;
	}

	free(dp);
	return r;
}

static VALUE
ncsb0(obj, s, o_cs_type, o_dir, o_weight_type, o_elemn, o_cno, o_cswmax)
VALUE obj, s, o_cs_type, o_dir, o_weight_type, o_elemn, o_cno, o_cswmax;
{
	struct wamdata *wamp;
	WAM *w;

	struct syminfo *q = NULL;
	int cs_type, dir, weight_type, elemn, cno, cswmax;
	struct syminfo_conv_cookie c;
	df_t len;
	size_t lp;

	struct cs_elem *cc;
	VALUE r;

	GetWAM2(obj, wamp, w);
	cs_type = NUM2INT(o_cs_type);
	dir = NUM2INT(o_dir);
	weight_type = NUM2INT(o_weight_type);
	elemn = NUM2INT(o_elemn);
	cno = NUM2INT(o_cno);
	cswmax = NUM2INT(o_cswmax);

	c.w = w;
	c.dir = dir;
	if (!(q = value2array(s, &lp, sizeof *q, value2syminfo_e_r, NULL, &c))) {
		return Qnil;
	}
	len = lp;

	cc = ncsb(q, len, cs_type, w, dir, weight_type, elemn, &cno, weight_type, cswmax);
	free_vecs(q, len);
	free(q);

	if (!cc) {
		return Qnil;
	}

	r = array2value(cc, cno, sizeof *cc, cc2value, &c);

	free(cc);
	return r;
}

static VALUE
fwam_wsh(obj, s, dir, weight_type, n, tak_iw_limit, obj2)
VALUE obj, s, dir, weight_type, n, tak_iw_limit, obj2;
{
	return wsh0(obj, s, Qnil, Qnil, dir, weight_type, n, tak_iw_limit, obj2);
}

static VALUE
fwam_bex_wsh(obj, s, e, dir, weight_type, n, tak_iw_limit, obj2)
VALUE obj, s, e, dir, weight_type, n, tak_iw_limit, obj2;
{
	return wsh0(obj, s, e, Qnil, dir, weight_type, n, tak_iw_limit, obj2);
}

static VALUE
fwam_fss_wsh(obj, s, f, dir, weight_type, n, tak_iw_limit, obj2)
VALUE obj, s, f, dir, weight_type, n, tak_iw_limit, obj2;
{
	return wsh0(obj, s, Qnil, f, dir, weight_type, n, tak_iw_limit, obj2);
}

static VALUE
fwam_ncsb(s, cs_type, d, mode, weight_type, elemn, cno, cswmax)
VALUE s, cs_type, d, mode, weight_type, elemn, cno, cswmax;
{
	return ncsb0(s, cs_type, d, mode, weight_type, elemn, cno, cswmax);
}

static VALUE
fwam_text2vec(obj, o_text, o_stemmer)
VALUE obj, o_text, o_stemmer;
{
	struct wamdata *wamp;
	WAM *w;
	char const *text, *stemmer;
	struct xr_vec *v;
	VALUE r;

	GetWAM2(obj, wamp, w);

	SafeStringValue(o_text);
	text = RSTRING_PTR(o_text);
	if (o_stemmer != Qnil) {
		SafeStringValue(o_stemmer);
		stemmer = RSTRING_PTR(o_stemmer);
	}
	else {
		stemmer = "auto";
	}
	if (!(v = text2vec(text, w, WAM_COL, stemmer))) {
		return Qnil;
	}

	r = vec2value(v, w, WAM_COL);
	free(v);
	return r;
}

static VALUE
fwam_get_last_nd(obj)
VALUE obj;
{
/*
	struct wamdata *wamp;
	WAM *w;
*/

	return INT2NUM(last_nd);
}

/**************************************************************************
 * UTILITIES (2)
 **************************************************************************/

static VALUE
vec2value(v, w, d)
struct xr_vec const *v;
WAM *w;
{
	VALUE r, rv;
	df_t i;

	r = rb_obj_alloc(cWAM_XR_VEC);
	rv = array2value((void *)v->elems, v->elem_num, sizeof *v->elems, elem2value, NULL);
	rb_iv_set(r, "@elems", rv);
	rb_iv_set(r, "@elem_num", INT2FIX(v->elem_num));
	rb_iv_set(r, "@freq_sum", INT2FIX(v->freq_sum));
	return r;
}

static int
elem2value(v, cookie, e)
void *v, *cookie;
VALUE *e;
{
	struct xr_elem *x = v;

	*e = rb_obj_alloc(cWAM_XR_ELEM);
	rb_iv_set(*e, "@id", INT2FIX(x->id));
	rb_iv_set(*e, "@freq", INT2FIX(x->freq));
	return 0;
}

static int
cc2value(v, cookie, e)
void *v, *cookie;
VALUE *e;
{
	VALUE rr;
	struct syminfo *dp;
	size_t nd;
	struct syminfo_conv_cookie *c;
	struct cs_elem *cci;
	int status = 0;

	cci = v;
	c = cookie;

	dp = cci->csa.s;
	nd = cci->csa.n;
	qsort(dp, nd, sizeof *dp, sym_compar);
	if ((*e = array2value(dp, nd, sizeof *dp, syminfo2value, c)) == Qnil) {
		status = -1;
	}
	free(dp);
	free(cci->csw.s);

	return status;
}

static int
charpp2value(v, cookie, e)
void *v, *cookie;
VALUE *e;
{
	char **s = v;

	if ((*e = rb_str_new2(*s)) == Qnil) {
		return -1;
	}
	return 0;
}

static VALUE
charpp2intern(v, cookie)
void *v, *cookie;
{
	char **s = v;

	return rb_intern(*s);
}

static int
syminfo2value(v, cookie, e)
void *v, *cookie;
VALUE *e;
{
	struct syminfo *x = v;
	struct syminfo_conv_cookie *c = cookie;
	char const *name;

	*e = rb_obj_alloc(cWAM_SYMINFO_ELEM);
	rb_iv_set(*e, "@v", Qnil);
	rb_iv_set(*e, "@id", INT2FIX(x->id));
	rb_iv_set(*e, "@TF", INT2FIX(x->TF));
	rb_iv_set(*e, "@TF_d", INT2FIX(x->TF_d));
	rb_iv_set(*e, "@DF", INT2FIX(x->DF));
	rb_iv_set(*e, "@DF_d", INT2FIX(x->DF_d));
	rb_iv_set(*e, "@weight", rb_float_new(x->weight));
	if (!(name = wam_id2name(c->w, c->dir, x->id))) {
		return -1;
	}
	rb_iv_set(*e, "@name", rb_str_new2(name));
	return 0;
}

static int
value2syminfo_e_r(v, e, cookie)
VALUE v;
void *e, *cookie;
{
	struct syminfo *x = e;
	struct syminfo_conv_cookie *c = cookie;
	VALUE s;

	x->v = NULL;
	x->id = 0;
	if ((s = rb_iv_get(v, "@v")) != Qnil) {
		if (!(x->v = value2vec(s, c))) {
			return -1;
		}
	}
	else if ((s = rb_iv_get(v, "@name")) != Qnil) {
		char *name;
		SafeStringValue(s);
		name = RSTRING_PTR(s);
		x->id = wam_name2id(c->w, c->dir, name);
	}
	else if ((s = rb_iv_get(v, "@id")) != Qnil) {
		x->id = NUM2INT(s);
	}
	else {
		return -1;
	}
	x->TF = NUM2INT(rb_iv_get(v, "@TF"));
	x->TF_d = NUM2INT(rb_iv_get(v, "@TF_d"));
	x->DF = NUM2INT(rb_iv_get(v, "@DF"));
	x->DF_d = NUM2INT(rb_iv_get(v, "@DF_d"));
	x->weight = NUM2DBL(rb_iv_get(v, "@weight"));

	return 0;
}

static struct bxue_t *
value2bxue_t(v, w, dir, n)
VALUE v;
WAM *w;
ssize_t *n;
{
	struct syminfo_conv_cookie c;
	c.w = w;
	c.dir = dir;
	return value2array(v, n, sizeof (struct bxue_t), value2bxue_t_r, NULL, &c);
}

static int
value2bxue_t_r(v, e, cookie)
VALUE v;
void *e, *cookie;
{
	struct bxue_t *x = e;
	struct syminfo_conv_cookie *c = cookie;
	VALUE s;

	x->id = 0;
	x->type = '"';
	if ((s = rb_iv_get(v, "@name")) != Qnil) {
		char *name;
		SafeStringValue(s);
		name = RSTRING_PTR(s);
		x->id = wam_name2id(c->w, c->dir, name);
	}
	else if ((s = rb_iv_get(v, "@id")) != Qnil) {
		x->id = NUM2INT(s);
	}
	else if (s = rb_iv_get(v, "@type")) {
		x->type = s ? NUM2INT(s) : 0;
	}
	else {
		return -1;
	}
	return 0;
}

static struct fss_query *
value2fss_query(v)
VALUE v;
{
	struct fss_query *r;

	struct fss_con_query *cq;
	size_t len;

	if (!(r = calloc(1, sizeof *r))) {
		perror("calloc");
		return NULL;
	}

	if (!(cq = value2array(v, &len, sizeof *cq, value2fss_con_query_r, NULL, NULL))) {
		return NULL;
	}

	r->query = cq;
	r->nq = len;
	return r;
}

static int
value2fss_con_query_r(v, e, cookie)
VALUE v;
void *e, *cookie;
{
	struct fss_con_query *x = e;

	struct fss_simple_query *sq;
	size_t len;

	if (!(sq = value2array(v, &len, sizeof *sq, value2fss_simple_query_r, NULL, NULL))) {
		return -1;
	}

	x->q = sq;
	x->n = len;
	return 0;
}

static int
value2fss_simple_query_r(v, e, cookie)
VALUE v;
void *e, *cookie;
{
	struct fss_simple_query *x = e;

	VALUE s;

	x->pattern = NULL;
	x->negativep = 0;
	x->segments = ANYSEG();
	x->options = 0;
	if ((s = rb_iv_get(v, "@pattern")) != Qnil) {
		SafeStringValue(s);
		x->pattern = RSTRING_PTR(s);
	}
	if ((s = rb_iv_get(v, "@negativep")) != Qnil) {
		x->negativep = FIX2INT(s);
	}
	if ((s = rb_iv_get(v, "@segments")) != Qnil) {
		SafeStringValue(s);
		if (parse_segs(RSTRING_PTR(s), &x->segments) == -1) {  
			return -1;
		}
	}
	if ((s = rb_iv_get(v, "@options")) != Qnil) {
		SafeStringValue(s);
		if (parse_opts(RSTRING_PTR(s), &x->options) == -1) {  
			return -1;
		}
	}
	return 0;
}

static struct xr_vec *
value2vec(v, c)
VALUE v;
struct syminfo_conv_cookie *c;
{
	size_t i, elem_num, freq_sum;
	size_t nm, fs;
	VALUE elems;
	struct xr_vec *r;

	elem_num = FIX2INT(rb_iv_get(v, "@elem_num"));
	freq_sum = FIX2INT(rb_iv_get(v, "@freq_sum"));
	elems = rb_iv_get(v, "@elems");

	nm = RARRAY_LEN(elems);
	assert(nm = elem_num);

	if (!(r = malloc(sizeof_dxr_vec(elem_num)))) {
		perror("malloc");
		return NULL;
	}
	r->elem_num = elem_num;
	for (i = 0, fs = 0; i < elem_num; i++) {
		VALUE f;
		f = rb_ary_entry(elems, (long)i);
		struct xr_elem *e = &r->elems[i];
		e->id = FIX2INT(rb_iv_get(f, "@id"));
		e->freq = FIX2INT(rb_iv_get(f, "@freq"));
		fs += e->freq;
	}
	assert(fs == freq_sum);
	r->freq_sum = freq_sum;
	return r;
}

static VALUE
array2value(base, nmemb, size, fn, cookie)
void *base, *cookie;
size_t nmemb, size;
int (*fn)(void *, void *, VALUE *);
{
	VALUE v;
	size_t i;

	if ((v = rb_ary_new2(nmemb)) == Qnil) {
		return Qnil;
	}
	for (i = 0; i < nmemb; i++) {
		VALUE e;
		if ((*fn)((void *)((char *)base) + i * size, cookie, &e) == -1) {
			return Qnil;
		}
		rb_ary_push(v, e);
	}
	return v;
}

static void *
value2array(v, nmemb, size, fn, cln, cookie)
VALUE v;
size_t *nmemb, size;
int (*fn)(VALUE, void *, void *);
void (*cln)(void *, void *);
void *cookie;
{
	void *base;
	size_t i;

	*nmemb = RARRAY_LEN(v);
	if (!(base = malloc(*nmemb * size))) {
		perror("malloc");
		return NULL;
	}
	for (i = 0; i < *nmemb; i++) {
		VALUE e;
		e = rb_ary_entry(v, (long)i);
		if ((*fn)(e, (void *)((char *)base) + i * size, cookie) == -1) {
			goto cleanup;
		}
	}

	return base;

cleanup:
	if (cln) {
		size_t j;
		for (j = 0; j < i; j++) {
			(*cln)((void *)((char *)base) + j * size, cookie);
		}
	}
	free(base);

	return NULL;
}

static void
free_vecs(q, len)
struct syminfo *q;
size_t len;
{
	size_t i;
	for (i = 0; i < len; i++) {
		free((void *)q[i].v);
	}
}

static void
free_fq(fq)
struct fss_query *fq;
{
	size_t i, j;
	for (i = 0; i < fq->nq; i++) {
		struct fss_con_query *c = &fq->query[i];
		for (j = 0; j < c->n; j++) {
			/*cleanup(c->q[j]);*/
		}
		free(c);
	}
	free(fq);
}

static size_t
iwl(terms, nkwords, tak_iw_limit, w, dir)
struct syminfo *terms;
size_t nkwords, tak_iw_limit;
WAM *w;
{
	size_t i, j;

	if (nkwords <= 0) {
		return nkwords;
	}
	/* note: terms[j..nkwords] ni gomi ga nokoru. */
	for (i = 0, j = 0; i < nkwords; i++) {
		if (wam_elem_num(w, dir, terms[i].id) <= tak_iw_limit) {
			terms[j++] = terms[i];
		}
	}
	return MAX(j, 1);
}

void
Init_rnwam()
{
	static char const *weight_types[] = WEIGHT_TYPE_NAMES;
	static char const *cs_types[] = CS_TYPE_NAMES;
	int i;
	VALUE r;

	cWAM = rb_define_class("WAM", rb_cObject);
	eWAMError = rb_define_class("WAMError", rb_eStandardError);

	rb_define_alloc_func(cWAM, fwam_alloc);
	rb_define_singleton_method(cWAM, "open", fwam_s_open, -1);

	cWAM_XR_VEC = rb_define_class("WAM_XR_VEC", rb_cObject);
	rb_define_attr(cWAM_XR_VEC, "elem_num", 1, 1);
	rb_define_attr(cWAM_XR_VEC, "freq_sum", 1, 1);
	rb_define_attr(cWAM_XR_VEC, "elems", 1, 1);

	cWAM_XR_ELEM = rb_define_class("WAM_XR_ELEM", rb_cObject);
	rb_define_attr(cWAM_XR_ELEM, "id", 1, 1);
	rb_define_attr(cWAM_XR_ELEM, "freq", 1, 1);

	cWAM_SYMINFO_ELEM = rb_define_class("WAM_SYMINFO_ELEM", rb_cObject);
	rb_define_attr(cWAM_SYMINFO_ELEM, "v", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "id", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "TF", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "TF_d", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "DF", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "DF_d", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "weight", 1, 1);
	rb_define_attr(cWAM_SYMINFO_ELEM, "name", 1, 1);

	cWAM_BXU_ELEM = rb_define_class("WAM_BXU_ELEM", rb_cObject);
	rb_define_attr(cWAM_BXU_ELEM, "name", 1, 1);
	rb_define_attr(cWAM_BXU_ELEM, "id", 1, 1);
	rb_define_attr(cWAM_BXU_ELEM, "type", 1, 1);

	cWAM_SIMPLE_QUERY = rb_define_class("WAM_SIMPLE_QUERY", rb_cObject);
	rb_define_attr(cWAM_SIMPLE_QUERY, "pattern", 1, 1);
	rb_define_attr(cWAM_SIMPLE_QUERY, "negativep", 1, 1);
	rb_define_attr(cWAM_SIMPLE_QUERY, "segments", 1, 1);
	rb_define_attr(cWAM_SIMPLE_QUERY, "options", 1, 1);

	rb_define_method(cWAM, "initialize", fwam_initialize, -1);
	rb_define_method(cWAM, "close", fwam_close, 0);
	rb_define_method(cWAM, "closed?", fwam_closed, 0);
	rb_define_method(cWAM, "size", fwam_size, 1);
	rb_define_method(cWAM, "elem_num", fwam_elem_num, 2);
	rb_define_method(cWAM, "freq_sum", fwam_freq_sum, 2);
	rb_define_method(cWAM, "max_elem_num", fwam_max_elem_num, 1);
	rb_define_method(cWAM, "max_freq_sum", fwam_max_freq_sum, 1);
	rb_define_method(cWAM, "total_elem_num", fwam_total_elem_num, 0);
	rb_define_method(cWAM, "total_freq_sum", fwam_total_freq_sum, 0);

	rb_define_method(cWAM, "id2name", fwam_id2name, 2);
	rb_define_method(cWAM, "name2id", fwam_name2id, 2);

	rb_define_method(cWAM, "prop_nentries", fwam_prop_nentries, 2);
	rb_define_method(cWAM, "prop_gets", fwam_prop_gets, 3);
	rb_define_method(cWAM, "fss_gets", fwam_fss_gets, 2);

	rb_define_method(cWAM, "get_vec", fwam_get_vec, 2);
	rb_define_method(cWAM, "wsh", fwam_wsh, 6);
	rb_define_method(cWAM, "bex_wsh", fwam_bex_wsh, 7);
	rb_define_method(cWAM, "fss_wsh", fwam_fss_wsh, 7);
	rb_define_method(cWAM, "ncsb", fwam_ncsb, 7);

/*
-- ??? xfss ???
*/

	rb_define_method(cWAM, "text2vec", fwam_text2vec, 2);
	rb_define_method(cWAM, "get_last_nd", fwam_get_last_nd, 0);

/*
-- ??? set_prob_mode ???
-- ??? wam_setopt ???
*/
	rb_define_const(cWAM, "ROW", INT2FIX(WAM_ROW));
	rb_define_const(cWAM, "COL", INT2FIX(WAM_COL));
	for (i = 0; i < nelems(weight_types); i++) {
		rb_define_const(cWAM, weight_types[i], INT2FIX(i));
	}
	r = array2value(weight_types, nelems(weight_types), sizeof *weight_types, charpp2value, NULL);
	assert(r != Qnil);
	rb_define_const(cWAM, "WEIGHT_TYPES", r);
	for (i = 0; i < nelems(cs_types); i++) {
		if (!*cs_types[i]) continue;
		rb_define_const(cWAM, cs_types[i], INT2FIX(i));
	}
	r = array2value(cs_types, nelems(cs_types), sizeof *cs_types, charpp2value, NULL);
	assert(r != Qnil);
	rb_define_const(cWAM, "CS_TYPES", r);
}
