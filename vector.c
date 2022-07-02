#define _GNU_SOURCE
#include <stdlib.h>
#include "vector.h"
#include "tools.h"

#define VECTOR_GROWTH_FACTOR 2

static void* vector_resize(Vector* v, size_t num);
static int vector_bsearch_(Vector* v, void* key, Vector_comp comp,
                           void* comp_data, int left, int right);

Vector* vector_init(Vector_free_val free_val)
{
	Vector* v = xcalloc(1, sizeof(*v));

	v->free_val = free_val;

	return v;
}

void vector_free(Vector* v)
{
	unsigned int i;

	if (v) {
		if (v->free_val)
			for (i = 0; i < v->size; ++i)
				v->free_val(v->val[i]);

		free(v->val);
		free(v);
	}
}

void vector_free_generic(void* v)
{
	vector_free(v);
}

static void* vector_resize(Vector* v, size_t num)
{
	void* tmp;

	tmp = xrealloc(v->val, num * sizeof(*v->val));

	v->val = tmp;
	v->capacity = num;

	return v->val;
}

size_t vector_size(const Vector* v)
{
	return v->size;
}

void vector_push(Vector* v, void* val)
{
	if (v->size == v->capacity)
		vector_resize(v, (v->capacity +1)*VECTOR_GROWTH_FACTOR);

	v->val[v->size] = val;
	v->size++;
}

void vector_insert(Vector* v, void* val, int pos)
{
	if (v->size == v->capacity)
		vector_resize(v, (v->capacity +1)*VECTOR_GROWTH_FACTOR);

	if (pos < v->size)
		for (int i = v->size -1; i >= pos; i--)
			v->val[i +1] = v->val[i];

	v->val[pos] = val;
	v->size++;
}

void* vector_pop(Vector* v)
{
	if (v->size == 0)
		return NULL;

	v->size--;
	return v->val[v->size];
}

void vector_del(Vector* v, int pos)
{
	if (v->free_val)
		v->free_val(v->val[pos]);

	for (int i = pos; i < v->size -1; ++i)
		v->val[i] = v->val[i +1];

	v->size--;
}

void* vector_extract(Vector* v, int pos)
{
	void* data = v->val[pos];

	for (int i = pos; i < v->size -1; ++i)
		v->val[i] = v->val[i +1];

	v->size--;

	return data;
}

/* Keep the first d elements */
void vector_trunc(Vector* v, int d)
{
	if (v->size > d) {
		for (int i = d; i < v->size; ++i)
			if (v->free_val)
				v->free_val(v->val[i -1]);

		v->size = d;
	}
}

void* vector_get(const Vector* v, int pos)
{
	if (pos >= v->size)
		return NULL;

	return v->val[pos];
}

void vector_concat(Vector* dst, const Vector* src)
{
	if (src == NULL)
		return;

	if (dst->capacity -dst->size < src->size)
		vector_resize(dst, dst->capacity +src->size);

	for (int i = 0; i < src->size; ++i) {
		dst->val[dst->size] = src->val[i];
		dst->size++;
	}
}

Vector* vector_dup(const Vector* src, Vector_copy_val copy_val)
{
	Vector* dst;
	int i;

	if (copy_val) {
		dst = vector_init(src->free_val);
		vector_resize(dst, src->size);

		for (i = 0; i < src->size; i++)
			dst->val[i] = copy_val(src->val[i]);
	}
	else {
		dst = vector_init(NULL);
		vector_resize(dst, src->size);

		for (i = 0; i < src->size; i++)
			dst->val[i] = src->val[i];
	}

	dst->size = src->size;

	return dst;
}

void vector_sort(Vector* v, Vector_comp comp, void* comp_data)
{
	qsort_r(v->val, v->size, sizeof(void*), comp, comp_data);
}

/* Removes duplicate entries. Run vector_sort first! */
void vector_rmdup(Vector* v, Vector_comp comp, void* comp_data)
{
	for (int i = 1; i < v->size; ++i) {
		if (!comp(v->val[i], v->val[i -1], comp_data)) {
			vector_del(v, i);
			i--;
		}
	}
}

/* Search for the specified "key" in vector "v" and return its index. If the
 * key is not found the closest match's index is returned */
int vector_bsearch(Vector* v, void* key, Vector_comp comp, void* comp_data)
{
	return vector_bsearch_(v, key, comp, comp_data, 0, v->size -1);
}

static int vector_bsearch_(Vector* v, void* key, Vector_comp comp,
                           void* comp_data, int left, int right)
{
	int mid, ret;

	if (right -left < 2)
		return right;

	mid = (left +right)/2;

	ret = comp(key, v->val[mid], comp_data);
	if (ret > 0)
		return vector_bsearch_(v, key, comp, comp_data, mid, right);

	if (ret < 0)
		return vector_bsearch_(v, key, comp, comp_data, left, mid);

	return mid;
}

