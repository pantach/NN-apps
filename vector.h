#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>
#include <stddef.h>

#define vector_prepend(v, val) vector_insert(v, val, 0)

typedef void  (*Vector_free_val)(void*);
typedef void* (*Vector_copy_val)(void*);
typedef int   (*Vector_comp)(const void*, const void*, void*);

typedef struct {
	void** val;
	size_t size;
	size_t capacity;

	Vector_free_val free_val;
} Vector;

Vector* vector_init    (Vector_free_val free_val);
void    vector_free    (Vector* v);
size_t  vector_size    (const Vector* v);
void*   vector_get     (const Vector* v, int pos);
void    vector_push    (Vector* v, void* val);
void    vector_insert  (Vector* v, void* val, int pos);
void*   vector_pop     (Vector* v);
void    vector_del     (Vector* v, int pos);
void*   vector_extract (Vector* v, int pos);
void    vector_trunc   (Vector* v, int n);
void    vector_concat  (Vector* dst, const Vector* src);
Vector* vector_dup     (const Vector* src, Vector_copy_val copy_val);
void    vector_sort    (Vector* v, Vector_comp comp, void* comp_data);
void    vector_rmdup   (Vector* v, Vector_comp comp, void* comp_data);
int     vector_bsearch (Vector* v, void* key, Vector_comp comp, void* comp_data);

void    vector_free_generic(void* v);
#endif
