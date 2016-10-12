#include "roadef.h"

struct machine_t {	
    int64_t location;
    int64_t neighbor;	
    size_t nres;
    int64_t *cap;
    int64_t *safecap;
    size_t nmach;
    int64_t *mov_cost;
};

machine_t *
machine_new(size_t nres, size_t nmach, char *line){
    assert(nres);
    assert(nmach);
    assert(strcmp("",line));
	
    machine_t *m = calloc(1,sizeof(machine_t));

    m->nres = nres;
    m->nmach = nmach;

    m->cap = calloc(nres, sizeof(int64_t));
    if (!m->cap) {
	machine_destroy(&m);
	return NULL;
    }
    m->safecap = calloc(nres, sizeof(int64_t));
    if (!m->safecap) {
	machine_destroy(&m);
	return NULL;
    }

    m->mov_cost = calloc(nmach, sizeof(int64_t));
    if (!m->mov_cost) {
	machine_destroy(&m);
	return NULL;
    }
	
    return m;
}
    
void
machine_destroy(machine_t **self_p){
    assert(self_p);
	
    if(!*self_p) return;

    machine_t *self = *self_p;
	
    if (self->cap) free(self->cap);
    if (self->safecap) free(self->safecap);
    if (self->mov_cost) free(self->mov_cost);
    free(self);
	
    *self_p=NULL;
}

int64_t
machine_location(machine_t *self){
    assert(self);

    return self->location;
}

int64_t
machine_neigh(machine_t *self){
    assert(self);

    return self->neighbor;
}
    
int64_t
machine_cap(machine_t *self, int64_t res){
    assert(self);
    assert(res < self->nres);

    return self->cap[res];
}
    
int64_t
machine_safecap(machine_t *self, int64_t res){
    assert(self);
    assert(res<self->nres);

    return self->safecap[res];
}
    
int64_t
machine_mvcost(machine_t *self, int64_t mach){
    assert(self);
    assert(mach < self->nmach);

    return self->mov_cost[mach];
}
    
