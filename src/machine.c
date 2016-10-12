#include "roadef.h"

#ifdef __cplusplus
extern "C" {
#endif

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
	machine_t *m = calloc(1,sizeof(machine_t));

	m->nres = nres;
	m->nmach = nmach;

	m->cap = calloc(nres, sizeof(int64_t));
	m->safecap = calloc(nres, sizeof(int64_t));

	m->mov_cost = calloc(nmach, sizeof(int64_t));
	
	return m;
    }
    
    void
    machine_destroy(machine_t **self_p){
	assert(self_p);
	if(!*self_p) return;

	machine_t *self = *self_p;
	
	if (self->cap) free(self->cap);
	if (self->safe_cap) free(self->safe_cap);
	if (self->mov_cost) free(self->mov_cost);
	free(self);
	*self_p=NULL;
    }

    int64_t
    machine_location(machine_t *self){
    }

    int64_t
    machine_neigh(machine_t *self){
    }
    
    int64_t
    machine_cap(machine_t *self, int64_t res){
    }
    
    int64_t
    machine_safecap(machine_t *self, int64_t res){

    }
    
    int64_t
    machine_mvcost(machine_t *self, int64_t mach){
    }
    

    
#ifdef __cplusplus
};
#endif
