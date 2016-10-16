#include "roadef.h"

struct machine_t {	
    int64_t location;
    int64_t neighbor;	
    int64_t nres;
    int64_t *cap;
    int64_t *safecap;
    int64_t nmach;
    int64_t *mov_cost;
};

machine_t *
machine_new(int64_t nres, int64_t nmach, char *line){
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

    char *endtok, *tok ;
    
    tok = strtok_r(line," ",&endtok);
    m->neighbor = strtoul(tok,NULL,10);
    tok = strtok_r(NULL," ",&endtok);
    m->location = strtoul(tok,NULL,10);

    for (int64_t i=0; i< nres; i++){
	tok = strtok_r(NULL," ",&endtok);
	m->cap[i] = strtoul(tok,NULL,10);
    }

    for (int64_t i=0; i< nres; i++){
	tok = strtok_r(NULL," ",&endtok);
	m->safecap[i] = strtoul(tok,NULL,10);
    }

    for (int64_t i=0; i< nmach; i++){
	tok = strtok_r(NULL," ",&endtok);
	m->mov_cost[i] = strtoul(tok,NULL,10);
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
    

void
machine_test(bool verbose){

    if (verbose) printf("  * machine: ");
    char *line;
    machine_t *m;

    line = strdup("1 2 10 1000 8 80 5 4 3 2 0");
    m = machine_new(2,4,line);
    free(line);
    assert(m);
    machine_destroy(&m);
    assert(!m);
    machine_destroy(&m);

    line = strdup("1 2 10 1000 8 80 5 4 3 2 0");
    m = machine_new(2,4,line);
    free(line);

    assert(1==machine_neigh(m));
    assert(2==machine_location(m));
    assert(10==machine_cap(m, 0));
    assert(1000==machine_cap(m, 1));
    assert(8==machine_safecap(m, 0));
    assert(80==machine_safecap(m, 1));
    assert(5==machine_mvcost(m, 0));
    assert(4==machine_mvcost(m, 1));
    assert(3==machine_mvcost(m, 2));
    assert(2==machine_mvcost(m, 3));

    machine_destroy(&m);
    
    if (verbose) printf("OK\n");
}
