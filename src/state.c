#include "roadef.h"

struct state_t {
    model_t *model;
    int64_t nproc;
    int64_t *procmap_begin;
    int64_t *procmap_curr;
    bool needs_update;
    int64_t nres;
    int64_t *obj_resources;
    int64_t nbalance;
    int64_t *obj_balance;
    int64_t obj_pmc;
    int64_t obj_smc;
    int64_t obj_mmc;

    int64_t **utilization;
    int64_t **tr_utilization;
};

state_t *
state_new(model_t *model) {
    assert(model);

    state_t * state = calloc (1 , sizeof(state_t));
    if (!state) return NULL;

    state->model = model;
    state->nproc = model_nproc(model);
    state->nres = model_nres(model);
    state->nbalance = model_nbalance(model);

    state->procmap_begin = calloc(state->nproc,sizeof(int64_t));
    state->procmap_curr = calloc(state->nproc,sizeof(int64_t));

    state->obj_resources=calloc(state->nres,sizeof(int64_t));
    state->obj_balance=calloc(state->nres,sizeof(int64_t));

    if (!state->procmap_begin ||
      !state->procmap_curr ||
      !state->obj_resources ||
      !state->obj_balance ||
      false)
	state_destroy(&state);
	
    return state;
}

void
state_destroy(state_t **self_p) {
    assert(self_p);

    if (!*self_p) return;

    state_t *self = *self_p;

    if (self->procmap_begin)
	free(self->procmap_begin);

    if (self->procmap_curr)
	free(self->procmap_curr);

    if (self->obj_resources)
	free(self->obj_resources);

    if (self->obj_balance)
	free(self->obj_balance);


    free(self);
    *self_p = NULL;

}

void
state_load(state_t *state, char *line){
    assert(state);
    assert(line);
    assert(strneq("",line));

    
}

static void
s_state_swap_procmap(state_t *state){
    int64_t *t = state->procmap_begin;
    state->procmap_begin = state->procmap_curr;
    state->procmap_curr = t;
}


void
state_step(state_t *state) {
    assert(state);
    s_state_swap_procmap(state);
    state->needs_update=true;
}

bool
state_validate(state_t *state) {
    assert(state);

    state->needs_update = false;
    return true;
}

int64_t
state_obj(state_t *state) {
    assert(state);
    if (state->needs_update && !state_validate(state))
	    return -1;
    
    int64_t obj = 0;
    for (int64_t i =0 ; i< state->nres;i++ ) 	
	obj+= resource_loadcost(model_resource(state->model,i)) *
	  (state->obj_resources[i]);

    for (int64_t i =0 ; i< state->nbalance;i++ ) 	
	obj+= balance_weightcost(model_balance(state->model,i)) *
	  (state->obj_balance[i]);

    obj += model_wpmc(state->model) * (state->obj_pmc);
    obj += model_wsmc(state->model) * (state->obj_smc);
    obj += model_wmmc(state->model) * (state->obj_mmc);
    
    return obj;
}

void
state_test(bool verbose) {

    if(verbose) printf("  * state: ");

    char *content = strdup(MODEL_EXAMPLE1);
    model_t *model = model_new(content);
    free(content);

    state_t *state = state_new(model);
    assert(state);
    state_destroy(&state);
    assert(!state);
    state_destroy(&state);

    state = state_new(model);

    
    state_destroy(&state);
    model_destroy(&model);
      

    if(verbose)
	printf("OK\n");
}

	  
