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
    int64_t *rutilization;	
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

    state->rutilization = calloc(state->nres,sizeof(int64_t));
    state->utilization = calloc(state->nres,sizeof(int64_t *));
    state->tr_utilization = calloc(state->nres,sizeof(int64_t *));

    for (int64_t i =0 ; i< state->nres; i++) {
	state->utilization[i] = calloc(model_nmach(state->model),sizeof(int64_t));
	state->tr_utilization[i] = calloc(model_nmach(state->model),sizeof(int64_t));
	if (!state->utilization[i] ||
	  !state->tr_utilization[i])
	    	state_destroy(&state);
	
    }	
    
    if (state && (
      !state->procmap_begin ||
      !state->procmap_curr ||
      !state->obj_resources ||
      !state->obj_balance ||
      !state->utilization ||
      !state->rutilization ||	
      !state->tr_utilization ||
      false))
	state_destroy(&state);
	
    return state;
}

void
state_destroy(state_t **self_p) {
    assert(self_p);

    if (!*self_p) return;

    state_t *self = *self_p;

    if (self->procmap_begin) free(self->procmap_begin);

    if (self->procmap_curr) free(self->procmap_curr);

    if (self->obj_resources) free(self->obj_resources);

    if (self->obj_balance) free(self->obj_balance);

    if (self->rutilization) free(self->rutilization);
    
    if (self->utilization) {
	for(int64_t i=0; i< self-> nres ; i++ ) {
	    if(self->utilization[i])
		free(self->utilization[i]);
	}
	free(self->utilization);	
    }

    if (self->tr_utilization) {
	for(int64_t i=0; i< self-> nres ; i++ ) {
	    if(self->tr_utilization[i])
		free(self->tr_utilization[i]);
	}
	free(self->tr_utilization);	
    }
    
    free(self);
    *self_p = NULL;

}

int
state_load(state_t *state, char *filename){
    assert(state);
    assert(filename);
    assert(strneq("",filename));

    FILE *fp = fopen(filename,"rb");
    if (!fp) return -1;

    fseek(fp,0L,SEEK_END);
    long fsz = ftell(fp);
    rewind(fp);

    char *filecontent=calloc(fsz+1, sizeof(char));

    if(1!=fread(filecontent, fsz, 1, fp)) {
	fclose(fp);
	free(filecontent);
	return -1;
    }

    fclose(fp);
    int rc= state_string(state, filecontent);
    free(filecontent);

    return rc;
    
}

int
state_string(state_t *state, char *filecontent) {
    assert(state);
    assert(filecontent);
    assert(strneq("",filecontent));

    char *endtok;
    int64_t i=0 ;
    for ( char *tok=strtok_r(filecontent, " \n",&endtok);
	 tok; tok = strtok_r(NULL, " \n",&endtok), i++) {
	state->procmap_begin[i] = strtoul(tok,NULL,10);
	state->procmap_curr[i] = strtoul(tok,NULL,10);
    }
    return 0;
}

int
state_move(state_t *state, int64_t pidx, int64_t midx) {
    assert(state);
    assert(pidx < state->nproc);
    assert(midx <model_nmach(state->model));

    for (int64_t r=0; r < model_nres(state->model); r++) {
	if (state->utilization[r][midx] +
	    process_requirement(model_process(state->model, pidx),r) >
	  machine_cap(model_machine(state->model,midx),r))
	    return -1;
    }
	  
    for (int64_t r=0; r < model_nres(state->model); r++) {
	state->utilization[r][state->procmap_curr[pidx]] -=
	  resource_transient(model_resource(state->model,r))? 0 :
	  process_requirement(model_process(state->model, pidx),r);

	state->utilization[r][midx] +=
	  process_requirement(model_process(state->model, pidx),r);
    }
    state->procmap_curr[pidx] = midx;
    state->needs_update=true;

    return 0;
}

int
state_swap(state_t *state, int64_t p1idx, int64_t p2idx) {
    assert(state);
    assert(p1idx < state->nproc);
    assert(p2idx < state->nproc);

    for (int64_t r=0; r < model_nres(state->model); r++) {
	int64_t midx = state->procmap_curr[p1idx];
	if(
	    state->utilization[r][midx] 
	    - resource_transient(model_resource(state->model,r))? 0 :
	      process_requirement(model_process(state->model, p1idx),r)
	    + process_requirement(model_process(state->model, p2idx),r) 
            > machine_cap(model_machine(state->model,midx),r)
	  )
	    return -1;

	midx = state->procmap_curr[p2idx];
	if (
	    state->utilization[r][midx] -
	    - resource_transient(model_resource(state->model,r))? 0 :
	      process_requirement(model_process(state->model, p2idx),r)
	    + process_requirement(model_process(state->model, p1idx),r) 
	    > machine_cap(model_machine(state->model,midx),r)
	    )
	  
	    return -1;

    }
    
    int64_t m = state->procmap_curr[p1idx];
    state->procmap_curr[p1idx] = state->procmap_curr[p2idx];
    state->procmap_curr[p2idx] = m;
    state->needs_update=true;

    return 0;
}


void
state_step(state_t *state) {
    assert(state);

    for (int64_t i=0; i< state->nproc; i++)
	if (state->procmap_begin[i] != state->procmap_curr[i])
	    state->procmap_begin[i] = state->procmap_curr[i];

    for (int64_t ridx=0; ridx < model_nres(state->model);ridx++)
	for (int64_t midx=0; midx < model_nmach(state->model);midx++)
	    state->utilization[ridx][midx] = 0;


    for (int64_t i=0; i< state->nproc; i++)
	for (int64_t ridx=0; ridx < model_nres(state->model);ridx++){
	    int64_t midx = state->procmap_begin[i];
	    state->utilization[ridx][midx] +=
	      process_requirement(model_process(state->model,i),ridx);
	}	
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

char *
state_current(state_t *state){
    assert(state);

    char *t =NULL;
    int64_t len = asprintf(&t,"%ld ",model_nmach(state->model));
    free(t);

    t = calloc(len * state->nproc + 1, sizeof(char));

    len=0;
    for( int p =0; p < state->nproc; p++)
	len += sprintf(t+len, "%ld ", state->procmap_curr[p]);

    t[len-1]='\0';
    return t;	
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

    content = strdup("0 3 0");
    state_string(state, content);
    free(content);


    assert(state_validate(state));

    content = state_current(state);
    assert(content);
    assert(streq("0 3 0",content));
    free(content);
    
    assert(0==state_move(state, 1, 2));
    
    assert(state_validate(state));

    assert(0==state_move(state, 2, 1));

    assert(state_validate(state));

    content = state_current(state);
    assert(content);

    assert(streq("0 2 1",content));

    free(content);

    state_destroy(&state);
    model_destroy(&model);
      

    if(verbose)
	printf("OK\n");
}

	  
