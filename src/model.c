#include "roadef.h"
    
struct model_t {
    int64_t nres;
    resource_t **resources;
    int64_t *wlc;
    int64_t nmach;
    machine_t **machines;
    int64_t nserv;
    service_t **services;
    int64_t nproc;
    process_t **processes;
    int64_t nbalance;
    balance_t **balance;
    int64_t wpmc;
    int64_t wsmc;
    int64_t wmmc;
    int64_t nlocations;
    int64_t **locations;
} ;

model_t *
model_load(char *filename){
    assert(filename);
    assert(strcmp("",filename));

    FILE *fp = fopen(filename,"rb");
    if (!fp ) return NULL;
	

    fseek(fp,0L,SEEK_END);
    long fsz = ftell(fp);
    rewind(fp);

    char *filecontent=calloc(fsz+1, sizeof(char));

    if(1!=fread(filecontent, fsz, 1, fp)) {
	fclose(fp);
	free(filecontent);
	return NULL;
    }

    fclose(fp);
    model_t * model=  model_new(filecontent);
    free(filecontent);
    return model;
}

model_t *
model_new(char *str){
    assert(str);
    assert(strneq("",str));

    model_t *model = calloc(1, sizeof(model_t));
    char *endstr;

    // Number of resources
    char *line = strtok_r(str, "\n",&endstr);
    model->nres = strtoul(line, NULL, 10);

    // Resources data
    model->resources = calloc(model->nres, sizeof(resource_t *));
    if(!model->resources) goto ERROR;
    for (int64_t i=0; i< model->nres;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	model->resources[i] = resource_new(line);
    }

    // Number of machines
    line = strtok_r(NULL,"\n",&endstr);
    model->nmach = strtoul(line, NULL, 10);

    model->nlocations=1;
    // machines data
    model->machines = calloc(model->nmach, sizeof(machine_t *));
    for (int64_t i=0; i< model->nmach;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	model->machines[i] = machine_new(model->nres,
	  model->nmach,
	  line);
	if (model->nlocations < machine_location(model->machines[i])+1)
	    model->nlocations = machine_location(model->machines[i])+1;
    }
    int64_t *loc=calloc(model->nlocations,sizeof(int64_t));
    for (int64_t i=0; i< model->nmach;i++) {	
	loc[machine_location(model->machines[i])]++;
    }
    model->locations = calloc(model->nlocations,sizeof(int64_t *));
    for (int64_t l=0; l<model->nlocations;l++) {
	model->locations[l]=calloc(loc[l]+1,sizeof(int64_t));
	model->locations[l][loc[l]--]=-1;
    }
    for (int64_t i=0; i< model->nmach;i++) {
	
	model->locations[machine_location(model->machines[i])][loc[machine_location(model->machines[i])]--]=i;
    }
    free(loc);
    
    
    // Number of services
    line = strtok_r(NULL,"\n",&endstr);
    model->nserv = strtoul(line, NULL, 10);

    // services data
    model->services = calloc(model->nmach, sizeof(service_t *));
    for (int64_t i=0; i< model->nserv;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	model->services[i] = service_new(line);
    }

    // Number of process
    line = strtok_r(NULL, "\n", &endstr);
    model->nproc = strtoul ( line, NULL, 10);

    // process data
    model->processes = calloc(model->nproc, sizeof(process_t *));
    for (int64_t i =0; i< model->nproc; i++) {
	line = strtok_r(NULL, "\n", &endstr);
	model->processes[i] = process_new(model->nres,line);
    }

    // Number of balance costs
    line = strtok_r(NULL, "\n", &endstr);
    model->nbalance = strtoul ( line, NULL, 10);
    
    // balance data
    model->balance = calloc(model->nbalance, sizeof(balance_t *));
    for (int64_t i =0; i< model->nbalance; i++) {
	line = strtok_r(NULL, "\n", &endstr);
	model->balance[i] = balance_new(line);
	line = strtok_r(NULL, "\n", &endstr);
	balance_set_cost(model->balance[i], line);
    }

    // Moves costs

    line = strtok_r(NULL, "\n", &endstr);

    char *tok, *endtok;
    tok = strtok_r(line, " ", &endtok);
    model->wpmc = strtoul(tok, NULL, 10);
    tok = strtok_r(NULL, " ", &endtok);
    model->wsmc = strtoul(tok, NULL, 10);
    tok = strtok_r(NULL, " ", &endtok);
    model->wmmc = strtoul(tok, NULL, 10);

    

    return model;
 ERROR:
    model_destroy(&model);
    return NULL;
}

void
model_destroy(model_t** self_p){
    assert(self_p);
    if(!*self_p) return;

    model_t *self = *self_p;

    for (int64_t i =0 ; i< self->nbalance; i++) {
	balance_destroy(&(self->balance[i]));
    }
    if(self->balance) free(self->balance);
    
    for (int64_t i =0 ; i< self->nproc; i++) {
	process_destroy(&(self->processes[i]));
    }
    if(self->processes) free(self->processes);

    for (int64_t i =0 ; i< self->nserv; i++) {
	service_destroy(&(self->services[i]));
    }
    if(self->services) free(self->services);

    for (int64_t i =0 ; i< self->nmach; i++) {
	machine_destroy(&(self->machines[i]));
    }
    if (self->machines) free(self->machines);

    for (int64_t i =0 ; i< self->nres; i++) {
	resource_destroy(&(self->resources[i]));
    }
    if(self->resources) free(self->resources);

    for (int64_t i =0 ; i< self->nlocations; i++) {
	free(self->locations[i]);
    }
    if(self->locations) free(self->locations);
    
    free(self);
    *self_p=NULL;
}

int64_t
model_nres(model_t *self) {
	assert(self);

	return self->nres;
}


int64_t
model_nmach(model_t *self) {
	assert(self);

	return self->nmach;
}


int64_t
model_nserv(model_t *self) {
	assert(self);

	return self->nserv;
}


int64_t
model_nproc(model_t *self) {
	assert(self);

	return self->nproc;
}


int64_t
model_nbalance(model_t *self) {
	assert(self);

	return self->nbalance;
}


int64_t
model_wpmc(model_t *self) {
	assert(self);

	return self->wpmc;
}


int64_t
model_wsmc(model_t *self) {
	assert(self);

	return self->wsmc;
}


int64_t
model_wmmc(model_t *self) {
	assert(self);

	return self->wmmc;
}

resource_t *
model_resource(model_t *self, int64_t idx) {
    assert(self);
    assert(idx < self->nres);

    return self->resources[idx];
}

machine_t *
model_machine(model_t *self, int64_t idx) {
    assert(self);
    assert(idx < self->nmach);

    return self->machines[idx];
}

process_t *
model_process(model_t *self, int64_t idx) {
    assert(self);
    assert(idx < self->nproc);

    return self->processes[idx];
}


balance_t *
model_balance(model_t *self, int64_t idx) {
    assert(self);
    assert(idx < self->nproc);

    return self->balance[idx];
}

// TODO
static bool
s_model_validate_capacity(model_t *self, int64_t sz, int64_t *procmap){
    if (self){}
    if (sz){}
    if (procmap){}
    return true;
}

// TODO
static bool
s_model_validate_conflict(model_t *self, int64_t sz, int64_t *procmap){
    if (self){}
    if (sz){}
    if (procmap){}
    return true;
}

// TODO
static bool
s_model_validate_spread(model_t *self, int64_t sz, int64_t *procmap){
    if (self){}
    if (sz) {}
    if (procmap){}
    return true;
}

// TODO
static bool
s_model_validate_dependency(model_t *self, int64_t sz, int64_t *procmap){
    if (self){}
    if (sz){}
    if (procmap){}
    return true;
}

// TODO
static bool
s_model_validate_transient(model_t *self, int64_t sz, int64_t *procmap){
    if (self){}
    if (sz){}
    if (procmap){}
    return true;
}

// TODO
bool
model_validate(model_t *self, int64_t sz, int64_t *procmap){
    assert(self);
    assert(procmap);
    assert(sz == self->nproc);

    if(!s_model_validate_capacity(self,sz,procmap))
	return false;

    if(!s_model_validate_conflict(self,sz,procmap))
	return false;

    if(!s_model_validate_spread(self,sz,procmap))
	return false;

    if(!s_model_validate_dependency(self,sz,procmap))
	return false;

    if(!s_model_validate_transient(self,sz,procmap))
	return false;

    return true;
    
}	

/*
static int64_t **
s_model_u_new(int64_t nmach, int64_t nres) {
    int64_t **u = calloc(nmach,sizeof(int64_t*));
    for (int64_t i=0 ; i < nmach; i++)
	u[i] = calloc(nres,sizeof(int64_t));
    return u;
}

static void 
s_model_u_free(int64_t **u,int64_t nmach) {
    for (int64_t i=0 ; i < nmach; i++)
	free(u[i]);
    free(u);
}
*/

// TODO
int64_t
model_calculate(model_t *self, int64_t procmapsz, int64_t *procmap, int64_t pobjsz, int64_t *pobj){
    assert(self);
    assert(procmap);
    assert(procmapsz == self->nproc);
    assert(!pobjsz || pobjsz == 5);
    assert((!pobjsz && !pobj) || (pobjsz && pobj));

    int64_t loadcost,balancecost,pmc,smc,mmc;
   
    loadcost = balancecost = pmc = smc = mmc = 0;

    int64_t u[model_nmach(self)][model_nres(self)];
    int64_t a[model_nmach(self)][model_nres(self)];
    
    for (int64_t i =0; i< model_nmach(self);i++) {
	for(int64_t j =0; j<model_nres(self); j++) {
	    u[i][j] = 0;
	}
    }
    
    for (int64_t pidx=0; pidx < model_nproc(self); pidx++) {
	int64_t midx = procmap[pidx];
	process_t *p = model_process(self,pidx);
	for(int64_t ridx =0; ridx<model_nres(self); ridx++) {
	    u[midx][ridx] += process_requirement(p,ridx);
	}
    }

    for(int64_t ridx =0; ridx<model_nres(self); ridx++) {
	resource_t *r = model_resource(self,ridx);
	int64_t rloadcost = 0;
	for (int64_t midx =0; midx< model_nmach(self);midx++) {
	    int64_t c = 0;
	    c = u[midx][ridx] - machine_safecap(model_machine(self,midx),ridx);

	    if (c<0) c=0;
	    rloadcost += c;
	}
	loadcost += resource_loadcost(r)*rloadcost;
    }

    for (int64_t i =0; i< model_nmach(self);i++) {
	for(int64_t j =0; j<model_nres(self); j++) {
	    a[i][j] = machine_cap(model_machine(self,i),j) - u[i][j];
	}
    }

    for (int64_t k = 0; k<model_nbalance(self);k++) {
	balance_t *b = model_balance(self,k);
	for (int64_t i =0; i< model_nmach(self);i++) {
	    int64_t p = (int64_t)balance_target(b) * a[i][balance_resource1(b)] - a[i][balance_resource2(b)];
	    if (p<0) p=0;
	    balancecost += balance_weightcost(b) * p;
	}
    }

    
    if(pobjsz) {
	pobj[0] = loadcost;
	pobj[1] = balancecost;
	pobj[2] = pmc;
	pobj[3] = smc;
	pobj[4] = mmc;
    }

    return loadcost +
      balancecost +
      model_wpmc(self)*pmc +
      model_wsmc(self)*smc +
      model_wmmc(self)*mmc +
      0;
}

void
model_test(bool verbose) {



    if (verbose) printf("  * model: ");
    char *line;
    model_t *inst;

    line = strdup(MODEL_EXAMPLE1);
    inst = model_new(line);
    free(line);
    assert(inst);
    model_destroy(&inst);
    assert(!inst);
    model_destroy(&inst);

    line = strdup(MODEL_EXAMPLE1);
    inst = model_new(line);
    free(line);

    assert(2 == model_nres(inst));
    assert(model_resource(inst,1));
    assert(10==resource_loadcost(model_resource(inst,1)));
    assert(4 == model_nmach(inst));
    assert(2 == model_nserv(inst));
    assert(3 == model_nproc(inst));
    assert(1 == model_nbalance(inst));
    assert(1 == model_wpmc(inst));
    assert(10 == model_wsmc(inst));
    assert(100 == model_wmmc(inst));

    assert(model_validate(inst,3,(int64_t[3]){0, 3, 0}));

    assert(4200==model_calculate(inst,3,(int64_t[3]){0, 3, 0},0,NULL));
    assert(model_validate(inst,3,(int64_t[3]){0, 2, 0}));
    assert(model_validate(inst,3,(int64_t[3]){0, 2, 1}));
    
    model_destroy(&inst);

    line = strdup(MODEL_EXAMPLE2);
    inst = model_new(line);
    free(line);
    assert(0 == model_nbalance(inst));
    assert(1 == model_wpmc(inst));
    assert(10 == model_wsmc(inst));
    assert(100 == model_wmmc(inst));
    
    model_destroy(&inst);
    
    if(verbose) printf("OK\n");
    
}
