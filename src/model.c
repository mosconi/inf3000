#include "roadef.h"
    
struct model_t {
    size_t nres;
    resource_t **resources;
    uint64_t *wlc;
    size_t nmach;
    machine_t **machines;
    size_t nserv;
    service_t **services;
    size_t nproc;
    process_t **processes;
    size_t nbalance;
    balance_t **balance;
    uint64_t wpmc;
    uint64_t wsmc;
    uint64_t wmmc;
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

    char filecontent[fsz+1];

    if(1!=fread(filecontent, fsz, 1, fp)) {
	fclose(fp);
	return NULL;
    }

    fclose(fp);
    return model_new(filecontent);
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
    for (uint64_t i=0; i< model->nres;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	model->resources[i] = resource_new(line);
    }

    // Number of machines
    line = strtok_r(NULL,"\n",&endstr);
    model->nmach = strtoul(line, NULL, 10);

    // machines data
    model->machines = calloc(model->nmach, sizeof(machine_t *));
    for (uint64_t i=0; i< model->nmach;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	model->machines[i] = machine_new(model->nres,
					    model->nmach,
					    line);
    }

    // Number of services
    line = strtok_r(NULL,"\n",&endstr);
    model->nserv = strtoul(line, NULL, 10);

    // services data
    model->services = calloc(model->nmach, sizeof(service_t *));
    for (uint64_t i=0; i< model->nserv;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	model->services[i] = service_new(line);
    }

    // Number of process
    line = strtok_r(NULL, "\n", &endstr);
    model->nproc = strtoul ( line, NULL, 10);

    // process data
    model->processes = calloc(model->nproc, sizeof(process_t *));
    for (uint64_t i =0; i< model->nproc; i++) {
	line = strtok_r(NULL, "\n", &endstr);
	model->processes[i] = process_new(model->nres,line);
    }

    // Number of balance costs
    line = strtok_r(NULL, "\n", &endstr);
    model->nbalance = strtoul ( line, NULL, 10);
    
    // balance data
    model->balance = calloc(model->nbalance, sizeof(balance_t *));
    for (uint64_t i =0; i< model->nbalance; i++) {
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

    for (size_t i =0 ; i< self->nbalance; i++) {
	balance_destroy(&(self->balance[i]));
    }
    if(self->balance) free(self->balance);
    
    for (size_t i =0 ; i< self->nproc; i++) {
	process_destroy(&(self->processes[i]));
    }
    if(self->processes) free(self->processes);

    for (size_t i =0 ; i< self->nserv; i++) {
	service_destroy(&(self->services[i]));
    }
    if(self->services) free(self->services);

    for (size_t i =0 ; i< self->nmach; i++) {
	machine_destroy(&(self->machines[i]));
    }
    if (self->machines) free(self->machines);

    for (size_t i =0 ; i< self->nres; i++) {
	resource_destroy(&(self->resources[i]));
    }
    if(self->resources) free(self->resources);
    
    free(self);
    *self_p=NULL;
}

uint64_t
model_nres(model_t *self) {
	assert(self);

	return self->nres;
}


uint64_t
model_nmach(model_t *self) {
	assert(self);

	return self->nmach;
}


uint64_t
model_nserv(model_t *self) {
	assert(self);

	return self->nserv;
}


uint64_t
model_nproc(model_t *self) {
	assert(self);

	return self->nproc;
}


uint64_t
model_nbalance(model_t *self) {
	assert(self);

	return self->nbalance;
}


uint64_t
model_wpmc(model_t *self) {
	assert(self);

	return self->wpmc;
}


uint64_t
model_wsmc(model_t *self) {
	assert(self);

	return self->wsmc;
}


uint64_t
model_wmmc(model_t *self) {
	assert(self);

	return self->wmmc;
}



void
model_test(bool verbose) {

#define MODEL_EXAMPLE1 ""			\
	"2\n"					\
	"1 100\n"				\
	"0 100\n"				\
	"4\n"					\
	"0 0 30 400 16 80 0 1 4 5\n"		\
	"0 0 10 240 8 160 1 0 3 4\n"		\
	"1 1 15 100 12 80 4 3 0 2\n"		\
	"1 2 10 100 8 80 5 4 2 0\n"		\
	"2\n"					\
	"2 0\n"					\
	"1 1 0\n"				\
	"3\n"					\
	"0 12 10 1000\n"			\
	"0 10 20 100\n"				\
	"1 16 200 1\n"				\
	"1\n"					\
	"0 1 20\n"				\
	"10\n"					\
	"1 10 100\n"

#define MODEL_EXAMPLE2 ""			\
	"2\n"					\
	"1 100\n"				\
	"0 100\n"				\
	"4\n"					\
	"0 0 30 400 16 80 0 1 4 5\n"		\
	"0 0 10 240 8 160 1 0 3 4\n"		\
	"1 1 15 100 12 80 4 3 0 2\n"		\
	"1 2 10 100 8 80 5 4 2 0\n"		\
	"2\n"					\
	"2 0\n"					\
	"1 1 0\n"				\
	"3\n"					\
	"0 12 10 1000\n"			\
	"0 10 20 100\n"				\
	"1 16 200 1\n"				\
	"0\n"					\
	"1 10 100\n"


    if (verbose) printf("  * model: ");
    char *line;
    model_t *inst;

    line = strdup(MODEL_EXAMPLE1);
    inst = model_new_string(line);
    free(line);
    assert(inst);
    model_destroy(&inst);
    assert(!inst);
    model_destroy(&inst);

    line = strdup(MODEL_EXAMPLE1);
    inst = model_new_string(line);
    free(line);

    assert(2 == model_nres(inst));
    assert(4 == model_nmach(inst));
    assert(2 == model_nserv(inst));
    assert(3 == model_nproc(inst));
    assert(1 == model_nbalance(inst));
    assert(1 == model_wpmc(inst));
    assert(10 == model_wsmc(inst));
    assert(100 == model_wmmc(inst));
    
    model_destroy(&inst);

    line = strdup(MODEL_EXAMPLE2);
    inst = model_new_string(line);
    free(line);
    assert(0 == model_nbalance(inst));
    assert(1 == model_wpmc(inst));
    assert(10 == model_wsmc(inst));
    assert(100 == model_wmmc(inst));
    
    model_destroy(&inst);

    
    if(verbose) printf("OK\n");
    
}
