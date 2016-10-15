#include "roadef.h"
    
struct instance_t {
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

instance_t *
instance_new(char *filename){
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
    return instance_new_string(filecontent);
}

instance_t *
instance_new_string(char *str){
    assert(str);
    assert(strneq("",str));

    instance_t *instance = calloc(1, sizeof(instance_t));
    char *endstr;

    // Number of resources
    char *line = strtok_r(str, "\n",&endstr);
    instance->nres = strtoul(line, NULL, 10);

    // Resources data
    instance->resources = calloc(instance->nres, sizeof(resource_t *));
    if(!instance->resources) goto ERROR;
    for (uint64_t i=0; i< instance->nres;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	instance->resources[i] = resource_new(line);
    }

    // Number of machines
    line = strtok_r(NULL,"\n",&endstr);
    instance->nmach = strtoul(line, NULL, 10);

    // machines data
    instance->machines = calloc(instance->nmach, sizeof(machine_t *));
    for (uint64_t i=0; i< instance->nmach;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	instance->machines[i] = machine_new(instance->nres,
					    instance->nmach,
					    line);
    }

    // Number of services
    line = strtok_r(NULL,"\n",&endstr);
    instance->nserv = strtoul(line, NULL, 10);

    // services data
    instance->services = calloc(instance->nmach, sizeof(service_t *));
    for (uint64_t i=0; i< instance->nserv;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	instance->services[i] = service_new(line);
    }

    // Number of process
    line = strtok_r(NULL, "\n", &endstr);
    instance->nproc = strtoul ( line, NULL, 10);

    // process data
    instance->processes = calloc(instance->nproc, sizeof(process_t *));
    for (uint64_t i =0; i< instance->nproc; i++) {
	line = strtok_r(NULL, "\n", &endstr);
	instance->processes[i] = process_new(instance->nres,line);
    }

    // Number of balance costs
    line = strtok_r(NULL, "\n", &endstr);
    instance->nbalance = strtoul ( line, NULL, 10);
    
    // balance data
    instance->balance = calloc(instance->nbalance, sizeof(balance_t *));
    for (uint64_t i =0; i< instance->nbalance; i++) {
	line = strtok_r(NULL, "\n", &endstr);
	instance->balance[i] = balance_new(line);
	line = strtok_r(NULL, "\n", &endstr);
	balance_set_cost(instance->balance[i], line);
    }

    // Moves costs

    line = strtok_r(NULL, "\n", &endstr);

    char *tok, *endtok;
    tok = strtok_r(line, " ", &endtok);
    instance->wpmc = strtoul(tok, NULL, 10);
    tok = strtok_r(NULL, " ", &endtok);
    instance->wsmc = strtoul(tok, NULL, 10);
    tok = strtok_r(NULL, " ", &endtok);
    instance->wmmc = strtoul(tok, NULL, 10);

    return instance;
 ERROR:
    instance_destroy(&instance);
    return NULL;
}

void
instance_destroy(instance_t** self_p){
    assert(self_p);
    if(!*self_p) return;

    instance_t *self = *self_p;

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
instance_nres(instance_t *self) {
	assert(self);

	return self->nres;
}


uint64_t
instance_nmach(instance_t *self) {
	assert(self);

	return self->nmach;
}


uint64_t
instance_nserv(instance_t *self) {
	assert(self);

	return self->nserv;
}


uint64_t
instance_nproc(instance_t *self) {
	assert(self);

	return self->nproc;
}


uint64_t
instance_nbalance(instance_t *self) {
	assert(self);

	return self->nbalance;
}


uint64_t
instance_wpmc(instance_t *self) {
	assert(self);

	return self->wpmc;
}


uint64_t
instance_wsmc(instance_t *self) {
	assert(self);

	return self->wsmc;
}


uint64_t
instance_wmmc(instance_t *self) {
	assert(self);

	return self->wmmc;
}



void
instance_test(bool verbose) {

#define INSTANCE_EXAMPLE1 ""			\
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

#define INSTANCE_EXAMPLE2 ""			\
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


    if (verbose) printf("  * instance: ");
    char *line;
    instance_t *inst;

    line = strdup(INSTANCE_EXAMPLE1);
    inst = instance_new_string(line);
    free(line);
    assert(inst);
    instance_destroy(&inst);
    assert(!inst);
    instance_destroy(&inst);

    line = strdup(INSTANCE_EXAMPLE1);
    inst = instance_new_string(line);
    free(line);

    assert(2 == instance_nres(inst));
    assert(4 == instance_nmach(inst));
    assert(2 == instance_nserv(inst));
    assert(3 == instance_nproc(inst));
    assert(1 == instance_nbalance(inst));
    assert(1 == instance_wpmc(inst));
    assert(10 == instance_wsmc(inst));
    assert(100 == instance_wmmc(inst));
    
    instance_destroy(&inst);

    line = strdup(INSTANCE_EXAMPLE2);
    inst = instance_new_string(line);
    free(line);
    
    instance_destroy(&inst);

    
    if(verbose) printf("OK\n");
    
}
