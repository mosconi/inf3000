#include "roadef.h"
    
struct instance_t {
    size_t nres;
    resource_t **resources;
    int64_t *wlc;
    size_t nmach;
    machine_t **machines;
	
    size_t nproc;
	
    size_t nserv;
    size_t nloc;
    size_t nneigh;
    size_t nbalance;
    balance_t **balance;
    int64_t wpmc;
    int64_t wsmc;
    int64_t wmmc;
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
    for (int i=0; i< instance->nres;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	instance->resources[i] = resource_new(line);
    }

    // Number of machines
    line = strtok_r(NULL,"\n",&endstr);
    instance->nmach = strtoul(line, NULL, 10);

    // machines data
    instance->machines = calloc(instance->nmach, sizeof(machine_t *));
    for (int i=0; i< instance->nmach;i++) {
	line = strtok_r(NULL,"\n",&endstr);
	instance->machines[i] = machine_new(instance->nres,
					    instance->nmach,
					    line);
    }

    // Number of services


    // services data

    // Number of process

    // process data

    // Number of balance costs
    
    // balance data

    // Moves costs

    return instance;
 ERROR:
    instance_destroy(&instance);
    return NULL;
}

void
instance_destroy(instance_t** self_p){
    assert(self_p);
    if(!self_p) return;

    instance_t *self = *self_p;


    free(self);
    *self_p=NULL;
}
