#include "roadef.h"

struct resource_t {
    bool transient;
    int64_t wlc;
};

resource_t *
resource_new(char *line){
    assert(line);
    assert(strneq("",line));

    resource_t *r = calloc(1, sizeof(resource_t));

    char *endtok;
    char *tok = strtok_r(line," ",&endtok);
    if (streq("1",tok)){
	r->transient= true;
    } else if (streq("0",tok)){
	r->transient= false;
    } else {
	free(r);
	return NULL;
    }

    tok = strtok_r(NULL," ", &endtok);

    r->wlc = strtoul(tok,NULL,10);

    return r;
}

void
resource_destroy(resource_t **self_p) {
    assert(self_p);
    if (!*self_p) return;
    resource_t *self = *self_p;
    free(self);
    *self_p=NULL;
}

bool
resource_transient(resource_t *self) {
    assert(self);
    return self->transient;
}

int64_t
resource_loadcost(resource_t *self){
    assert(self);
    return self->wlc;
}

void
resource_test(bool verbose){

    if (verbose) printf("  * resource: ");
    #define RESOURCE_TEST1 "1 100"
    char *line = strdup(RESOURCE_TEST1);
    resource_t *r = resource_new(line);
    free(line);
    assert(r);
    resource_destroy(&r);
    assert(!r);
    resource_destroy(&r);
    line = strdup("1 100");
    r = resource_new(line);
    free(line);
    assert(resource_transient(r));
    assert(100 == resource_loadcost(r));
    resource_destroy(&r);
    line = strdup("0 100");
    r = resource_new(line);
    free(line);
    assert(!resource_transient(r));
    resource_destroy(&r);

    line = strdup("2 100");
    r = resource_new(line);
    free(line);
    assert(!r);
    
    if (verbose) printf("OK\n");
}
