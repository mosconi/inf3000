#include "roadef.h"


struct service_t {
    uint64_t spread;
    size_t ndep;
    uint64_t *deps;
};

service_t *
service_new(char *line) {
    assert(line);
    assert(strneq("",line));

    service_t *s = calloc(1,sizeof(service_t));
    if (!s) return NULL;

    char *endtok, *tok;
    tok = strtok_r(line, " " , &endtok);
    s->spread = strtoul(tok, NULL, 10);

    tok = strtok_r(NULL, " " , &endtok);
    s->ndep = strtoul(tok, NULL, 10);

    s->deps = calloc(s->ndep,sizeof(uint64_t));

    for (uint64_t i = 0; i<s->ndep; i++) {
	tok = strtok_r(NULL, " " , &endtok);
	s->deps[i] = strtoul(tok, NULL, 10);
    }

    return s;
}

void
service_destroy(service_t **self_p) {
    assert(self_p);

    if(!*self_p) return;

    service_t *s = *self_p;

    if (s->deps) free(s->deps);
    
    free(s);
    *self_p = NULL;
}

uint64_t
service_spread(service_t *self){
    assert(self);

    return self->spread;
}

uint64_t
service_ndeps(service_t *self){
    assert(self);

    return self->ndep;
}

uint64_t
service_dep(service_t *self, uint64_t idx){
    assert(self);
    assert(idx < self->ndep);

    return self->deps[idx];
}


     
void
service_test(bool verbose){

    if (verbose) printf("  * service: ");
    char *line;
    service_t *s;

    line = strdup("2 0");
    s = service_new(line);
    free(line);
    assert(s);
    service_destroy(&s);
    assert(!s);
    service_destroy(&s);

    line = strdup("2 0");
    s = service_new(line);
    free(line);

    assert(2 == service_spread(s));
    assert(0 == service_ndeps(s));

    service_destroy(&s);

    line = strdup("1 1 0");
    s = service_new(line);
    free(line);

    assert(1 == service_spread(s));
    assert(1 == service_ndeps(s));
    assert(0 == service_dep(s, 0));

    service_destroy(&s);

    if (verbose) printf("OK\n");
}
