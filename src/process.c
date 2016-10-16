#include "roadef.h"

struct process_t {
    int64_t service;
    int64_t nres;
    int64_t *req;
    int64_t pmc;
};

process_t *
process_new(int64_t nres, char *line) {
    assert(line);
    assert(strneq("",line));

    process_t *p = calloc (1,sizeof(process_t));

    char *tok,*endtok;

    tok = strtok_r(line," ",&endtok);
    p->service = strtoul(tok, NULL, 10);
    
    p->nres = nres;
    p->req = calloc(nres, sizeof(int64_t));

    if (!p->req) {
	process_destroy(&p);
	return NULL;
    }

    for (int64_t i=0; i < nres; i++) {
	tok = strtok_r(NULL," ",&endtok);
	p->req[i] = strtoul(tok, NULL, 10);
    }
    
    tok = strtok_r(NULL," ",&endtok);
    p->pmc = strtoul(tok, NULL, 10);

    return p;
}

void
process_destroy(process_t **self_p) {
    assert(self_p);

    if (!*self_p) return;

    process_t *p = *self_p;

    if(p->req) free(p->req);

    free(p);
    *self_p = NULL;
}

int64_t
process_service(process_t *self) {
    assert(self);

    return self->service;
}

int64_t
process_movecost(process_t *self) {
    assert(self);

    return self->pmc;
}

int64_t
process_requirement(process_t *self, int64_t idx) {
    assert(self);
    assert(idx < self->nres);

    return self->req[idx];
}

void
process_test(bool verbose) {

    if (verbose) printf("  * process: ");
    char *line = strdup("0 12 10 1000");
    process_t *p = process_new(2,line);
    free(line);
    assert(p);
    process_destroy(&p);
    assert(!p);
    process_destroy(&p);
    line = strdup("0 12 10 1000");
    p = process_new(2,line);
    free(line);
    assert(0==process_service(p));
    assert(12==process_requirement(p,0));
    assert(10==process_requirement(p,1));
    assert(1000==process_movecost(p));
    process_destroy(&p);

    if(verbose) printf("OK\n");
}
    
