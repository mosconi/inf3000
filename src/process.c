#include "roadef.h"

struct process_t {
    uint64_t service;
    size_t nres;
    uint64_t *req;
    uint64_t pmc;
};

process_t *
process_new(size_t nres, char *line) {
    assert(line);
    assert(strneq("",line));

    process_t *p = calloc (1,sizeof(process_t));

    char *tok,*endtok;

    tok = strtok_r(line," ",&endtok);
    p->service = strtoul(tok, NULL, 10);
    
    p->nres = nres;
    p->req = calloc(nres, sizeof(uint64_t));

    if (!p->req) {
	process_destroy(&p);
	return NULL;
    }

    for (int i=0; i< nres; i++) {
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
