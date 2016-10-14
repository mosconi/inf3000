#include "roadef.h"

struct resource_t {
    bool transient;
    uint64_t wlc;
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
