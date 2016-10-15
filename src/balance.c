#include "roadef.h"

struct balance_t {
    uint64_t resource1;
    uint64_t resource2;
    uint64_t target;
    uint64_t wbc;//weight for balance cost
};

balance_t *
balance_new(char *line){
    assert(line);
    assert(strneq("",line));

    balance_t *b = calloc(1,sizeof(balance_t));
    if (!b) return NULL;
    
    char *tok, *endtok;
    tok = strtok_r(line," ", &endtok);
    b->resource1 = strtoul(tok,NULL, 10);
    tok = strtok_r(NULL, " ", &endtok);
    b->resource2 = strtoul(tok,NULL, 10);
    tok = strtok_r(NULL, " ", &endtok);
    b->target = strtoul(tok,NULL, 10);


    return b;
}

void
balance_destroy(balance_t **self_p){
    assert(self_p);

    if(!*self_p) return ;

    balance_t *self = *self_p;

    free(self);
    
    *self_p=NULL;
}

void
balance_set_cost(balance_t *self, char *line){
    assert(self);
    assert(line);
    assert(strneq("",line));

    self->wbc =  strtoul(line,NULL, 10);

}
    

int
balance_cmp(balance_t *self, uint64_t k1, uint64_t k2){
    assert(self);

    if (self->resource1 < k1) return -1;
    if (self->resource1 > k1) return 1;

    if (self->resource2 < k2) return -1;
    if (self->resource2 > k2) return 1;

    return 0;
}

uint64_t
balance_target(balance_t* self) {
    assert(self);

    return self->target;
}

uint64_t
balance_weightcost(balance_t* self) {
    assert(self);

    return self->wbc;
}

uint64_t
balance_resource1(balance_t *self) {
    assert(self);
    return self->resource1;
}
uint64_t
balance_resource2(balance_t *self) {
    assert(self);
    return self->resource2;
}

void
balance_test(bool verbose){

    char *line;
    balance_t *b;

    line = strdup("0 1 20");
    b = balance_new(line);
    free(line);
    assert(b);
    balance_destroy(&b);
    assert(!b);
    balance_destroy(&b);

    line = strdup("0 1 20");
    b = balance_new(line);
    free(line);
    assert(0 == balance_resource1(b));
    assert(1 == balance_resource2(b));
    assert(20 == balance_target(b));

    assert(balance_cmp(b, 0, 0) > 0);
    assert(balance_cmp(b, 0, 2) < 0);
    assert(balance_cmp(b, 1, 0) < 0);
    assert(balance_cmp(b, 0, 1) ==0);

    balance_set_cost(b,"10");
    assert(10 == balance_weightcost(b));
    balance_destroy(&b);
    

}
	
