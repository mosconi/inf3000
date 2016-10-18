#include "roadef.h"


int
main(int argc, char **argv) {

    int opt;

    char *modelfile=NULL;
    char *instance=NULL;

    while ((opt = getopt(argc, argv, "i:m:")) != -1 ){
	switch (opt) {
	case 'i':
	    if (instance) free(instance);
	    instance = strdup(optarg);
	    break;
	case 'm':
	    if (modelfile) free(modelfile);
	    modelfile = strdup(optarg);
	    break;
	default:
	    if (instance) free(instance);
	    if (modelfile) free(modelfile);
	    exit(-1);
	}
    }

    if (!modelfile || streq("",modelfile))
	exit(-2);
    
    model_t *model = model_load(modelfile);

    printf("    memoria (current): %zu B\n", getCurrentRSS() );
    printf("    memoria (Peak):    %zu B\n", getPeakRSS() );
    
    model_destroy(&model);
    if (instance) free(instance);
    if (modelfile) free(modelfile);
    
}
