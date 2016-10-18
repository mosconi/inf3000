#ifndef __MODEL_H__
#define __MODEL_H__

#ifdef __cplpusplus
extern "C" {
#endif

    typedef struct model_t model_t;

    model_t *
    model_new(char *);

    model_t *
    model_load(char *);

    void
    model_destroy(model_t **);

    int64_t
    model_nres(model_t *);

    int64_t
    model_nmach(model_t *);

    int64_t
    model_nserv(model_t *);

    int64_t
    model_nproc(model_t *);

    int64_t
    model_nbalance(model_t *);

    int64_t
    model_wpmc(model_t *);

    int64_t
    model_wsmc(model_t *);

    int64_t
    model_wmmc(model_t *);

    resource_t *
    model_resource(model_t*, int64_t);

    machine_t *
    model_machine(model_t*, int64_t);

    process_t *
    model_process(model_t*, int64_t);

    service_t *
    model_service(model_t*, int64_t);

    balance_t *
    model_balance(model_t*, int64_t);

    bool
    model_validate(model_t *, int64_t , int64_t *);

    int64_t
    model_calculate(model_t *, int64_t , int64_t *, int64_t, int64_t *);

    
    void
    model_test(bool);

    #define MODEL_EXAMPLE1 ""			\
	"2\n"					\
	"1 100\n"				\
	"0 10\n"				\
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
	"1 6 200 1\n"				\
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

    
#ifdef __cplpusplus
}
#endif

#endif
