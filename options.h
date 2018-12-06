#define DEF_LENGTH_THRESHOLD 500
#define DEF_GC_THRESHOLD 0.0
#define DEF_OBSEXP_THRESHOLD 0.0

void usage(char *prog);
int options(int argc, char *argv[]);

extern int length_threshold;
extern double gc_threshold;
extern double observed_expected_threshold;
extern int ace_output;
