#include "sarlacc.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(count_gaps_by_base, 3),
    REGISTER(count_gaps_by_align, 3),

    REGISTER(mask_bad_bases, 3),    
    REGISTER(unmask_bases, 2),

    REGISTER(create_consensus_basic, 3), 
    REGISTER(create_consensus_basic_loop, 3), 
    REGISTER(create_consensus_quality, 3), 
    REGISTER(create_consensus_quality_loop, 3), 

    REGISTER(umi_group, 2),
    REGISTER(descending_graph_cluster, 1),
    REGISTER(compute_lev_masked, 1),

    REGISTER(find_homopolymers, 1),
    REGISTER(match_homopolymers, 2),
    REGISTER(find_errors, 2),

    REGISTER(get_aligned_sequence, 5),
    {NULL, NULL, 0}
};

void attribute_visible R_init_sarlacc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

