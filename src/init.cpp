#include "sarlacc.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(adaptor_align, 8),
    REGISTER(adaptor_align_score_only, 6),
    REGISTER(barcode_align, 6),
    REGISTER(general_align, 7),

    REGISTER(mask_bad_bases, 4),
    REGISTER(unmask_alignment, 2),

    REGISTER(create_consensus_basic, 3), 
    REGISTER(create_consensus_basic_loop, 3), 
    REGISTER(create_consensus_quality, 4), 
    REGISTER(create_consensus_quality_loop, 4), 

    REGISTER(umi_group, 5),
    REGISTER(fast_levdist_test, 3),
    REGISTER(cluster_umis_test, 1), 
    REGISTER(compute_lev_masked, 1),

    REGISTER(find_homopolymers, 1),
    REGISTER(match_homopolymers, 2),
    REGISTER(find_errors, 2),

    REGISTER(quick_msa, 7),

    {NULL, NULL, 0}
};

void attribute_visible R_init_sarlacc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
