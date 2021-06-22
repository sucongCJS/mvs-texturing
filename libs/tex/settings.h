#ifndef TEX_SETTINGS_HEADER
#define TEX_SETTINGS_HEADER

#include <sstream>
#include <vector>
#include <string>

#include "defines.h"

TEX_NAMESPACE_BEGIN

enum DataTerm {
    DATA_TERM_AREA = 0,
    DATA_TERM_GMI = 1
};

enum SmoothnessTerm {
    SMOOTHNESS_TERM_POTTS = 0
};

enum OutlierRemoval {
    OUTLIER_REMOVAL_NONE = 0,
    OUTLIER_REMOVAL_GAUSS_DAMPING = 1,
    OUTLIER_REMOVAL_GAUSS_CLAMPING = 2
};

/** Enum representing tone mapping choice. 色调映射 */
enum ToneMapping {
    TONE_MAPPING_NONE = 0,
    TONE_MAPPING_GAMMA = 1
};

struct Settings{

    DataTerm data_term = DATA_TERM_GMI;
    SmoothnessTerm smoothness_term = SMOOTHNESS_TERM_POTTS;
    OutlierRemoval outlier_removal = OUTLIER_REMOVAL_NONE;
    ToneMapping tone_mapping = TONE_MAPPING_NONE;

    bool geometric_visibility_test = true;
    bool global_seam_leveling = true;
    bool local_seam_leveling = true;
};

TEX_NAMESPACE_END

#endif /* TEX_SETTINGS_HEADER */