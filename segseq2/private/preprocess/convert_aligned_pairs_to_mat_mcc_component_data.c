/*
 * MATLAB Compiler: 4.7 (R2007b)
 * Date: Mon May 11 10:13:49 2009
 * Arguments: "-B" "macro_default" "-v" "-m" "-W" "main" "-T" "link:exe" "-w"
 * "enable" "convert_aligned_pairs_to_mat.m" 
 */

#include "mclmcr.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_convert_aligned_pairs_to_mat_session_key[] = {
    '3', 'F', 'E', '5', '7', '6', 'C', '5', 'C', '8', 'C', 'A', '5', '0', 'B',
    '4', '1', '2', '2', '0', 'C', 'D', '3', '2', '1', '5', 'C', '0', '9', '3',
    'C', 'E', '3', 'A', '4', '2', '2', 'A', '2', 'B', '2', 'E', '5', '1', 'B',
    'F', '0', 'D', '2', '0', 'B', 'C', '8', '1', 'F', 'C', 'F', '8', 'C', 'E',
    '2', '3', 'B', 'A', '8', '8', 'F', 'F', 'C', '3', '7', 'D', 'C', '7', 'F',
    '1', '5', '2', '9', '8', '5', '2', '1', '4', 'C', 'A', 'A', '3', '7', '3',
    '7', '4', '7', '3', '4', '8', '5', '4', 'A', '2', '8', 'F', '2', 'A', '5',
    'D', 'F', 'B', 'D', '0', 'F', '7', '3', 'C', 'C', '2', 'C', '1', '9', '8',
    '8', 'A', 'F', 'A', 'F', 'B', 'D', '1', '9', 'D', '5', '2', '4', '7', 'E',
    '5', '9', '8', 'F', '7', '2', '6', 'F', 'A', '3', '8', '3', '4', 'B', 'A',
    'D', '9', '9', '5', '7', '0', 'F', '0', '5', '8', '3', '9', '7', 'A', 'D',
    '0', '5', '0', '2', '7', '7', '2', '6', 'D', '2', 'D', 'F', 'F', '9', 'E',
    '5', '4', '3', 'C', 'B', '0', '5', 'C', 'E', '1', '3', '3', '3', 'A', '9',
    '3', 'B', '8', '9', '0', 'E', '0', '9', '5', '4', 'A', '1', 'A', '0', 'D',
    '4', '2', '1', '1', '4', 'E', '4', '2', '6', '1', 'B', '1', '4', 'F', '5',
    '8', 'B', 'B', 'D', 'F', 'F', '5', '3', '6', '6', 'B', '9', '6', '8', '4',
    '1', 'C', '2', '3', '2', '6', '1', '7', '5', '1', '0', '2', '6', '1', '2',
    'A', '\0'};

const unsigned char __MCC_convert_aligned_pairs_to_mat_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_convert_aligned_pairs_to_mat_matlabpath_data[] = 
  { "convert_aligned_pairs_to_mat/", "toolbox/compiler/deploy/",
    "$TOOLBOXMATLABDIR/general/", "$TOOLBOXMATLABDIR/ops/",
    "$TOOLBOXMATLABDIR/lang/", "$TOOLBOXMATLABDIR/elmat/",
    "$TOOLBOXMATLABDIR/elfun/", "$TOOLBOXMATLABDIR/specfun/",
    "$TOOLBOXMATLABDIR/matfun/", "$TOOLBOXMATLABDIR/datafun/",
    "$TOOLBOXMATLABDIR/polyfun/", "$TOOLBOXMATLABDIR/funfun/",
    "$TOOLBOXMATLABDIR/sparfun/", "$TOOLBOXMATLABDIR/scribe/",
    "$TOOLBOXMATLABDIR/graph2d/", "$TOOLBOXMATLABDIR/graph3d/",
    "$TOOLBOXMATLABDIR/specgraph/", "$TOOLBOXMATLABDIR/graphics/",
    "$TOOLBOXMATLABDIR/uitools/", "$TOOLBOXMATLABDIR/strfun/",
    "$TOOLBOXMATLABDIR/imagesci/", "$TOOLBOXMATLABDIR/iofun/",
    "$TOOLBOXMATLABDIR/audiovideo/", "$TOOLBOXMATLABDIR/timefun/",
    "$TOOLBOXMATLABDIR/datatypes/", "$TOOLBOXMATLABDIR/verctrl/",
    "$TOOLBOXMATLABDIR/codetools/", "$TOOLBOXMATLABDIR/helptools/",
    "$TOOLBOXMATLABDIR/demos/", "$TOOLBOXMATLABDIR/timeseries/",
    "$TOOLBOXMATLABDIR/hds/", "$TOOLBOXMATLABDIR/guide/",
    "$TOOLBOXMATLABDIR/plottools/", "toolbox/local/" };

static const char * MCC_convert_aligned_pairs_to_mat_classpath_data[] = 
  { "" };

static const char * MCC_convert_aligned_pairs_to_mat_libpath_data[] = 
  { "" };

static const char * MCC_convert_aligned_pairs_to_mat_app_opts_data[] = 
  { "" };

static const char * MCC_convert_aligned_pairs_to_mat_run_opts_data[] = 
  { "" };

static const char * MCC_convert_aligned_pairs_to_mat_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_convert_aligned_pairs_to_mat_component_data = { 

  /* Public key data */
  __MCC_convert_aligned_pairs_to_mat_public_key,

  /* Component name */
  "convert_aligned_pairs_to_mat",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_convert_aligned_pairs_to_mat_session_key,

  /* Component's MATLAB Path */
  MCC_convert_aligned_pairs_to_mat_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  34,

  /* Component's Java class path */
  MCC_convert_aligned_pairs_to_mat_classpath_data,
  /* Number of directories in the Java class path */
  0,

  /* Component's load library path (for extra shared libraries) */
  MCC_convert_aligned_pairs_to_mat_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_convert_aligned_pairs_to_mat_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_convert_aligned_pairs_to_mat_run_opts_data,
  /* Number of MCR global runtime options */
  0,
  
  /* Component preferences directory */
  "convert_aligned_pairs_to_mat_058E0D32F78E2CA222940AE4C89C8220",

  /* MCR warning status data */
  MCC_convert_aligned_pairs_to_mat_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


