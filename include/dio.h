
#ifndef WORLD_DIO_H_
#define WORLD_DIO_H_

#include "common.h"

BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Struct for DIO
//-----------------------------------------------------------------------------
typedef struct {
    float f0_floor;
    float f0_ceil;
    float channels_in_octave;
    float frame_period; // msec
    int speed; // (1, 2, ..., 12)
    float allowed_range; // Threshold used for fixing the F0 contour.
} DioOption;

//-----------------------------------------------------------------------------
// DIO
//
// Input:
//   x                    : Input signal
//   x_length             : Length of x
//   fs                   : Sampling frequency
//   option               : Struct to order the parameter for DIO
//
// Output:
//   temporal_positions   : Temporal positions.
//   f0                   : F0 contour.
//-----------------------------------------------------------------------------
void Dio(const float* x, int x_length, int fs, const DioOption* option,
    float* temporal_positions, float* f0);

//-----------------------------------------------------------------------------
// InitializeDioOption allocates the memory to the struct and sets the
// default parameters.
//
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
void InitializeDioOption(DioOption* option);

//-----------------------------------------------------------------------------
// GetSamplesForDIO() calculates the number of samples required for Dio().
//
// Input:
//   fs             : Sampling frequency [Hz]
//   x_length       : Length of the input signal [Sample].
//   frame_period   : Frame shift [msec]
//
// Output:
//   The number of samples required to store the results of Dio()
//-----------------------------------------------------------------------------
int GetSamplesForDIO(int fs, int x_length, float frame_period);

//-----------------------------------------------------------------------------
// StoneMask() refines the estimated F0 by Dio()
//
// Input:
//   x                      : Input signal
//   x_length               : Length of the input signal
//   fs                     : Sampling frequency
//   time_axis              : Temporal information
//   f0                     : f0 contour
//   f0_length              : Length of f0
//
// Output:
//   refined_f0             : Refined F0
//-----------------------------------------------------------------------------
void StoneMask(const float* x, int x_length, int fs,
    const float* temporal_positions, const float* f0, int f0_length,
    float* refined_f0);

END_C_DECLS

#endif // WORLD_DIO_H_
