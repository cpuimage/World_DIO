//-----------------------------------------------------------------------------
// F0 estimation based on DIO (Distributed Inline-filter Operation).
//-----------------------------------------------------------------------------
#include "dio.h"

#include <math.h>

#include "common.h"
#include "matlabfunctions.h"

//-----------------------------------------------------------------------------
// struct for GetFourZeroCrossingIntervals()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
//-----------------------------------------------------------------------------
typedef struct {
    float* negative_interval_locations;
    float* negative_intervals;
    int number_of_negatives;
    float* positive_interval_locations;
    float* positive_intervals;
    int number_of_positives;
    float* peak_interval_locations;
    float* peak_intervals;
    int number_of_peaks;
    float* dip_interval_locations;
    float* dip_intervals;
    int number_of_dips;
} ZeroCrossings;

namespace {
//-----------------------------------------------------------------------------
// DesignLowCutFilter() calculates the coefficients the filter.
//-----------------------------------------------------------------------------
static void DesignLowCutFilter(int N, int fft_size, float* low_cut_filter)
{
    for (int i = 1; i <= N; ++i)
        low_cut_filter[i - 1] = 0.5f - 0.5f * cosf(i * 2.0f * kPi / (N + 1));
    for (int i = N; i < fft_size; ++i)
        low_cut_filter[i] = 0.0f;
    float sum_of_amplitude = 0.0;
    for (int i = 0; i < N; ++i)
        sum_of_amplitude += low_cut_filter[i];
    for (int i = 0; i < N; ++i)
        low_cut_filter[i] = -low_cut_filter[i] / sum_of_amplitude;
    for (int i = 0; i < (N - 1) / 2; ++i)
        low_cut_filter[fft_size - (N - 1) / 2 + i] = low_cut_filter[i];
    for (int i = 0; i < N; ++i)
        low_cut_filter[i] = low_cut_filter[i + (N - 1) / 2];
    low_cut_filter[0] += 1.0;
}

//-----------------------------------------------------------------------------
// GetSpectrumForEstimation() calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
//-----------------------------------------------------------------------------
static void GetSpectrumForEstimation(const float* x, int x_length,
    int y_length, float actual_fs, int fft_size, int decimation_ratio,
    cmplx* y_spectrum)
{
    float* y = new float[fft_size];

    // Initialization
    for (int i = 0; i < fft_size; ++i)
        y[i] = 0.0;

    // Downsampling
    if (decimation_ratio != 1)
        decimate(x, x_length, decimation_ratio, y);
    else
        for (int i = 0; i < x_length; ++i)
            y[i] = x[i];

    // Removal of the DC component (y = y - mean value of y)
    float mean_y = 0.0;
    for (int i = 0; i < y_length; ++i)
        mean_y += y[i];
    mean_y /= y_length;
    for (int i = 0; i < y_length; ++i)
        y[i] -= mean_y;
    for (int i = y_length; i < fft_size; ++i)
        y[i] = 0.0;

    stb_fft_real_plan* forwardFFT = stb_fft_real_plan_dft_1d(fft_size);
    stb_fft_r2c_exec(forwardFFT, y, y_spectrum);

    // Low cut filtering (from 0.1.4)
    int cutoff_in_sample = matlab_round(actual_fs / 50.0f); // Cutoff is 50.0 Hz
    DesignLowCutFilter(cutoff_in_sample * 2 + 1, fft_size, y);

    cmplx* filter_spectrum = new cmplx[fft_size];
    stb_fft_r2c_exec(forwardFFT, y, filter_spectrum);
    float tmp = 0;
    for (int i = 0; i <= fft_size / 2; ++i) {
        // Complex number multiplications.
        tmp = y_spectrum[i].real * filter_spectrum[i].real - y_spectrum[i].imag * filter_spectrum[i].imag;
        y_spectrum[i].imag = y_spectrum[i].real * filter_spectrum[i].imag + y_spectrum[i].imag * filter_spectrum[i].real;
        y_spectrum[i].real = tmp;
    }

    free(forwardFFT);
    delete[] y;
    delete[] filter_spectrum;
}

//-----------------------------------------------------------------------------
// GetBestF0Contour() calculates the best f0 contour based on scores of
// all candidates. The F0 with highest score is selected.
//-----------------------------------------------------------------------------
static void GetBestF0Contour(int f0_length,
    const float* const* f0_candidates, const float* const* f0_scores,
    int number_of_bands, float* best_f0_contour)
{
    float tmp;
    for (int i = 0; i < f0_length; ++i) {
        tmp = f0_scores[0][i];
        best_f0_contour[i] = f0_candidates[0][i];
        for (int j = 1; j < number_of_bands; ++j) {
            if (tmp > f0_scores[j][i]) {
                tmp = f0_scores[j][i];
                best_f0_contour[i] = f0_candidates[j][i];
            }
        }
    }
}

//-----------------------------------------------------------------------------
// FixStep1() is the 1st step of the postprocessing.
// This function eliminates the unnatural change of f0 based on allowed_range.
//-----------------------------------------------------------------------------
static void FixStep1(const float* best_f0_contour, int f0_length,
    int voice_range_minimum, float allowed_range, float* f0_step1)
{
    float* f0_base = new float[f0_length];
    // Initialization
    for (int i = 0; i < voice_range_minimum; ++i)
        f0_base[i] = 0.0;
    for (int i = voice_range_minimum; i < f0_length - voice_range_minimum; ++i)
        f0_base[i] = best_f0_contour[i];
    for (int i = f0_length - voice_range_minimum; i < f0_length; ++i)
        f0_base[i] = 0.0;

    // Processing to prevent the jumping of f0
    for (int i = 0; i < voice_range_minimum; ++i)
        f0_step1[i] = 0.0f;
    for (int i = voice_range_minimum; i < f0_length; ++i)
        f0_step1[i] = fabsf((f0_base[i] - f0_base[i - 1]) / (kMySafeGuardMinimum + f0_base[i])) < allowed_range ? f0_base[i] : 0.0f;

    delete[] f0_base;
}

//-----------------------------------------------------------------------------
// FixStep2() is the 2nd step of the postprocessing.
// This function eliminates the suspected f0 in the anlaut and auslaut.
//-----------------------------------------------------------------------------
static void FixStep2(const float* f0_step1, int f0_length,
    int voice_range_minimum, float* f0_step2)
{
    for (int i = 0; i < f0_length; ++i)
        f0_step2[i] = f0_step1[i];

    int center = (voice_range_minimum - 1) / 2;
    for (int i = center; i < f0_length - center; ++i) {
        for (int j = -center; j <= center; ++j) {
            if (f0_step1[i + j] == 0) {
                f0_step2[i] = 0.0;
                break;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// GetNumberOfVoicedSections() counts the number of voiced sections.
//-----------------------------------------------------------------------------
static void GetNumberOfVoicedSections(const float* f0, int f0_length,
    int* positive_index, int* negative_index, int* positive_count,
    int* negative_count)
{
    *positive_count = *negative_count = 0;
    for (int i = 1; i < f0_length; ++i)
        if (f0[i] == 0 && f0[i - 1] != 0)
            negative_index[(*negative_count)++] = i - 1;
        else if (f0[i - 1] == 0 && f0[i] != 0)
            positive_index[(*positive_count)++] = i;
}

//-----------------------------------------------------------------------------
// SelectOneF0() corrects the f0[current_index] based on
// f0[current_index + sign].
//-----------------------------------------------------------------------------
static float SelectBestF0(float current_f0, float past_f0,
    const float* const* f0_candidates, int number_of_candidates,
    int target_index, float allowed_range)
{
    float reference_f0 = (current_f0 * 3.0f - past_f0) * 0.5f;

    float minimum_error = fabsf(reference_f0 - f0_candidates[0][target_index]);
    float best_f0 = f0_candidates[0][target_index];

    float current_error;
    for (int i = 1; i < number_of_candidates; ++i) {
        current_error = fabsf(reference_f0 - f0_candidates[i][target_index]);
        if (current_error < minimum_error) {
            minimum_error = current_error;
            best_f0 = f0_candidates[i][target_index];
        }
    }
    if (fabsf(1.0f - best_f0 / reference_f0) > allowed_range)
        return 0.0;
    return best_f0;
}

//-----------------------------------------------------------------------------
// FixStep3() is the 3rd step of the postprocessing.
// This function corrects the f0 candidates from backward to forward.
//-----------------------------------------------------------------------------
static void FixStep3(const float* f0_step2, int f0_length,
    const float* const* f0_candidates, int number_of_candidates,
    float allowed_range, const int* negative_index, int negative_count,
    float* f0_step3)
{
    for (int i = 0; i < f0_length; i++)
        f0_step3[i] = f0_step2[i];

    int limit;
    for (int i = 0; i < negative_count; ++i) {
        limit = i == negative_count - 1 ? f0_length - 1 : negative_index[i + 1];
        for (int j = negative_index[i]; j < limit; ++j) {
            f0_step3[j + 1] = SelectBestF0(f0_step3[j], f0_step3[j - 1], f0_candidates,
                number_of_candidates, j + 1, allowed_range);
            if (f0_step3[j + 1] == 0)
                break;
        }
    }
}

//-----------------------------------------------------------------------------
// FixStep4() is the 4th step of the postprocessing.
// This function corrects the f0 candidates from forward to backward.
//-----------------------------------------------------------------------------
static void FixStep4(const float* f0_step3, int f0_length,
    const float* const* f0_candidates, int number_of_candidates,
    float allowed_range, const int* positive_index, int positive_count,
    float* f0_step4)
{
    for (int i = 0; i < f0_length; ++i)
        f0_step4[i] = f0_step3[i];

    int limit;
    for (int i = positive_count - 1; i >= 0; --i) {
        limit = i == 0 ? 1 : positive_index[i - 1];
        for (int j = positive_index[i]; j > limit; --j) {
            f0_step4[j - 1] = SelectBestF0(f0_step4[j], f0_step4[j + 1], f0_candidates,
                number_of_candidates, j - 1, allowed_range);
            if (f0_step4[j - 1] == 0)
                break;
        }
    }
}

//-----------------------------------------------------------------------------
// FixF0Contour() calculates the definitive f0 contour based on all f0
// candidates. There are four steps.
//-----------------------------------------------------------------------------
static void FixF0Contour(float frame_period, int number_of_candidates,
    int fs, const float* const* f0_candidates,
    const float* best_f0_contour, int f0_length, float f0_floor,
    float allowed_range, float* fixed_f0_contour)
{
    int voice_range_minimum = static_cast<int>(0.5f + 1000.0f / frame_period / f0_floor) * 2 + 1;

    if (f0_length <= voice_range_minimum)
        return;

    float* f0_tmp1 = new float[f0_length];
    float* f0_tmp2 = new float[f0_length];

    FixStep1(best_f0_contour, f0_length, voice_range_minimum,
        allowed_range, f0_tmp1);
    FixStep2(f0_tmp1, f0_length, voice_range_minimum, f0_tmp2);

    int positive_count, negative_count;
    int* positive_index = new int[f0_length];
    int* negative_index = new int[f0_length];
    GetNumberOfVoicedSections(f0_tmp2, f0_length, positive_index,
        negative_index, &positive_count, &negative_count);
    FixStep3(f0_tmp2, f0_length, f0_candidates, number_of_candidates,
        allowed_range, negative_index, negative_count, f0_tmp1);
    FixStep4(f0_tmp1, f0_length, f0_candidates, number_of_candidates,
        allowed_range, positive_index, positive_count, fixed_f0_contour);

    delete[] f0_tmp1;
    delete[] f0_tmp2;
    delete[] positive_index;
    delete[] negative_index;
}

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in RawEventByDio()
//-----------------------------------------------------------------------------
static void GetFilteredSignal(int half_average_length, int fft_size,
    const cmplx* y_spectrum, int y_length, float* filtered_signal)
{
    float* low_pass_filter = new float[fft_size];
    // Nuttall window is used as a low-pass filter.
    // Cutoff frequency depends on the window length.
    NuttallWindow(half_average_length * 4, low_pass_filter);
    for (int i = half_average_length * 4; i < fft_size; ++i)
        low_pass_filter[i] = 0.0;

    cmplx* low_pass_filter_spectrum = new cmplx[fft_size];

    stb_fft_real_plan* forwardFFT = stb_fft_real_plan_dft_1d(fft_size);
    stb_fft_r2c_exec(forwardFFT, low_pass_filter, low_pass_filter_spectrum);

    // Convolution
    float tmp = y_spectrum[0].real * low_pass_filter_spectrum[0].real - y_spectrum[0].imag * low_pass_filter_spectrum[0].imag;
    low_pass_filter_spectrum[0].imag = y_spectrum[0].real * low_pass_filter_spectrum[0].imag + y_spectrum[0].imag * low_pass_filter_spectrum[0].real;
    low_pass_filter_spectrum[0].real = tmp;
    for (int i = 1; i <= fft_size / 2; ++i) {
        tmp = y_spectrum[i].real * low_pass_filter_spectrum[i].real - y_spectrum[i].imag * low_pass_filter_spectrum[i].imag;
        low_pass_filter_spectrum[i].imag = y_spectrum[i].real * low_pass_filter_spectrum[i].imag + y_spectrum[i].imag * low_pass_filter_spectrum[i].real;
        low_pass_filter_spectrum[i].real = tmp;
        low_pass_filter_spectrum[fft_size - i - 1].real = low_pass_filter_spectrum[i].real;
        low_pass_filter_spectrum[fft_size - i - 1].imag = low_pass_filter_spectrum[i].imag;
    }

    stb_fft_c2r_exec(forwardFFT, low_pass_filter_spectrum, filtered_signal);
    // Compensation of the delay.
    int index_bias = half_average_length * 2;
    for (int i = 0; i < y_length; ++i)
        filtered_signal[i] = filtered_signal[i + index_bias];

    free(forwardFFT);
    delete[] low_pass_filter_spectrum;
    delete[] low_pass_filter;
}

//-----------------------------------------------------------------------------
// CheckEvent() returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
//-----------------------------------------------------------------------------
static inline int CheckEvent(int x)
{
    return x > 0 ? 1 : 0;
}

//-----------------------------------------------------------------------------
// ZeroCrossingEngine() calculates the zero crossing points from positive to
// negative. Thanks to Custom.Maid http://custom-made.seesaa.net/ (2012/8/19)
//-----------------------------------------------------------------------------
static int ZeroCrossingEngine(const float* filtered_signal, int y_length,
    float fs, float* interval_locations, float* intervals)
{
    int* negative_going_points = new int[y_length];

    for (int i = 0; i < y_length - 1; ++i)
        negative_going_points[i] = 0.0 < filtered_signal[i] && filtered_signal[i + 1] <= 0.0 ? i + 1 : 0;
    negative_going_points[y_length - 1] = 0;

    int* edges = new int[y_length];
    int count = 0;
    for (int i = 0; i < y_length; ++i)
        if (negative_going_points[i] > 0)
            edges[count++] = negative_going_points[i];

    if (count < 2) {
        delete[] edges;
        delete[] negative_going_points;
        return 0;
    }

    float* fine_edges = new float[count];
    for (int i = 0; i < count; ++i)
        fine_edges[i] = edges[i] - filtered_signal[edges[i] - 1] / (filtered_signal[edges[i]] - filtered_signal[edges[i] - 1]);

    for (int i = 0; i < count - 1; ++i) {
        intervals[i] = fs / (fine_edges[i + 1] - fine_edges[i]);
        interval_locations[i] = (fine_edges[i] + fine_edges[i + 1]) * 0.5f / fs;
    }

    delete[] fine_edges;
    delete[] edges;
    delete[] negative_going_points;
    return count - 1;
}

//-----------------------------------------------------------------------------
// GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
//-----------------------------------------------------------------------------
static void GetFourZeroCrossingIntervals(float* filtered_signal, int y_length,
    float actual_fs, ZeroCrossings* zero_crossings)
{
    // x_length / 4 (old version) is fixed at 2013/07/14
    const int kMaximumNumber = y_length;
    zero_crossings->negative_interval_locations = new float[kMaximumNumber];
    zero_crossings->positive_interval_locations = new float[kMaximumNumber];
    zero_crossings->peak_interval_locations = new float[kMaximumNumber];
    zero_crossings->dip_interval_locations = new float[kMaximumNumber];
    zero_crossings->negative_intervals = new float[kMaximumNumber];
    zero_crossings->positive_intervals = new float[kMaximumNumber];
    zero_crossings->peak_intervals = new float[kMaximumNumber];
    zero_crossings->dip_intervals = new float[kMaximumNumber];

    zero_crossings->number_of_negatives = ZeroCrossingEngine(filtered_signal,
        y_length, actual_fs,
        zero_crossings->negative_interval_locations,
        zero_crossings->negative_intervals);

    for (int i = 0; i < y_length; ++i)
        filtered_signal[i] = -filtered_signal[i];
    zero_crossings->number_of_positives = ZeroCrossingEngine(filtered_signal,
        y_length, actual_fs,
        zero_crossings->positive_interval_locations,
        zero_crossings->positive_intervals);

    for (int i = 0; i < y_length - 1; ++i)
        filtered_signal[i] = filtered_signal[i] - filtered_signal[i + 1];
    zero_crossings->number_of_peaks = ZeroCrossingEngine(filtered_signal,
        y_length - 1, actual_fs,
        zero_crossings->peak_interval_locations,
        zero_crossings->peak_intervals);

    for (int i = 0; i < y_length - 1; ++i)
        filtered_signal[i] = -filtered_signal[i];
    zero_crossings->number_of_dips = ZeroCrossingEngine(filtered_signal,
        y_length - 1, actual_fs,
        zero_crossings->dip_interval_locations,
        zero_crossings->dip_intervals);
}

//-----------------------------------------------------------------------------
// GetF0CandidateContourSub() calculates the f0 candidates and deviations.
// This is the sub-function of GetF0Candidates() and assumes the calculation.
//-----------------------------------------------------------------------------
static void GetF0CandidateContourSub(
    const float* const* interpolated_f0_set, int f0_length, float f0_floor,
    float f0_ceil, float boundary_f0, float* f0_candidate,
    float* f0_score)
{
    for (int i = 0; i < f0_length; ++i) {
        f0_candidate[i] = (interpolated_f0_set[0][i] + interpolated_f0_set[1][i] + interpolated_f0_set[2][i] + interpolated_f0_set[3][i]) * 0.25f;

        f0_score[i] = sqrt(((interpolated_f0_set[0][i] - f0_candidate[i]) * (interpolated_f0_set[0][i] - f0_candidate[i]) + (interpolated_f0_set[1][i] - f0_candidate[i]) * (interpolated_f0_set[1][i] - f0_candidate[i]) + (interpolated_f0_set[2][i] - f0_candidate[i]) * (interpolated_f0_set[2][i] - f0_candidate[i]) + (interpolated_f0_set[3][i] - f0_candidate[i]) * (interpolated_f0_set[3][i] - f0_candidate[i])) * (1.0f / 3.0f));

        if (f0_candidate[i] > boundary_f0 || f0_candidate[i] < boundary_f0 * 0.5f || f0_candidate[i] > f0_ceil || f0_candidate[i] < f0_floor) {
            f0_candidate[i] = 0.0;
            f0_score[i] = kMaximumValue;
        }
    }
}

//-----------------------------------------------------------------------------
// GetF0CandidateContour() calculates the F0 candidates based on the
// zero-crossings.
//-----------------------------------------------------------------------------
static void GetF0CandidateContour(const ZeroCrossings* zero_crossings,
    float boundary_f0, float f0_floor, float f0_ceil,
    const float* temporal_positions, int f0_length,
    float* f0_candidate, float* f0_score)
{
    if (0 == CheckEvent(zero_crossings->number_of_negatives - 2) * CheckEvent(zero_crossings->number_of_positives - 2) * CheckEvent(zero_crossings->number_of_peaks - 2) * CheckEvent(zero_crossings->number_of_dips - 2)) {
        for (int i = 0; i < f0_length; ++i) {
            f0_score[i] = kMaximumValue;
            f0_candidate[i] = 0.0;
        }
        return;
    }

    float* interpolated_f0_set[4];
    for (int i = 0; i < 4; ++i)
        interpolated_f0_set[i] = new float[f0_length];

    interp1(zero_crossings->negative_interval_locations,
        zero_crossings->negative_intervals,
        zero_crossings->number_of_negatives,
        temporal_positions, f0_length, interpolated_f0_set[0]);
    interp1(zero_crossings->positive_interval_locations,
        zero_crossings->positive_intervals,
        zero_crossings->number_of_positives,
        temporal_positions, f0_length, interpolated_f0_set[1]);
    interp1(zero_crossings->peak_interval_locations,
        zero_crossings->peak_intervals, zero_crossings->number_of_peaks,
        temporal_positions, f0_length, interpolated_f0_set[2]);
    interp1(zero_crossings->dip_interval_locations,
        zero_crossings->dip_intervals, zero_crossings->number_of_dips,
        temporal_positions, f0_length, interpolated_f0_set[3]);

    GetF0CandidateContourSub(interpolated_f0_set, f0_length, f0_floor,
        f0_ceil, boundary_f0, f0_candidate, f0_score);
    for (int i = 0; i < 4; ++i)
        delete[] interpolated_f0_set[i];
}

//-----------------------------------------------------------------------------
// DestroyZeroCrossings() frees the memory of array in the struct
//-----------------------------------------------------------------------------
static void DestroyZeroCrossings(ZeroCrossings* zero_crossings)
{
    delete[] zero_crossings->negative_interval_locations;
    delete[] zero_crossings->positive_interval_locations;
    delete[] zero_crossings->peak_interval_locations;
    delete[] zero_crossings->dip_interval_locations;
    delete[] zero_crossings->negative_intervals;
    delete[] zero_crossings->positive_intervals;
    delete[] zero_crossings->peak_intervals;
    delete[] zero_crossings->dip_intervals;
}

//-----------------------------------------------------------------------------
// GetF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
//-----------------------------------------------------------------------------
static void GetF0CandidateFromRawEvent(float boundary_f0, float fs,
    const cmplx* y_spectrum, int y_length, int fft_size, float f0_floor,
    float f0_ceil, const float* temporal_positions, int f0_length,
    float* f0_score, float* f0_candidate)
{
    float* filtered_signal = new float[fft_size];
    GetFilteredSignal(matlab_round(fs / boundary_f0 * 0.5f), fft_size, y_spectrum,
        y_length, filtered_signal);

    ZeroCrossings zero_crossings = { 0 };
    GetFourZeroCrossingIntervals(filtered_signal, y_length, fs,
        &zero_crossings);

    GetF0CandidateContour(&zero_crossings, boundary_f0, f0_floor, f0_ceil,
        temporal_positions, f0_length, f0_candidate, f0_score);

    DestroyZeroCrossings(&zero_crossings);
    delete[] filtered_signal;
}

//-----------------------------------------------------------------------------
// GetF0CandidatesAndScores() calculates all f0 candidates and their scores.
//-----------------------------------------------------------------------------
static void GetF0CandidatesAndScores(const float* boundary_f0_list,
    int number_of_bands, float actual_fs, int y_length,
    const float* temporal_positions, int f0_length,
    const cmplx* y_spectrum, int fft_size, float f0_floor,
    float f0_ceil, float** raw_f0_candidates, float** raw_f0_scores)
{
    float* f0_candidate = new float[f0_length];
    float* f0_score = new float[f0_length];

    // Calculation of the acoustics events (zero-crossing)
    for (int i = 0; i < number_of_bands; ++i) {
        GetF0CandidateFromRawEvent(boundary_f0_list[i], actual_fs, y_spectrum,
            y_length, fft_size, f0_floor, f0_ceil, temporal_positions, f0_length,
            f0_score, f0_candidate);
        for (int j = 0; j < f0_length; ++j) {
            // A way to avoid zero division
            raw_f0_scores[i][j] = f0_score[j] / (f0_candidate[j] + kMySafeGuardMinimum);
            raw_f0_candidates[i][j] = f0_candidate[j];
        }
    }

    delete[] f0_candidate;
    delete[] f0_score;
}

//-----------------------------------------------------------------------------
// DioGeneralBody() estimates the F0 based on Distributed Inline-filter
// Operation.
//-----------------------------------------------------------------------------
static void DioGeneralBody(const float* x, int x_length, int fs,
    float frame_period, float f0_floor, float f0_ceil,
    float channels_in_octave, int speed, float allowed_range,
    float* temporal_positions, float* f0)
{
    int number_of_bands = 1 + static_cast<int>(log(f0_ceil / f0_floor) / kLog2 * channels_in_octave);
    float* boundary_f0_list = new float[number_of_bands];
    for (int i = 0; i < number_of_bands; ++i)
        boundary_f0_list[i] = f0_floor * powf(2.0f, (i + 1) / channels_in_octave);

    // normalization
    int decimation_ratio = MAX(MIN(speed, 12), 1);
    int y_length = (1 + static_cast<int>(x_length / decimation_ratio));
    float actual_fs = static_cast<float>(fs) / decimation_ratio;
    int fft_size = GetSuitableFFTSize(y_length + (4 * static_cast<int>(1.0f + actual_fs / boundary_f0_list[0] * 0.5f)));

    // Calculation of the spectrum used for the f0 estimation
    cmplx* y_spectrum = new cmplx[fft_size];
    GetSpectrumForEstimation(x, x_length, y_length, actual_fs, fft_size,
        decimation_ratio, y_spectrum);

    float** f0_candidates = new float*[number_of_bands];
    float** f0_scores = new float*[number_of_bands];
    int f0_length = GetSamplesForDIO(fs, x_length, frame_period);
    for (int i = 0; i < number_of_bands; ++i) {
        f0_candidates[i] = new float[f0_length];
        f0_scores[i] = new float[f0_length];
    }

    for (int i = 0; i < f0_length; ++i)
        temporal_positions[i] = i * frame_period / 1000.0f;

    GetF0CandidatesAndScores(boundary_f0_list, number_of_bands,
        actual_fs, y_length, temporal_positions, f0_length, y_spectrum,
        fft_size, f0_floor, f0_ceil, f0_candidates, f0_scores);

    // Selection of the best value based on fundamental-ness.
    // This function is related with SortCandidates() in MATLAB.
    float* best_f0_contour = new float[f0_length];
    GetBestF0Contour(f0_length, f0_candidates, f0_scores,
        number_of_bands, best_f0_contour);

    // Postprocessing to find the best f0-contour.
    FixF0Contour(frame_period, number_of_bands, fs, f0_candidates,
        best_f0_contour, f0_length, f0_floor, allowed_range, f0);

    delete[] best_f0_contour;
    delete[] y_spectrum;
    for (int i = 0; i < number_of_bands; ++i) {
        delete[] f0_scores[i];
        delete[] f0_candidates[i];
    }
    delete[] f0_scores;
    delete[] f0_candidates;
    delete[] boundary_f0_list;
}

//-----------------------------------------------------------------------------
// GetBaseIndex() calculates the temporal positions for windowing.
// Since the result includes negative value and the value that exceeds the
// length of the input signal, it must be modified appropriately.
//-----------------------------------------------------------------------------
static void GetBaseIndex(float current_position, const float* base_time,
    int base_time_length, int fs, int* index_raw)
{
    for (int i = 0; i < base_time_length; ++i)
        index_raw[i] = matlab_round((current_position + base_time[i]) * fs);
}

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
static void GetMainWindow(float current_position, const int* index_raw,
    int base_time_length, int fs, float window_length_in_time,
    float* main_window)
{
    float tmp = 0.0;
    for (int i = 0; i < base_time_length; ++i) {
        tmp = (index_raw[i] - 1.0f) / fs - current_position;
        main_window[i] = 0.42f + 0.5f * cosf(2.0f * kPi * tmp / window_length_in_time) + 0.08f * cosf(4.0f * kPi * tmp / window_length_in_time);
    }
}

//-----------------------------------------------------------------------------
// GetDiffWindow() generates the differentiated window.
// Diff means differential.
//-----------------------------------------------------------------------------
static void GetDiffWindow(const float* main_window, int base_time_length,
    float* diff_window)
{
    diff_window[0] = -main_window[1] * 0.5f;
    for (int i = 1; i < base_time_length - 1; ++i)
        diff_window[i] = -(main_window[i + 1] - main_window[i - 1]) * 0.5f;
    diff_window[base_time_length - 1] = main_window[base_time_length - 2] * 0.5f;
}

//-----------------------------------------------------------------------------
// GetSpectra() calculates two spectra of the waveform windowed by windows
// (main window and diff window).
//-----------------------------------------------------------------------------
static void GetSpectra(const float* x, int x_length, int fft_size,
    const int* index_raw, const float* main_window, const float* diff_window,
    int base_time_length, const ForwardRealFFT* forward_real_fft,
    cmplx* main_spectrum, cmplx* diff_spectrum)
{
    int* index = new int[base_time_length];

    for (int i = 0; i < base_time_length; ++i)
        index[i] = MAX(0, MIN(x_length - 1, index_raw[i] - 1));
    for (int i = 0; i < base_time_length; ++i)
        forward_real_fft->waveform[i] = x[index[i]] * main_window[i];
    for (int i = base_time_length; i < fft_size; ++i)
        forward_real_fft->waveform[i] = 0.0;

    stb_fft_r2c_exec(forward_real_fft->forward_fft, forward_real_fft->waveform, forward_real_fft->spectrum);
    for (int i = 0; i <= fft_size / 2; ++i) {
        main_spectrum[i].real = forward_real_fft->spectrum[i].real;
        main_spectrum[i].imag = forward_real_fft->spectrum[i].imag;
    }

    for (int i = 0; i < base_time_length; ++i)
        forward_real_fft->waveform[i] = x[index[i]] * diff_window[i];
    for (int i = base_time_length; i < fft_size; ++i)
        forward_real_fft->waveform[i] = 0.0;
    stb_fft_r2c_exec(forward_real_fft->forward_fft, forward_real_fft->waveform, forward_real_fft->spectrum);
    for (int i = 0; i <= fft_size / 2; ++i) {
        diff_spectrum[i].real = forward_real_fft->spectrum[i].real;
        diff_spectrum[i].imag = forward_real_fft->spectrum[i].imag;
    }

    delete[] index;
}

//-----------------------------------------------------------------------------
// FixF0() fixed the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
static float FixF0(const float* power_spectrum, const float* numerator_i,
    int fft_size, int fs, float initial_f0, int number_of_harmonics)
{
    float* amplitude_list = new float[number_of_harmonics];
    float* instantaneous_frequency_list = new float[number_of_harmonics];
    int index;
    for (int i = 0; i < number_of_harmonics; ++i) {
        index = matlab_round(initial_f0 * fft_size / fs * (i + 1));
        instantaneous_frequency_list[i] = power_spectrum[index] == 0.0f ? 0.0f : static_cast<float>(index) * fs / fft_size + numerator_i[index] / power_spectrum[index] * fs * 0.5f / kPi;
        amplitude_list[i] = sqrt(power_spectrum[index]);
    }
    float denominator = 0.0;
    float numerator = 0.0;
    for (int i = 0; i < number_of_harmonics; ++i) {
        numerator += amplitude_list[i] * instantaneous_frequency_list[i];
        denominator += amplitude_list[i] * (i + 1);
    }
    delete[] amplitude_list;
    delete[] instantaneous_frequency_list;
    return numerator / (denominator + kMySafeGuardMinimum);
}

//-----------------------------------------------------------------------------
// GetTentativeF0() calculates the F0 based on the instantaneous frequency.
//-----------------------------------------------------------------------------
static float GetTentativeF0(const float* power_spectrum,
    const float* numerator_i, int fft_size, int fs, float initial_f0)
{
    float tentative_f0 = FixF0(power_spectrum, numerator_i, fft_size, fs, initial_f0, 2);

    // If the fixed value is too large, the result will be rejected.
    if (tentative_f0 <= 0.0 || tentative_f0 > initial_f0 * 2)
        return 0.0;

    return FixF0(power_spectrum, numerator_i, fft_size, fs, tentative_f0, 6);
}

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------
static float GetMeanF0(const float* x, int x_length, int fs,
    float current_position, float initial_f0, int fft_size,
    float window_length_in_time, const float* base_time,
    int base_time_length)
{
    ForwardRealFFT forward_real_fft = { 0 };
    InitializeForwardRealFFT(fft_size, &forward_real_fft);
    cmplx* main_spectrum = new cmplx[fft_size];
    cmplx* diff_spectrum = new cmplx[fft_size];

    int* index_raw = new int[base_time_length];
    float* main_window = new float[base_time_length];
    float* diff_window = new float[base_time_length];

    GetBaseIndex(current_position, base_time, base_time_length, fs, index_raw);
    GetMainWindow(current_position, index_raw, base_time_length, fs,
        window_length_in_time, main_window);
    GetDiffWindow(main_window, base_time_length, diff_window);
    GetSpectra(x, x_length, fft_size, index_raw, main_window, diff_window,
        base_time_length, &forward_real_fft, main_spectrum, diff_spectrum);

    float* power_spectrum = new float[fft_size / 2 + 1];
    float* numerator_i = new float[fft_size / 2 + 1];
    for (int j = 0; j <= fft_size / 2; ++j) {
        numerator_i[j] = main_spectrum[j].real * diff_spectrum[j].imag - main_spectrum[j].imag * diff_spectrum[j].real;
        power_spectrum[j] = main_spectrum[j].real * main_spectrum[j].real + main_spectrum[j].imag * main_spectrum[j].imag;
    }

    float tentative_f0 = GetTentativeF0(power_spectrum, numerator_i,
        fft_size, fs, initial_f0);

    delete[] diff_spectrum;
    delete[] diff_window;
    delete[] main_window;
    delete[] index_raw;
    delete[] numerator_i;
    delete[] power_spectrum;
    delete[] main_spectrum;
    DestroyForwardRealFFT(&forward_real_fft);

    return tentative_f0;
}

//-----------------------------------------------------------------------------
// GetRefinedF0() fixes the F0 estimated by Dio(). This function uses
// instantaneous frequency.
//-----------------------------------------------------------------------------
static float GetRefinedF0(const float* x, int x_length, int fs,
    float current_potision, float initial_f0)
{
    if (initial_f0 <= kFloorF0StoneMask || initial_f0 > fs / 12.0)
        return 0.0;

    int half_window_length = static_cast<int>(1.5f * fs / initial_f0 + 1.0f);
    float window_length_in_time = (2.0f * half_window_length + 1.0f) / fs;
    float* base_time = new float[half_window_length * 2 + 1];
    for (int i = 0; i < half_window_length * 2 + 1; i++)
        base_time[i] = static_cast<float>(-half_window_length + i) / fs;
    int fft_size = static_cast<int>(powf(2.0f, 2.0f + static_cast<int>(log(half_window_length * 2.0f + 1.0f) / kLog2)));

    float mean_f0 = GetMeanF0(x, x_length, fs, current_potision,
        initial_f0, fft_size, window_length_in_time, base_time,
        half_window_length * 2 + 1);

    // If amount of correction is overlarge (20 %), initial F0 is employed.
    if (fabsf(mean_f0 - initial_f0) / initial_f0 > 0.2f)
        mean_f0 = initial_f0;

    delete[] base_time;

    return mean_f0;
}

} // namespace

int GetSamplesForDIO(int fs, int x_length, float frame_period)
{
    return static_cast<int>(1000.0 * x_length / fs / frame_period) + 1;
}

void Dio(const float* x, int x_length, int fs, const DioOption* option,
    float* temporal_positions, float* f0)
{
    DioGeneralBody(x, x_length, fs, option->frame_period, option->f0_floor,
        option->f0_ceil, option->channels_in_octave, option->speed,
        option->allowed_range, temporal_positions, f0);
}

void InitializeDioOption(DioOption* option)
{
    // You can change default parameters.
    option->channels_in_octave = 2.0f;
    option->f0_ceil = kCeilF0;
    option->f0_floor = kFloorF0;
    option->frame_period = 5;

    // You can use the value from 1 to 12.
    // Default value 11 is for the fs of 44.1 kHz.
    // The lower value you use, the better performance you can obtain.
    option->speed = 1;

    // You can give a positive real number as the threshold.
    // The most strict value is 0, and there is no upper limit.
    // On the other hand, I think that the value from 0.02 to 0.2 is reasonable.
    option->allowed_range = 0.1f;
}

void StoneMask(const float* x, int x_length, int fs,
    const float* temporal_positions, const float* f0, int f0_length,
    float* refined_f0)
{
    for (int i = 0; i < f0_length; i++)
        refined_f0[i] = GetRefinedF0(x, x_length, fs, temporal_positions[i], f0[i]);
}
