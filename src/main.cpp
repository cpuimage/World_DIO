#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if (defined(__WIN32__) || defined(_WIN32)) && !defined(__MINGW32__)
#include <conio.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
#pragma warning(disable : 4996)
#endif
#if (defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__))

#include <stdint.h>
#include <sys/time.h>

#endif

// For .wav input/output functions.
#include "audioio.h"

// WORLD core functions.
// Note: win.sln uses an option in Additional Include Directories.
// To compile the program, the option "-I $(SolutionDir)..\src" was set.
#include "dio.h"
#include "matlabfunctions.h"

#if (defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
// Linux porting section: implement timeGetTime() by gettimeofday(),
#ifndef DWORD
#define DWORD uint32_t
#endif

DWORD timeGetTime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    DWORD ret = static_cast<DWORD>(tv.tv_usec / 1000 + tv.tv_sec * 1000);
    return ret;
}

#endif

namespace {
void DisplayInformation(int fs, int input_length)
{
    printf("File information\n");
    printf("Sampling : %d Hz\n", fs);
    printf("Length %d [sample]\n", input_length);
    printf("Length %f [sec]\n", static_cast<float>(input_length) / fs);
}

//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
    float frame_period;
    int fs;

    float* f0;
    float* time_axis;
    int f0_length;

    int fft_size;
} DIO_Parameters;
#define DEBUG 1

void F0EstimationDio(float* input, int input_length,
    DIO_Parameters* world_parameters)
{
    DioOption option = { 0 };
    InitializeDioOption(&option);

    // Modification of the option
    option.frame_period = world_parameters->frame_period;

    // Valuable option.speed represents the ratio for downsampling.
    // The signal is downsampled to fs / speed Hz.
    // If you want to obtain the accurate result, speed should be set to 1.
    option.speed = 1;

    // You can set the f0_floor below world::kFloorF0.
    option.f0_floor = 40.0;

    // You can give a positive real number as the threshold.
    // Most strict value is 0, but almost all results are counted as unvoiced.
    // The value from 0.02 to 0.2 would be reasonable.
    option.allowed_range = 0.1;

    // Parameters setting and memory allocation.
    world_parameters->f0_length = GetSamplesForDIO(world_parameters->fs,
        input_length, world_parameters->frame_period);
    world_parameters->f0 = new float[world_parameters->f0_length];
    world_parameters->time_axis = new float[world_parameters->f0_length];
    float* refined_f0 = new float[world_parameters->f0_length];

    printf("\nAnalysis\n");
    DWORD elapsed_time = timeGetTime();
    Dio(input, input_length, world_parameters->fs, &option, world_parameters->time_axis,
        world_parameters->f0);
    printf("DIO: %d [msec]\n", timeGetTime() - elapsed_time);

    // StoneMask is carried out to improve the estimation performance.
    elapsed_time = timeGetTime();
    StoneMask(input, input_length, world_parameters->fs, world_parameters->time_axis,
        world_parameters->f0, world_parameters->f0_length, refined_f0);
    printf("StoneMask: %d [msec]\n", timeGetTime() - elapsed_time);

    for (int i = 0; i < world_parameters->f0_length; ++i)
        world_parameters->f0[i] = refined_f0[i];
#if DEBUG
    printf("f0 length %d: \n", world_parameters->f0_length);
    for (int i = 0; i < world_parameters->f0_length; ++i) {
        printf("f0 %f: time_axis %f: \n", world_parameters->f0[i], world_parameters->time_axis[i]);
    }
#endif
    delete[] refined_f0;
}

void DestroyMemory(DIO_Parameters* world_parameters)
{
    delete[] world_parameters->time_axis;
    delete[] world_parameters->f0;
}
} // namespace

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc < 2) {
        printf("usage : %s input.wav \n", argv[0]);
        return -2;
    }

    const char* inFile = argv[1];
    int input_length = 0;

    int fs;
    float* input = wavread(inFile, &fs, &input_length);
    if (input == NULL) {
        printf("error: The file is not .wav format.\n");
        return -1;
    }
    DisplayInformation(fs, input_length);

    //---------------------------------------------------------------------------
    // Analysis part
    //---------------------------------------------------------------------------
    DIO_Parameters parameters = { 0 };
    // You must set fs and frame_period before analysis/synthesis.
    parameters.fs = fs;
    // 5.0 ms is the default value.
    parameters.frame_period = 5.0;

    // F0 estimation
    // DIO
    F0EstimationDio(input, input_length, &parameters);

    free(input);
    DestroyMemory(&parameters);

    printf("complete.\n");
    return 0;
}