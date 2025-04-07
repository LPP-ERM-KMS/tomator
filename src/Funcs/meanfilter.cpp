//   meanfilter.cpp - impelementation of
//   1D and 2D mean filter routines
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#include "meanfilter.h"
#include <memory.h>

//   1D MEAN FILTER implementation
//     signal - input signal
//     result - output signal
//     N      - length of the signal
void _meanfilter(const element *signal, element *result, int N) {
    //   Move window through all elements of the signal
    for (int i = 2; i < N - 2; ++i)
        //   Take the average
        result[i - 2] = (signal[i - 2] +
                         signal[i - 1] +
                         signal[i] +
                         signal[i + 1] +
                         signal[i + 2]) /
                        5;
}

//   1D MEAN FILTER wrapper
//     signal - input signal
//     result - output signal
//     N      - length of the signal
void meanfilter(element *signal, element *result, int N) {
    //   Check arguments
    if (!signal || N < 1)
        return;
    //   Treat special case N = 1
    if (N == 1) {
        if (result)
            result[0] = signal[0];
        return;
    }
    //   Allocate memory for signal extension
    element *extension = new element[N + 4];
    //   Check memory allocation
    if (!extension)
        return;
    //   Create signal extension
    memcpy(extension + 2, signal, N * sizeof(element));
    for (int i = 0; i < 2; ++i) {
        extension[i] = signal[1 - i];
        extension[N + 2 + i] = signal[N - 1 - i];
    }
    //   Call mean filter implementation
    _meanfilter(extension, result ? result : signal, N + 4);
    //   Free memory
    delete[] extension;
}

//   2D MEAN FILTER implementation
//     image  - input image
//     result - output image
//     N      - width of the image
//     M      - height of the image
void _meanfilter(const element *image, element *result, int N, int M) {
    //   Move window through all elements of the image
    for (int m = 1; m < M - 1; ++m)
        for (int n = 1; n < N - 1; ++n)
            //   Take the average
            result[(m - 1) * (N - 2) + n - 1] = (image[(m - 1) * N + n - 1] +
                                                 image[(m - 1) * N + n] +
                                                 image[(m - 1) * N + n + 1] +
                                                 image[m * N + n - 1] +
                                                 image[m * N + n] +
                                                 image[m * N + n + 1] +
                                                 image[(m + 1) * N + n - 1] +
                                                 image[(m + 1) * N + n] +
                                                 image[(m + 1) * N + n + 1]) /
                                                9;
}

//   2D MEAN FILTER wrapper
//     image  - input image
//     result - output image
//     N      - width of the image
//     M      - height of the image
void meanfilter(element *image, element *result, int N, int M) {
    //   Check arguments
    if (!image || N < 1 || M < 1)
        return;
    //   Allocate memory for signal extension
    element *extension = new element[(N + 2) * (M + 2)];
    //   Check memory allocation
    if (!extension)
        return;
    //   Create image extension
    for (int i = 0; i < M; ++i) {
        memcpy(extension + (N + 2) * (i + 1) + 1,
               image + N * i,
               N * sizeof(element));
        extension[(N + 2) * (i + 1)] = image[N * i];
        extension[(N + 2) * (i + 2) - 1] = image[N * (i + 1) - 1];
    }
    //   Fill first line of image extension
    memcpy(extension, extension + N + 2, (N + 2) * sizeof(element));
    //   Fill last line of image extension
    memcpy(extension + (N + 2) * (M + 1), extension + (N + 2) * M, (N + 2) * sizeof(element));
    //   Call mean filter implementation
    _meanfilter(extension, result ? result : image, N + 2, M + 2);
    //   Free memory
    delete[] extension;
}
