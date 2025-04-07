//   meanfilter.h - declarations for
//   1D and 2D mean filter routines
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#ifndef _MEANFILTER_H_
#define _MEANFILTER_H_

//   Signal/image element type
typedef double element;

//   1D MEAN FILTER, window size 5
//     signal - input signal
//     result - output signal, NULL for inplace processing
//     N      - length of the signal
void meanfilter(element *signal, element *result, int N);

//   2D MEAN FILTER, window size 3x3
//     image  - input image
//     result - output image, NULL for inplace processing
//     N      - width of the image
//     M      - height of the image
void meanfilter(element *image, element *result, int N, int M);

#endif
