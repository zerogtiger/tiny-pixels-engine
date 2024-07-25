# Graphics Processor
A C++ based image and graphics processing library. 

## Features
- [ ] Fractional scaling (via `Image::f_scale()`)
    - [ ] Bicubic
    - [ ] Fourier-transform
    - [ ] Edge-directed interpolation
    - [ ] Lanczos
    - [ ] hqx
    - [ ] Mipmap
- [ ] Framing
- [ ] Masking
- [ ] Procedural texture generation

## Enhancements
- [ ] Edge bleeding in edge detection via `Image::edge()` function due to Fast Fourier transform limitations
- [ ] Sobel edge detection gradient for single channel images
- [ ] Make readme build/run more efficient
- [x] Refactor `src/` and move library files into `src/lib/`

