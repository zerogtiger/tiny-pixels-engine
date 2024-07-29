# Graphics Processor
A C++ based image and graphics processing library implemented from scratch.

## Todo
- Color
    - [x] Invert
    - [x] Histogram
    - [ ] Color correction
    - [x] Color balance (lift, gamma, gain)
    - [ ] Exposure
    - [x] Gamma
    - [ ] HSV adjustments
    - [ ] RGB curves
    - [ ] Tone map
    - [ ] White balance
    - [ ] Mixing options
        - [ ] Alpha over
        - [x] Separate colors (RGBA channels)
        - [x] Combine colors (RGBA channels)
        - [ ] Color ramp
            - [x] Constant
            - [x] Linear
            - [ ] Bezier
            - [ ] B-Spline
- Filtering
    - [ ] Color reduce
        - [x] 3-bit color
        - [x] 3-bit Floyd-Steinberg error diffusion
        - [ ] 8-bit color
        - [ ] 8-bit Floyd-Steinberg error diffusion
        - [ ] 8/16-bit color quantization
    - [ ] Blur
        - [ ] Box
        - [ ] Bokeh
        - [ ] Gaussian (3x3 case implemented)
    - [ ] Anti-aliasing
    - [ ] Masking
- Transform
    - [ ] Rotate
    - [ ] Scale
        - [x] Nearest
        - [x] Bilinear
        - [ ] Bicubic
        - [ ] Fourier
        - [ ] Edge
        - [ ] HQX
        - [ ] Mipmap
    - [x] Translate
    - [x] Crop
    - [x] Flip
    - [ ] Lens distortion
- Misc.
    - [ ] Frame
    - [x] Text overlay
    - [ ] Procedural texture generation
        - [ ] Brick
        - [ ] Gradient
        - [ ] Musgrave
        - [ ] Noise
        - [ ] Voronoi
        - [ ] Wave
- Documentation (not started)

## High Priority Enhancement/Bug Fixes
- [ ] Bilinear scaling (`Image::preview_color_ramp` then scale to a large image via `Image::f_scale` with `TwoDimInterp::Biliea` as interpolation method

## Medium Priority Enhancement/Bug Fixes
- [ ] Ensure uniformity of logic when using data from another image (first channel vs. average grayscale)
- [ ] Unnecessary memory use/leaks
- [x] Rename `enum` to generalize for 1-D and 2-D interpolation methods

## Low Priority Bug Fixes
- [ ] Use `uint_t` and `size_t` where necessary
- [ ] Edge bleeding in edge detection via `Image::edge()` function due to Fast Fourier transform limitations
- [ ] Sobel edge detection gradient for single channel images
- [ ] Make makefile build/run more efficient
- [x] Refactor `src/` and move library files into `src/lib/`

