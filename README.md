# Graphics Processor
A C++ based image and graphics processing library implemented from scratch.

## To-do
- Color
    - [x] Invert
    - [x] Histogram
    - [ ] Hue correct
    - [ ] Tonal correction
    - [x] Color balance (lift, gamma, gain)
    - [ ] Exposure
    - [x] Gamma
    - [x] HSV adjustments
    - [ ] RGB curves
    - [ ] Tone map
    - [ ] White balance
    - [ ] False color
    - [ ] Mixing options
        - [ ] Alpha over
        - [x] Separate colors (RGBA channels)
        - [x] Combine colors (RGBA channels)
        - [ ] Color ramp
            - [x] Constant
            - [x] Linear
            - [ ] Bézier
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

## High Priority Enhancements/Bug Fixes
- [x] Bilinear scaling (`Image::preview_color_ramp` then scale to a large image via `Image::f_scale` with `TwoDimInterp::Biliea` as interpolation method

## Medium Priority Enhancements/Bug Fixes
- [ ] Support separated channel images in `Image::HSV`
- [ ] Color class HSV integration
- [ ] Demo images in `README.md`
- [ ] Ensure uniformity of logic when using data from another image (first channel vs. average grayscale vs. separated channel image)
- [ ] Unnecessary memory use/leaks
- [x] Rename `enum` to generalize for 1-D and 2-D interpolation methods

## Low Priority Enhancements/Bug Fixes
- [ ] Clean up code in `Image::false_color` function
- [ ] Use `uint_t` and `size_t` where necessary
- [ ] Edge bleeding in edge detection via `Image::edge()` function due to Fast Fourier transform limitations
- [ ] Sobel edge detection gradient for single channel images
- [ ] Make makefile build/run more efficient
- [x] Refactor `src/` and move library files into `src/lib/`

