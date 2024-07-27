# Graphics Processor
A C++ based image and graphics processing library implemented from scratch.

## Todo
- Color
    - [x] Invert
    - [ ] Color correction/balance
    - [ ] Exposure
    - [x] Gamma
    - [ ] HSV adjustments
    - [ ] RGB curves
    - [ ] Tone map
    - [ ] White balance
    - [ ] Mixing options
        - [ ] Alpha over
        - [x] Separate colors (RGBA channels)
        - [ ] Combine colors (RGBA channels)
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
    - [ ] Translate
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

## Medium Priority Bug Fixes
- [ ] Unnecessary memory use/leaks
- [ ] Rename `enum` to generalize for 1-D and 2-D interpolation methods

## Low Priority Bug Fixes
- [ ] Edge bleeding in edge detection via `Image::edge()` function due to Fast Fourier transform limitations
- [ ] Sobel edge detection gradient for single channel images
- [ ] Make makefile build/run more efficient
- [x] Refactor `src/` and move library files into `src/lib/`

