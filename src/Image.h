#include "schrift.h"
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <stdint.h>

#define _USE_MATH_DEFINES

enum ImageType { PNG, JPG, BMP, TGA };

struct Font;

struct Image {
    uint8_t* data = NULL; // Butes of image
    size_t size = 0;      // size of data
    int w;
    int h;
    int channels;

    Image(const char* filename);
    Image(int w, int h, int channels);
    Image(const Image& img);
    ~Image(); // destructor

    bool read(const char* filename);
    bool write(const char* filename);

    ImageType getFileType(const char* filename);

    Image& grayscale_avg();
    Image& grayscale_lum();

    Image& color_mask(float r, float g, float b);

    Image& encodemessage(const char* message);

    Image& decodemessage(char* buffer, size_t* messageLength);

    Image& diffmap(Image& img);
    Image& diffmap_scale(Image& img, uint8_t scl = 0);

    Image& std_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& std_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& std_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    Image& flip_x();
    Image& flip_y();

    Image& overlay(const Image& src, int x, int y);

    Image& overlay_text(const char* txt, const Font& font, int x, int y, uint8_t r = 255, uint8_t g = 255, uint8_t b = 255, uint8_t a = 255);

    Image& crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch);
    static uint32_t rev(uint32_t n, uint32_t a); // n: 2^a
    static void bit_rev(uint32_t n, std::complex<double> a[], std::complex<double>* A);

    static void fft(uint32_t n, std::complex<double> x[], std::complex<double>* X);
    static void ifft(uint32_t n, std::complex<double> x[], std::complex<double>* X);
    static void dft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double>* X);
    static void idft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double>* X);

    static void pad_kernel(uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc, uint32_t pw, uint32_t ph, std::complex<double>* pad_ker);
    static inline void pointwise_product(uint64_t l, std::complex<double> a[], std::complex<double> b[], std::complex<double>* p);

    Image& fd_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& fd_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& fd_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    Image& convolve_linear(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    Image& brightness(uint8_t channel, double brightness_delta);
    Image& contrast(uint8_t channel, double contrast_delta);
};

struct Font {
    SFT sft = {NULL, 12, 12, 0, 0, SFT_DOWNWARD_Y | SFT_RENDER_IMAGE};
    Font(const char* fontfile, uint16_t size) {
        if ((sft.font = sft_loadfile(fontfile)) == NULL) {
            printf("\e[31m[ERROR] Failed to load %s\e[0m\n", fontfile);
            return;
        }
        setSize(size);
    }

    ~Font() { sft_freefont(sft.font); }

    void setSize(uint16_t size) {
        sft.xScale = size;
        sft.yScale = size;
    }
};
