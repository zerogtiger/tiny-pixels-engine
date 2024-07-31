#include "lib/schrift.h"
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <stdint.h>
#include <utility>

#define _USE_MATH_DEFINES

enum ImageType { PNG, JPG, BMP, TGA };
enum TwoDimInterp { Nearest, Bilinear, Bicubic, Fourier, Edge, HQX, Mipmap, RSamp, RShear, RAM };
enum OneDimInterp { Constant, Linear, BSpline };

struct Adjustment {
    double brightness = 0, contrast = 0, hue = 0, saturation = 0, value = 0, lift = 1, gamma = 1, gain = 1;

    Adjustment(double brightness, double contrast, double hue, double saturation, double value, double lift, double gamma, double gain)
        : brightness(brightness), contrast(contrast), hue(hue), saturation(saturation), value(value), lift(lift), gamma(gamma), gain(gain) {}

    ~Adjustment() {}

    Adjustment& create_adj_bcs(double brightness, double contrast, double hue, double saturation, double value) {
        Adjustment* ret = new Adjustment(brightness, contrast, hue, saturation, value, 1, 1, 1);
        return *ret;
    }

    Adjustment& create_adj_lgg(double lift, double gamma, double gain) {
        Adjustment* ret = new Adjustment(0, 0, 0, 0, 0, lift, gamma, gain);
        return *ret;
    }
};

struct Color {
    double r, g, b, a = 255;
    Color(double r, double g, double b) : r(r), g(g), b(b) {}
    Color(double r, double g, double b, double a) : r(r), g(g), b(b), a(a) {}
    ~Color() {}

    void hsv_to_rgb(double h_deg, double s, double v) {
        if (abs(s) > 1 || abs(v) > 1) {
            throw std::invalid_argument("saturation or value not in the range [0, 1]");
        }
        const int mod = 360;
        h_deg = fmod(fmod(h_deg, mod) + mod, mod);

        // hsl -> rgb
        double c = s * v;
        double x = c * (1 - abs(fmod((h_deg / 60), 2) - 1));
        double m = v - c;

        double rt, gt, bt;
        rt = gt = bt = 0;
        if (h_deg < 60) {
            rt = c;
            gt = x;
        } else if (h_deg < 120) {
            rt = x;
            gt = c;
        } else if (h_deg < 180) {
            gt = c;
            bt = x;
        } else if (h_deg < 240) {
            gt = x;
            bt = c;
        } else if (h_deg < 300) {
            bt = c;
            rt = x;
        } else {
            bt = x;
            rt = c;
        }

        r = std::clamp(round(255 * (rt + m)), 0.0, 255.0);
        g = std::clamp(round(255 * (gt + m)), 0.0, 255.0);
        b = std::clamp(round(255 * (bt + m)), 0.0, 255.0);
    }
    Color& rgb_to_hsv(double rr, double gg, double bb) {
        rr = rr / 255.0;
        gg = gg / 255.0;
        bb = bb / 255.0;
        double c_max = fmax(rr, fmax(gg, bb)), c_min = fmin(rr, fmin(gg, bb)), delta = c_max - c_min;
        Color* ret = new Color(0, 0, 0);
        if (delta == 0) {
            ret->r = 0;
        } else if (c_max == rr) {
            ret->r = 60.0 * fmod(fmod(((gg - bb) / delta), 6) + 6, 6);
        } else if (c_max == gg) {
            ret->r = 60.0 * (fmod(fmod((bb - rr) / delta, 6) + 6, 6) + 2);
        } else if (c_max == bb) {
            ret->r = 60.0 * (fmod(fmod((rr - gg) / delta, 6) + 6, 6) + 4);
        }

        if (c_max == 0) {
            ret->g = 0;
        } else {
            ret->g = delta / c_max;
        }
        ret->b = c_max;
        return *ret;
    }

    double get(int col) {
        switch (col) {
        case 0:
            return r;
            break;
        case 1:
            return g;
            break;
        case 2:
            return b;
            break;
        case 3:
            return a;
            break;
        default:
            return -1;
            break;
        }
    }
    bool set(int col, double val) {
        switch (col) {
        case 0:
            r = val;
            return true;
            break;
        case 1:
            g = val;
            return true;
            break;
        case 2:
            b = val;
            return true;
            break;
        case 3:
            a = val;
            return true;
            break;
        default:
            return false;
            break;
        }
    }

    bool operator<(const Color& other) const {
        return r == other.r ? (g == other.g ? b < other.b : g < other.g) : r < other.r;
    }

    Color operator+(const Color& other) const { return Color(r + other.r, g + other.g, b + other.b, a + other.a); }

    double luminance(double r, double g, double b) { return 0.2126 * r + 0.7152 * g + 0.0722 * b; }

    double luminance() { return 0.2126 * r + 0.7152 * g + 0.0722 * b; }

    static double luminance(const Color& c) { return 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b; }

    Color& apply_adj_rgb(Adjustment adj) {
        // double brightness = 0, contrast = 0, hue = 0, saturation = 0, value = 0, lift = 1, gamma = 1, gain = 1;

        double F = 259.0 * (adj.contrast + 255) / (255.0 * (259 - adj.contrast));
        r = std::clamp(F * (r - 128) + 128 + adj.brightness, 0.0, 255.0);
        g = std::clamp(F * (g - 128) + 128 + adj.brightness, 0.0, 255.0);
        b = std::clamp(F * (b - 128) + 128 + adj.brightness, 0.0, 255.0);

        *this = rgb_to_hsv(r, g, b);
        hsv_to_rgb(r + adj.hue, std::clamp(g + adj.saturation, 0.0, 1.1), std::clamp(b + adj.value, 0.0, 1.0));


        r = pow(adj.gain * (r/255.0 + adj.lift * (1 - r/255.0)), 1.0 / adj.gamma) * 255.0;
        g = pow(adj.gain * (g/255.0 + adj.lift * (1 - g/255.0)), 1.0 / adj.gamma) * 255.0;
        b = pow(adj.gain * (b/255.0 + adj.lift * (1 - b/255.0)), 1.0 / adj.gamma) * 255.0;

    }
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

struct Image {
    uint8_t* data = NULL; // Bytes of image
    size_t size = 0;      // size of data
    int w;
    int h;
    int channels;

    Image(const char* filename);
    Image(int w, int h, int channels);
    Image(const Image& img);
    // Image(const Image& img, int channel);
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

    Image& std_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                      uint32_t cc, bool normalize = false);
    Image& std_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                        uint32_t cc, bool normalize = false);
    Image& std_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc,
                               bool normalize = false);

    // void std_convolve_clamp_to_zero_double(double* new_data, uint8_t channel, uint32_t ker_w, uint32_t ker_h,
    // double ker[], uint32_t cr, uint32_t cc);
    Image& flip_x();
    Image& flip_y();

    Image& overlay(const Image& src, int x, int y);

    Image& overlay_text(const char* txt, const Font& font, int x, int y, uint8_t r = 255, uint8_t g = 255,
                        uint8_t b = 255, uint8_t a = 255);

    Image& crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch);
    static uint32_t rev(uint32_t n, uint32_t a); // n: 2^a
    static void bit_rev(uint32_t n, std::complex<double> a[], std::complex<double>* A);

    static void fft(uint32_t n, std::complex<double> x[], std::complex<double>* X);
    static void ifft(uint32_t n, std::complex<double> x[], std::complex<double>* X);
    static void dft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double>* X);
    static void idft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double>* X);

    static void pad_kernel(uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc, uint32_t pw,
                           uint32_t ph, std::complex<double>* pad_ker);
    static inline void pointwise_product(uint64_t l, std::complex<double> a[], std::complex<double> b[],
                                         std::complex<double>* p);

    std::complex<double>* fd_convolve_clamp_to_zero_raw(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[],
                                                        uint32_t cr, uint32_t cc);
    Image& fd_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                     uint32_t cc, bool normalize = false);
    std::complex<double>* fd_convolve_clamp_to_border_raw(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[],
                                                          uint32_t cr, uint32_t cc);
    Image& fd_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                       uint32_t cc, bool normalize = false);
    std::complex<double>* fd_convolve_cyclic_raw(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[],
                                                 uint32_t cr, uint32_t cc);
    Image& fd_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc,
                              bool normalize = false);

    Image& convolve_linear(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc,
                           bool normalize = false);
    Image& convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                    uint32_t cc, bool normalize = false);
    Image& convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc,
                           bool normalize = false);

    Image& brightness(uint8_t channel, double brightness_delta);
    Image& contrast(uint8_t channel, double contrast_delta);
    Image& saturation(uint8_t channel, double contrast_delta);

    Image& shade_h();
    Image& shade_v();
    Image& shade();
    Image& edge(bool gradient = false, double detail_threshold = 0.09);

    Image& f_scale(uint32_t new_w, uint32_t new_h, bool linked = false, TwoDimInterp method = TwoDimInterp::Nearest);
    Image& translate(int x, int y, Color fill = Color(0, 0, 0, 255));
    Image& rotate(double degrees); // unimplemented

    Image& invert_color(uint8_t channel);
    Image& gamma(uint8_t channel, double gamma_delta);

    Image& color_reduce(bool error_diffusion = true);

    Image& color_ramp(std::vector<std::pair<double, Color>> points, OneDimInterp method = OneDimInterp::Linear);
    Image& preview_color_ramp(std::vector<std::pair<double, Color>> points,
                              OneDimInterp method = OneDimInterp::Linear) const;

    std::vector<Image*> seperate_channels();
    Image& combine_channels(std::vector<Image*> imgs, bool resize_to_fit = false,
                            TwoDimInterp method = TwoDimInterp::Bilinear);
    Image& set_alpha(Image& alph, bool resize_to_fit = false, TwoDimInterp method = TwoDimInterp::Bilinear);

    Image& color_balance(Color lift, Color gamma, Color gain);

    Image& histogram(bool inc_lum = true, int channel = -1);
    Image& histogram_lum();

    Image& HSV(double hue_delta, double saturation_delta, double value_delta);

    Image& false_color(bool overwrite = false);

    Image& tone_correct(uint8_t midtones_start, uint8_t midtones_end, Adjustment shadow, Adjustment midtone, Adjustment highlight);
};
