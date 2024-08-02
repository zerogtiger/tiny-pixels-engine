#include "Enums.h"
#include "Interpolation.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <utility>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#define BYTE_BOUND(value) ((value) > 255 ? 255 : ((value) < 0 ? 0 : (value)))
#define SGN(value) ((value > 0) - (value < 0))
#define STEG_HEADER_SIZE sizeof(uint32_t) * 8

#include "Image.h"
#include "lib/schrift.h"
#include "lib/stb_image.h"
#include "lib/stb_image_write.h"
#include "string.h"

Image::Image(const char* filename) {
    if (read(filename)) {
        printf("Read %s\n", filename);
        size = w * h * channels;
    } else {
        printf("Failed to read %s\n", filename);
    }
}
Image::Image(int w, int h, int channels)
    : w(w), h(h), channels(channels) // initializer list
{
    size = w * h * channels;
    data = new uint8_t[size];
    memset(data, 0, size);
    // all black "image"
}
Image::Image(const Image& img) : Image(img.w, img.h, img.channels) { memcpy(data, img.data, img.size); }
Image::~Image() { stbi_image_free(data); }

bool Image::read(const char* filename) {
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
        data = stbi_load(filename, &w, &h, &channels, 0);
        return data != NULL;
    } else {
        throw std::invalid_argument(std::string{"Unable to find file "} + filename);
        return false;
    }
}
bool Image::write(const char* filename) {
    ImageType type = getFileType(filename);
    int success;
    switch (type) {
    case ImageType::PNG:
        success = stbi_write_png(filename, w, h, channels, data, w * channels);
        break;
    case ImageType::BMP:
        success = stbi_write_bmp(filename, w, h, channels, data);
        break;
    case ImageType::JPG:
        success = stbi_write_jpg(filename, w, h, channels, data, 100);
        break;
    case ImageType::TGA:
        success = stbi_write_tga(filename, w, h, channels, data);
        break;
    }
    return success != 0;
}

ImageType Image::getFileType(const char* filename) {
    const char* ext = strrchr(filename, '.');
    if (ext != nullptr) {
        if (strcmp(ext, ".png") == 0) {
            return ImageType::PNG;
        } else if (strcmp(ext, ".jpg") == 0) {
            return ImageType::JPG;
        } else if (strcmp(ext, ".BMP") == 0) {
            return ImageType::BMP;
        } else if (strcmp(ext, ".tga") == 0) {
            return ImageType::TGA;
        }
    }
    return PNG;
}

uint8_t Image::get(uint32_t row, uint32_t col, uint32_t channel) { return data[(row * w + col) * channels + channel]; }
uint8_t Image::get_or_default(uint32_t row, uint32_t col, uint32_t channel, uint8_t fallback) {
    if (row < 0 || row >= w || col < 0 || col >= h) {
        return fallback;
    } else {
        return get(row, col, channel);
    }
}
Color Image::get_color(uint32_t row, uint32_t col) {
    Color ret(0, 0, 0);
    for (int i = 0; i < 3; i++) {
        ret.set(i, get(row, col, (channels < 3 ? 0 : i)));
    }
    if (channels >= 4) {
        ret.a = get(row, col, 3);
    }
    return ret;
}
Color Image::get_color_or_default(int row, int col, Color fallback) {
    if (row < 0 || row >= h || col < 0 || col >= w) {
        return fallback;
    } else {
        return get_color(row, col);
    }
}

Image& Image::grayscale_avg() {
    if (channels < 3) {
        printf("Image %p has less than 3 channels, assumed to already be "
               "grayscale",
               this);
    } else {
        for (int i = 0; i < size; i += channels) {
            int gray = (data[i] + data[i + 1] + data[i + 2]) / 3;
            memset(data + i, gray, 3);
        }
    }
    return *this;
}
Image& Image::grayscale_lum() {
    if (channels < 3) {
        printf("Image %p has less than 3 channels, assumed to already be "
               "grayscale",
               this);
    } else {
        for (int i = 0; i < size; i += channels) {
            int gray = 0.2126 * data[i] + 0.7152 * data[i + 1] + 0.0722 * data[i + 2];
            memset(data + i, gray, 3);
        }
    }
    return *this;
}
Image& Image::color_mask(float r, float g, float b) {
    if (channels < 3) {
        printf("\e[31m[ERROR] Color mask requries at least 3 channels, but "
               "this image has %d channels \e[0m\n",
               channels);
    } else {
        for (int i = 0; i < size; i += channels) {
            data[i] *= r;
            data[i + 1] *= g;
            data[i + 2] *= b;
        }
    }
    return *this;
}

Image& Image::encodemessage(const char* message) {
    uint32_t len = strlen(message) * 8;

    if (len + STEG_HEADER_SIZE > size) {
        printf("\e[31m[ERROR] This message is too large (%lu bits / %zu "
               "bits)\e[0m\n",
               len + STEG_HEADER_SIZE, size);
        return *this;
    }
    printf("LENGTH: %d\n", len);
    for (uint8_t i = 0; i < STEG_HEADER_SIZE; i++) {
        data[i] &= 0xFE;
        data[i] |= (len >> (STEG_HEADER_SIZE - 1 - i)) & 1UL;
    }

    for (uint32_t i = 0; i < len; i++) {
        data[i + STEG_HEADER_SIZE] &= 0xFE;
        data[i + STEG_HEADER_SIZE] |= (message[i / 8] >> ((len - 1 - i) % 8)) & 1;
    }
    return *this;
}

Image& Image::decodemessage(char* buffer, size_t* messageLength) {
    uint32_t len = 0;
    for (uint8_t i = 0; i < STEG_HEADER_SIZE; i++) {
        len = (len << 1) | (data[i] & 1);
    }

    *messageLength = len / 8;

    for (uint32_t i = 0; i < len; i++) {
        buffer[i / 8] = (buffer[i / 8] << 1) | (data[i + STEG_HEADER_SIZE] & 1);
    }

    return *this;
}

Image& Image::diffmap_scale(Image& img, uint8_t scl) {
    int compare_width = fmin(w, img.w);
    int compare_height = fmin(h, img.h);
    int compare_channels = fmin(channels, img.channels);
    uint8_t largest = 0;
    for (uint32_t i = 0; i < compare_height; i++) {
        for (uint32_t j = 0; j < compare_width; j++) {
            for (uint8_t k = 0; k < compare_channels; k++) {
                data[(i * w + j) * channels + k] =
                    BYTE_BOUND(abs(data[(i * w + j) * channels + k] - img.data[(i * img.w + j) * img.channels + k]));
                largest = fmax(largest, data[(i * w + j) * channels + k]);
            }
        }
    }
    scl = 255 / fmax(1, fmax(scl, largest));
    for (int i = 0; i < size; i++) {
        data[i] *= scl;
    }
    return *this;
}
Image& Image::diffmap(Image& img) {
    int compare_width = fmin(w, img.w);
    int compare_height = fmin(h, img.h);
    int compare_channels = fmin(channels, img.channels);
    for (uint32_t i = 0; i < compare_height; i++) {
        for (uint32_t j = 0; j < compare_width; j++) {
            for (uint8_t k = 0; k < compare_channels; k++) {
                data[(i * w + j) * channels + k] =
                    BYTE_BOUND(abs(data[(i * w + j) * channels + k] - img.data[(i * img.w + j) * img.channels + k]));
            }
        }
    }
    return *this;
}

Image& Image::std_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                         uint32_t cc, bool normalize) {
    double new_data[w * h];
    uint64_t center = cr * ker_w + cc;
    for (uint64_t k = channel; k < size; k += channels) {
        double c = 0;
        for (long i = -((long)cr); i < (long)ker_h - cr; ++i) {
            long row = ((long)k / channels) / w - i;
            if (row < 0 || row > h - 1) {
                continue;
            }
            for (long j = -((long)cc); j < (long)ker_w - cc; ++j) {
                long col = ((long)k / channels) % w - j;
                if (col < 0 || col > w - 1) {
                    continue;
                }
                c += ker[center + i * (long)ker_w + j] * data[(row * w + col) * channels + channel];
            }
        }
        new_data[k / channels] = c;
    }
    if (normalize) {
        double mx = -INFINITY, mn = INFINITY;
        for (uint64_t k = channel; k < size; k += channels) {
            mx = fmax(mx, new_data[k / channels]);
            mn = fmin(mn, new_data[k / channels]);
        }
        for (uint64_t k = channel; k < size; k += channels) {
            new_data[k / channels] = (255.0 * (new_data[k / channels] - mn) / (mx - mn));
        }
    }

    for (uint64_t k = channel; k < size; k += channels) {
        data[k] = (uint8_t)BYTE_BOUND(round(new_data[k / channels]));
    }
    return *this;
}

Image& Image::std_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                           uint32_t cc, bool normalize) {
    double new_data[w * h];
    uint64_t center = cr * ker_w + cc;
    for (uint64_t k = channel; k < size; k += channels) {
        double c = 0;
        for (long i = -((long)cr); i < (long)ker_h - cr; i++) {
            long row = ((long)k / channels) / w - i;
            if (row < 0) {
                row = 0;
            } else if (row > h - 1) {
                row = h - 1;
            }
            for (long j = -((long)cc); j < (long)ker_w - cc; j++) {
                long col = ((long)k / channels) % w - j;
                if (col < 0) {
                    col = 0;
                } else if (col > w - 1) {
                    col = w - 1;
                }
                c += ker[center + i * (long)ker_w + j] * data[(row * w + col) * channels + channel];
            }
        }
        new_data[k / channels] = c;
    }
    if (normalize) {
        double mx = -INFINITY, mn = INFINITY;
        for (uint64_t k = channel; k < size; k += channels) {
            mx = fmax(mx, new_data[k / channels]);
            mn = fmin(mn, new_data[k / channels]);
        }
        for (uint64_t k = channel; k < size; k += channels) {
            new_data[k / channels] = (255.0 * (new_data[k / channels] - mn) / (mx - mn));
        }
    }

    for (uint64_t k = channel; k < size; k += channels) {
        data[k] = (uint8_t)BYTE_BOUND(round(new_data[k / channels]));
    }
    return *this;
}

Image& Image::std_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                  uint32_t cc, bool normalize) {
    double new_data[w * h];
    uint64_t center = cr * ker_w + cc;
    for (uint64_t k = channel; k < size; k += channels) {
        double c = 0;
        for (long i = -((long)cr); i < (long)ker_h - cr; i++) {
            long row = ((long)k / channels) / w - i;
            if (row < 0) {
                row = row % h + h;
            } else if (row > h - 1) {
                row %= h;
            }
            for (long j = -((long)cc); j < (long)ker_w - cc; j++) {
                long col = ((long)k / channels) % w - j;
                if (col < 0) {
                    col = col % w + w;
                } else if (col > w - 1) {
                    col %= w;
                }
                c += ker[center + i * (long)ker_w + j] * data[(row * w + col) * channels + channel];
            }
        }
        new_data[k / channels] = c;
    }
    if (normalize) {
        double mx = -INFINITY, mn = INFINITY;
        for (uint64_t k = channel; k < size; k += channels) {
            mx = fmax(mx, new_data[k / channels]);
            mn = fmin(mn, new_data[k / channels]);
        }
        for (uint64_t k = channel; k < size; k += channels) {
            new_data[k / channels] = (255.0 * (new_data[k / channels] - mn) / (mx - mn));
        }
    }

    for (uint64_t k = channel; k < size; k += channels) {
        data[k] = (uint8_t)BYTE_BOUND(round(new_data[k / channels]));
    }
    return *this;
}

Image& Image::flip_x() {
    uint8_t tmp[channels], *px1, *px2;
    for (int y = 0; y < h / 2; y++) {
        for (int x = 0; x < w; x++) {
            px1 = &data[(y * w + x) * channels];
            px2 = &data[((h - y - 1) * w + x) * channels];
            memcpy(tmp, px1, channels);
            memcpy(px1, px2, channels);
            memcpy(px2, tmp, channels);
        }
    }
    return *this;
}
Image& Image::flip_y() {
    uint8_t tmp[channels], *px1, *px2;
    for (int x = 0; x < w / 2; x++) {
        for (int y = 0; y < h; y++) {
            px1 = &data[(y * w + x) * channels];
            px2 = &data[(y * w + (w - 1 - x)) * channels];
            memcpy(tmp, px1, channels);
            memcpy(px1, px2, channels);
            memcpy(px2, tmp, channels);
        }
    }
    return *this;
}

Image& Image::overlay(const Image& src, int x, int y) {
    uint8_t *srcPx, *dstPx;
    for (int sy = 0; sy < src.h; ++sy) {
        if (sy + y < 0) {
            continue;
        } else if (sy + y >= h) {
            break;
        }
        for (int sx = 0; sx < src.w; ++sx) {
            if (sx + x < 0) {
                continue;
            } else if (sx + x >= w) {
                break;
            }
            srcPx = &src.data[(sx + sy * src.w) * src.channels];
            dstPx = &data[(sx + x + (sy + y) * w) * channels];

            float srcAlpha = src.channels < 4 ? 1 : srcPx[3] / 255.f;
            float dstAlpha = src.channels < 4 ? 1 : dstPx[3] / 255.f;

            if (srcAlpha > .99 && dstAlpha > .99) {
                if (src.channels >= channels)
                    memcpy(dstPx, srcPx, channels); // might have different channels
                else
                    memset(dstPx, srcPx[0],
                           channels); // requires fix for src with k>=2 channels and dest with n>k
            } else {
                float outAlpha = srcAlpha + dstAlpha * (1 - srcAlpha);
                if (outAlpha < .01) {
                    memset(dstPx, 0, channels);
                } else {
                    for (int chnl = 0; chnl < channels; chnl++) {
                        dstPx[chnl] = (uint8_t)BYTE_BOUND(
                            (srcPx[chnl] / 255.f * srcAlpha + dstPx[chnl] / 255.f * dstAlpha * (1 - srcAlpha)) /
                            outAlpha * 255.f);
                    }
                    if (channels > 3) {
                        dstPx[3] = (uint8_t)BYTE_BOUND(outAlpha * 255.f);
                    }
                }
            }
        }
    }
    return *this;
}

Image& Image::overlay_text(const char* txt, const Font& font, int x, int y, uint8_t r, uint8_t g, uint8_t b,
                           uint8_t a) {
    size_t len = strlen(txt);
    SFT_Char c;
    int32_t dx, dy;
    uint8_t *dstPx, srcPx, color[4] = {r, g, b, a};

    for (size_t i = 0; i < len; i++) {
        if (sft_char(&font.sft, txt[i], &c) != 0) {
            printf("\e[31m[ERROR] Font is missing character '%c'\e[0m\n", txt[i]);
            continue;
        }
        for (uint16_t sy = 0; sy < c.height; sy++) {
            dy = sy + y + c.y;
            if (dy < 0)
                continue;
            else if (dy >= h)
                break;
            for (uint16_t sx = 0; sx < c.width; sx++) {
                dx = sx + x + c.x;
                if (dx < 0)
                    continue;
                else if (dx >= w)
                    break;
                dstPx = &data[(dx + dy * w) * channels];
                srcPx = c.image[sx + sy * c.width];

                if (srcPx != 0) {
                    float srcAlpha = (srcPx / 255.f) * (a / 255.f);
                    float dstAlpha = channels < 4 ? 1 : dstPx[3] / 255.f;
                    if (srcAlpha > .99 && dstAlpha > .99) {
                        memcpy(dstPx, color, channels); // might have different channels
                    } else {
                        float outAlpha = srcAlpha + dstAlpha * (1 - srcAlpha);
                        if (outAlpha < .01) {
                            memset(dstPx, 0, channels);
                        } else {
                            for (int chnl = 0; chnl < channels; chnl++) {
                                dstPx[chnl] = (uint8_t)BYTE_BOUND(
                                    (color[chnl] / 255.f * srcAlpha + dstPx[chnl] / 255.f * dstAlpha * (1 - srcAlpha)) /
                                    outAlpha * 255.f);
                            }
                            if (channels > 3) {
                                dstPx[3] = (uint8_t)BYTE_BOUND(outAlpha * 255.f);
                            }
                        }
                    }
                }
            }
        }
        x += c.advance;
        free(c.image);
    }
    return *this;
}

Image& Image::crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch) {
    size = w * h * channels;
    uint8_t* croppedImage = new uint8_t[cw * ch * channels];

    memset(croppedImage, 0, size);
    for (uint16_t y = 0; y < ch; y++) {
        if (y + cy >= h)
            break;

        for (uint16_t x = 0; x < cw; x++) {
            if (x + cx >= w)
                break;
            memcpy(&croppedImage[(x + y * cw) * channels], &data[(x + cx + (y + cy) * w) * channels], channels);
        }
    }

    w = cw;
    h = ch;
    delete data;
    data = croppedImage;
    croppedImage = nullptr;
    return *this;
}

uint32_t Image::rev(uint32_t n, uint32_t a) {
    uint8_t max_bits = (uint8_t)ceil(log2(n));
    uint32_t reversed_a = 0;
    for (uint8_t i = 0; i < max_bits; i++) {
        if (a & (1 << i)) {
            reversed_a |= (1 << (max_bits - 1 - i));
        }
    }
    return reversed_a;
}
void Image::bit_rev(uint32_t n, std::complex<double> a[], std::complex<double>* A) {
    for (uint32_t i = 0; i < n; i++) {
        A[rev(n, i)] = a[i];
    }
}

void Image::fft(uint32_t n, std::complex<double> x[], std::complex<double>* X) {
    // x in standard order
    if (x != X) {
        memcpy(X, x, n * sizeof(std::complex<double>));
    }

    // Gentleman-Sande butterfly
    uint32_t sub_probs = 1;
    uint32_t sub_prob_size = n;
    uint32_t half;
    uint32_t i;
    uint32_t j_begin;
    uint32_t j_end;
    uint32_t j;
    std::complex<double> w_step;
    std::complex<double> w;
    std::complex<double> tmp1, tmp2;
    while (sub_prob_size > 1) {
        half = sub_prob_size >> 1;
        w_step = std::complex<double>(cos(-2 * M_PI / sub_prob_size), sin(-2 * M_PI / sub_prob_size));
        for (i = 0; i < sub_probs; ++i) {
            j_begin = i * sub_prob_size;
            j_end = j_begin + half;
            w = std::complex<double>(1, 0);
            for (j = j_begin; j < j_end; ++j) {
                tmp1 = X[j];
                tmp2 = X[j + half];
                X[j] = tmp1 + tmp2;
                X[j + half] = (tmp1 - tmp2) * w;
                w *= w_step;
            }
        }
        sub_probs <<= 1;
        sub_prob_size = half;
    }
    // X in bit reversed order
}
void Image::ifft(uint32_t n, std::complex<double> X[], std::complex<double>* x) {
    // X in bit reversed order
    if (X != x) {
        memcpy(x, X, n * sizeof(std::complex<double>));
    }

    // Cooley-Tukey butterfly
    uint32_t sub_probs = n >> 1;
    uint32_t sub_prob_size;
    uint32_t half = 1;
    uint32_t i;
    uint32_t j_begin;
    uint32_t j_end;
    uint32_t j;
    std::complex<double> w_step;
    std::complex<double> w;
    std::complex<double> tmp1, tmp2;
    while (half < n) {
        sub_prob_size = half << 1;
        w_step = std::complex<double>(cos(2 * M_PI / sub_prob_size), sin(2 * M_PI / sub_prob_size));
        for (i = 0; i < sub_probs; i++) {
            j_begin = i * sub_prob_size;
            j_end = j_begin + half;
            w = std::complex<double>(1, 0);
            for (j = j_begin; j < j_end; ++j) {
                tmp1 = x[j];
                tmp2 = w * x[j + half];
                x[j] = tmp1 + tmp2;
                x[j + half] = tmp1 - tmp2;
                w *= w_step;
            }
        }
        sub_probs >>= 1;
        half = sub_prob_size;
    }
    for (uint32_t i = 0; i < n; ++i) {
        x[i] /= n;
    }
    // x in standard order
}
void Image::dft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double>* X) {
    // x in row-major & standard order
    std::complex<double>* intermediate = new std::complex<double>[m * n];

    // rows
    for (uint32_t i = 0; i < m; i++) {
        fft(n, x + i * n, intermediate + i * n);
    }

    // cols
    for (uint32_t j = 0; j < n; j++) {
        for (uint32_t i = 0; i < m; i++) {
            X[j * m + i] = intermediate[i * n + j];
        }
        fft(m, X + j * m, X + j * m);
    }

    delete[] intermediate;
    // X in column-major & bit-reversed (int rows then columns)
}
void Image::idft_2D(uint32_t m, uint32_t n, std::complex<double> X[], std::complex<double>* x) {
    // X in column-major & bit-reversed (in rows then columns)
    std::complex<double>* intermediate = new std::complex<double>[m * n];
    // cols
    for (uint32_t j = 0; j < n; ++j) {
        ifft(m, X + j * m, intermediate + j * m);
    }
    // rows
    for (uint32_t i = 0; i < m; ++i) {
        for (uint32_t j = 0; j < n; ++j) {
            x[i * n + j] = intermediate[j * m + i]; // row-major <-- col-major
        }
        ifft(n, x + i * n, x + i * n);
    }
    delete[] intermediate;
    // x in row-major & standard order
}

void Image::pad_kernel(uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc, uint32_t pw, uint32_t ph,
                       std::complex<double>* pad_ker) {
    // padded so center of kernal is at top left
    for (long i = -((long)cr); i < (long)ker_h - cr; i++) {
        uint32_t r = (i < 0) ? i + ph : i;
        for (long j = -((long)cc); j < (long)ker_w - cc; j++) {
            uint32_t c = (j < 0) ? j + pw : j;
            pad_ker[r * pw + c] = std::complex<double>(ker[(i + cr) * ker_w + (j + cc)], 0);
        }
    }
}
void Image::pointwise_product(uint64_t l, std::complex<double> a[], std::complex<double> b[], std::complex<double>* p) {
    for (uint64_t k = 0; k < l; ++k) {
        p[k] = a[k] * b[k];
    }
}

std::complex<double>* Image::fd_convolve_clamp_to_zero_raw(uint8_t channel, uint32_t ker_w, uint32_t ker_h,
                                                           double ker[], uint32_t cr, uint32_t cc) {
    // calculate paddina
    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + ker_w - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + ker_h - 1)));
    uint64_t psize = pw * ph;

    // pad image
    std::complex<double>* pad_img = new std::complex<double>[psize];
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            pad_img[i * pw + j] = std::complex<double>(data[(i * w + j) * channels + channel], 0);
        }
    }

    // pad kernal
    std::complex<double>* pad_ker = new std::complex<double>[psize];
    pad_kernel(ker_w, ker_h, ker, cr, cc, pw, ph, pad_ker);

    // convolution
    dft_2D(ph, pw, pad_img, pad_img);
    dft_2D(ph, pw, pad_ker, pad_ker);
    pointwise_product(psize, pad_img, pad_ker, pad_img);
    idft_2D(ph, pw, pad_img, pad_img);

    return pad_img;
}
Image& Image::fd_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                        uint32_t cc, bool normalize) {
    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + ker_w - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + ker_h - 1)));
    std::complex<double>* pad_img = fd_convolve_clamp_to_zero_raw(channel, ker_w, ker_h, ker, cr, cc);
    if (normalize) {
        double mx = -INFINITY, mn = INFINITY;
        for (uint32_t i = 0; i < h; i++) {
            for (uint32_t j = 0; j < w; j++) {
                mx = fmax(mx, pad_img[i * pw + j].real());
                mn = fmin(mn, pad_img[i * pw + j].real());
            }
        }
        for (uint32_t i = 0; i < h; i++) {
            for (uint32_t j = 0; j < w; j++) {
                pad_img[i * pw + j].real(255.0 * (pad_img[i * pw + j].real() - mn) / (mx - mn));
            }
        }
    }
    // update pixel data
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] = (uint8_t)BYTE_BOUND(round(pad_img[i * pw + j].real()));
        }
    }

    delete[] pad_img;
    return *this;
}

std::complex<double>* Image::fd_convolve_clamp_to_border_raw(uint8_t channel, uint32_t ker_w, uint32_t ker_h,
                                                             double ker[], uint32_t cr, uint32_t cc) {
    // calculate padding
    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + ker_w - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + ker_h - 1)));
    uint64_t psize = pw * ph;

    // pad image
    std::complex<double>* pad_img = new std::complex<double>[psize];
    for (uint32_t i = 0; i < ph; i++) {
        uint32_t r = (i < h) ? i : ((i < h + cr ? h - 1 : 0));
        for (uint32_t j = 0; j < pw; j++) {
            uint32_t c = (j < w) ? j : ((j < w + cc ? w - 1 : 0));
            pad_img[i * pw + j] = std::complex<double>(data[(r * w + c) * channels + channel], 0);
        }
    }

    // pad kernal
    std::complex<double>* pad_ker = new std::complex<double>[psize];
    pad_kernel(ker_w, ker_h, ker, cr, cc, pw, ph, pad_ker);

    // convolution
    dft_2D(ph, pw, pad_img, pad_img);
    dft_2D(ph, pw, pad_ker, pad_ker);
    pointwise_product(psize, pad_img, pad_ker, pad_img);
    idft_2D(ph, pw, pad_img, pad_img);

    delete[] pad_ker;

    return pad_img;
}

Image& Image::fd_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                          uint32_t cc, bool normalize) {
    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + ker_w - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + ker_h - 1)));
    std::complex<double>* pad_img = fd_convolve_clamp_to_border_raw(channel, ker_w, ker_h, ker, cr, cc);
    if (normalize) {
        double mx = -INFINITY, mn = INFINITY;
        for (uint32_t i = 0; i < h; i++) {
            for (uint32_t j = 0; j < w; j++) {
                mx = fmax(mx, pad_img[i * pw + j].real());
                mn = fmin(mn, pad_img[i * pw + j].real());
            }
        }
        for (uint32_t i = 0; i < h; i++) {
            for (uint32_t j = 0; j < w; j++) {
                pad_img[i * pw + j].real(255.0 * (pad_img[i * pw + j].real() - mn) / (mx - mn));
            }
        }
    }
    // update pixel data
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] = (uint8_t)BYTE_BOUND(round(pad_img[i * pw + j].real()));
        }
    }

    delete[] pad_img;
    return *this;
}

std::complex<double>* Image::fd_convolve_cyclic_raw(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[],
                                                    uint32_t cr, uint32_t cc) {
    // calculate padding
    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + ker_w - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + ker_h - 1)));
    uint64_t psize = pw * ph;

    // pad image
    std::complex<double>* pad_img = new std::complex<double>[psize];
    for (uint32_t i = 0; i < ph; i++) {
        uint32_t r = (i < h) ? i : (i < h + cr ? i % h : h - ph + i);
        for (uint32_t j = 0; j < pw; j++) {
            uint32_t c = (j < w) ? j : (j < w + cc ? j % w : w - pw + j);
            pad_img[i * pw + j] = std::complex<double>(data[(i * w + j) * channels + channel], 0);
        }
    }

    // pad kernal
    std::complex<double>* pad_ker = new std::complex<double>[psize];
    pad_kernel(ker_w, ker_h, ker, cr, cc, pw, ph, pad_ker);

    // convolution
    dft_2D(ph, pw, pad_img, pad_img);
    dft_2D(ph, pw, pad_ker, pad_ker);
    pointwise_product(psize, pad_img, pad_ker, pad_img);
    idft_2D(ph, pw, pad_img, pad_img);

    return pad_img;
}
Image& Image::fd_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                 uint32_t cc, bool normalize) {
    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + ker_w - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + ker_h - 1)));
    std::complex<double>* pad_img = fd_convolve_cyclic_raw(channel, ker_w, ker_h, ker, cr, cc);

    if (normalize) {
        double mx = -INFINITY, mn = INFINITY;
        for (uint32_t i = 0; i < h; i++) {
            for (uint32_t j = 0; j < w; j++) {
                mx = fmax(mx, pad_img[i * pw + j].real());
                mn = fmin(mn, pad_img[i * pw + j].real());
            }
        }
        for (uint32_t i = 0; i < h; i++) {
            for (uint32_t j = 0; j < w; j++) {
                pad_img[i * pw + j].real(255.0 * (pad_img[i * pw + j].real() - mn) / (mx - mn));
            }
        }
    }

    // update pixel data
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] = BYTE_BOUND((uint8_t)round(pad_img[i * pw + j].real()));
        }
    }

    delete[] pad_img;
    return *this;
}

Image& Image::convolve_linear(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc,
                              bool normalize) {
    if (ker_w * ker_h > 224) {
        return fd_convolve_clamp_to_zero(channel, ker_w, ker_h, ker, cr, cc, normalize);
    } else {
        return std_convolve_clamp_to_zero(channel, ker_w, ker_h, ker, cr, cc, normalize);
    }
}
Image& Image::convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr,
                                       uint32_t cc, bool normalize) {
    if (ker_w * ker_h > 224) {
        return fd_convolve_clamp_to_border(channel, ker_w, ker_h, ker, cr, cc, normalize);
    } else {
        return std_convolve_clamp_to_border(channel, ker_w, ker_h, ker, cr, cc, normalize);
    }
}
Image& Image::convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc,
                              bool normalize) {
    if (ker_w * ker_h > 224) {
        return fd_convolve_cyclic(channel, ker_w, ker_h, ker, cr, cc, normalize);
    } else {
        return std_convolve_cyclic(channel, ker_w, ker_h, ker, cr, cc, normalize);
    }
}
Image& Image::brightness(uint8_t channel, double brightness_delta) {
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] =
                (uint8_t)BYTE_BOUND(data[(i * w + j) * channels + channel] - 128 + 128 + brightness_delta);
        }
    }
    return *this;
}
Image& Image::contrast(uint8_t channel, double contrast_delta) {
    double F = 259.0 * (contrast_delta + 255) / (255.0 * (259 - contrast_delta));
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] =
                (uint8_t)BYTE_BOUND(round(F * (data[(i * w + j) * channels + channel] - 128) + 128));
        }
    }
    return *this;
}
Image& Image::shade_h() {
    double gaussian_blur[] = {1 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 4 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 1 / 16.0};
    double scharr_y[] = {47, 162, 47, 0, 0, 0, -47, -162, -47};
    grayscale_avg();
    convolve_linear(0, 3, 3, gaussian_blur, 1, 1);
    convolve_linear(0, 3, 3, scharr_y, 1, 1, true);
    if (channels >= 3) {
        for (int i = 0; i < size; i += channels) {
            data[i + 1] = data[i + 2] = data[i];
        }
    }
    return *this;
}
Image& Image::shade_v() {
    double gaussian_blur[] = {1 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 4 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 1 / 16.0};
    double scharr_x[] = {47, 0, -47, 162, 0, -162, 47, 0, -47};
    grayscale_avg();
    convolve_clamp_to_border(0, 3, 3, gaussian_blur, 1, 1);
    convolve_clamp_to_border(0, 3, 3, scharr_x, 1, 1, true);
    if (channels >= 3) {
        for (int i = 0; i < size; i += channels) {
            data[i + 1] = data[i + 2] = data[i];
        }
    }
    return *this;
}
Image& Image::shade() {

    grayscale_avg();
    Image itmd_h(w, h, 1), itmd_v(w, h, 1);
    for (uint64_t k = 0; k < w * h; k++) {
        itmd_h.data[k] = data[k * channels];
        itmd_v.data[k] = data[k * channels];
    }
    itmd_h.shade_h();
    itmd_v.shade_v();
    double h, v;
    for (uint64_t i = 0; i < size; i += channels) {
        h = itmd_h.data[i / channels];
        v = itmd_v.data[i / channels];
        data[i] = (uint8_t)BYTE_BOUND(sqrt(h * h + v * v));
        if (channels >= 3) {
            data[i + 1] = data[i + 2] = data[i];
        }
    }
    return *this;
}
Image& Image::edge(bool gradient, double detail_threshold) {

    uint32_t pw = 1 << ((uint8_t)ceil(log2(w + 3 - 1)));
    uint32_t ph = 1 << ((uint8_t)ceil(log2(h + 3 - 1)));

    double gaussian_blur[] = {1 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 4 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 1 / 16.0};
    double scharr_x[] = {47, 0, -47, 162, 0, -162, 47, 0, -47};
    double scharr_y[] = {47, 162, 47, 0, 0, 0, -47, -162, -47};
    grayscale_avg();
    convolve_linear(0, 3, 3, gaussian_blur, 1, 1);
    // Image itmd_y(w, h, 1);
    // for (uint64_t k = 0; k < w * h; k++) {
    //     itmd_y.data[k] = data[k * channels];
    // }
    Image itmd_y(*this);
    std::complex<double>* gx = fd_convolve_clamp_to_border_raw(0, 3, 3, scharr_x, 1, 1);
    std::complex<double>* gy = itmd_y.fd_convolve_clamp_to_border_raw(0, 3, 3, scharr_y, 1, 1);

    // Normalization
    double* g = new double[w * h];
    double* theta = new double[w * h];
    for (uint64_t i = 0; i < h; i++) {
        for (uint64_t j = 0; j < w; j++) {
            g[i * w + j] =
                sqrt(gx[i * pw + j].real() * gx[i * pw + j].real() + gy[i * pw + j].real() * gy[i * pw + j].real());
            theta[i * w + j] = atan2(gy[i * pw + j].real(), gx[i * pw + j].real());
        }
    }

    // make images
    double mx = -INFINITY, mn = INFINITY;
    for (uint64_t i = 2; i < h - 2; i++) {
        for (uint64_t j = 2; j < w - 2; j++) {
            mx = fmax(mx, g[i * w + j]);
            mn = fmin(mn, g[i * w + j]);
        }
    }
    double v;
    for (uint64_t i = 0; i < h; i++) {
        for (uint64_t j = 0; j < w; j++) {
            if (mn == mx) {
                v = 0;
            } else {
                v = (g[i * w + j] - mn) / (mx - mn);
                v = (v < detail_threshold || v > 1 ? 0 : v);
            }
            data[(i * w + j) * channels] = (uint8_t)(255 * v);
        }
    }
    if (!gradient && channels >= 3) {
        for (uint64_t i = 0; i < h; i++) {
            for (uint64_t j = 0; j < w; j++) {
                data[(i * w + j) * channels + 1] = data[(i * w + j) * channels + 2] = data[(i * w + j) * channels];
            }
        }
    } else if (gradient) {
        std::cout << channels << "\n";
        Color c(0.0, 0.0, 0.0);
        for (uint64_t i = 0; i < h; i++) {
            for (uint64_t j = 0; j < w; j++) {
                double h, s;
                h = theta[i * w + j] * 180.0 / M_PI + 180.0;
                v = (g[i * w + j] - mn) / (mx - mn);
                v = (v < detail_threshold || v > 1 ? 0 : v);
                c.hsv_to_rgb(h, v, v);
                data[(i * w + j) * channels] = (uint8_t)BYTE_BOUND(round(c.r));
                data[(i * w + j) * channels + 1] = (uint8_t)BYTE_BOUND(round(c.g));
                data[(i * w + j) * channels + 2] = (uint8_t)BYTE_BOUND(round(c.b));
            }
        }
    }
    delete[] gx;
    delete[] gy;
    delete[] g;
    delete[] theta;

    return *this;
}

Image& Image::f_scale(uint32_t new_w, uint32_t new_h, bool linked, TwoDimInterp method) {
    printf("Channels: %d", channels);
    if (linked) {
        new_h = (uint32_t)round(((double)h) / w * new_w);
    }
    uint8_t* new_data = new uint8_t[new_w * new_h * channels];
    double r_old, c_old;
    if (method == TwoDimInterp::Nearest) {
        for (int r = 0; r < new_h; r++) {
            for (int c = 0; c < new_w; c++) {
                r_old = (double)r * h / new_h;
                c_old = (double)c * w / new_w;
                for (int cd = 0; cd < channels; cd++) {
                    new_data[(r * new_w + c) * channels + cd] =
                        data[((uint32_t)round(r_old) * w + (uint32_t)round(c_old)) * channels + cd];
                }
            }
        }
    } else if (method == TwoDimInterp::Bilinear) {
        double r_diff, c_diff;
        for (int r = 0; r < new_h; r++) {
            for (int c = 0; c < new_w; c++) {
                r_old = (double)r * h / new_h;
                c_old = (double)c * w / new_w;
                // r_diff = r_old - floor(r_old) < ceil(r_old) - r_old ? r_old - floor(r_old) : r_old - ceil(r_old);
                // c_diff = c_old - floor(c_old) < ceil(c_old) - c_old ? c_old - floor(c_old) : c_old - ceil(c_old);
                for (int cd = 0; cd < channels; cd++) {
                    if ((r_old == floor(r_old) && c_old == floor(c_old)) || (ceil(r_old) == h && ceil(c_old) == w) ||
                        (r_old == floor(r_old) && ceil(c_old) == w) || (ceil(r_old) == h && c_old == floor(c_old))) {
                        new_data[(r * new_w + c) * channels + cd] =
                            data[(uint32_t)round(floor(r_old) * w + floor(c_old)) * channels + cd];
                    } else if (c_old == floor(c_old) || ceil(c_old) == w) {
                        uint32_t y1 = (uint32_t)round(floor(r_old)), y2 = (uint32_t)round(ceil(r_old));
                        new_data[(r * new_w + c) * channels + cd] =
                            data[(uint32_t)round(y1 * w + floor(c_old)) * channels + cd] * (double)(y2 - r_old) /
                                (y2 - y1) +
                            data[(uint32_t)round(y2 * w + floor(c_old)) * channels + cd] * (double)(r_old - y1) /
                                (y2 - y1);
                    } else if (r_old == floor(r_old) || ceil(r_old) == h) {
                        uint32_t x1 = (uint32_t)round(floor(c_old)), x2 = (uint32_t)round(ceil(c_old));
                        new_data[(r * new_w + c) * channels + cd] =
                            data[(uint32_t)round(floor(r_old) * w + x1) * channels + cd] * (double)(x2 - c_old) /
                                (x2 - x1) +
                            data[(uint32_t)round(floor(r_old) * w + x2) * channels + cd] * (double)(c_old - x1) /
                                (x2 - x1);
                    }
                    // else if (abs(r_diff) < 0.005 && abs(c_diff) < 0.005) {
                    //     new_data[(r * new_w + c) * channels + cd] =
                    //     data[(((uint32_t)round(floor(ceil(r_old) - r_diff))) * w +
                    //     ((uint32_t)round(floor(ceil(c_old) - c_diff))))*channels + cd];
                    // }
                    else {
                        uint32_t y1 = floor(r_old), y2 = ceil(r_old), x1 = floor(c_old), x2 = ceil(c_old);
                        new_data[(r * new_w + c) * channels + cd] =
                            (data[(y1 * w + x1) * channels + cd] * (double)(x2 - c_old) * (double)(y2 - r_old) +
                             (double)data[(y2 * w + x1) * channels + cd] * (double)(c_old - x1) * (double)(y2 - r_old) +
                             (double)data[(y1 * w + x2) * channels + cd] * (double)(x2 - c_old) * (double)(r_old - y1) +
                             (double)data[(y2 * w + x2) * channels + cd] * (double)(c_old - x1) *
                                 (double)(r_old - y1)) /
                            ((double)(x2 - x1) * (y2 - y1));
                    }
                }
            }
        }
    } else {
        throw std::invalid_argument("The scale method specified is not yet supported\n");
    }
    w = new_w;
    h = new_h;
    delete[] data;
    data = new_data;
    size = new_w * new_h * channels;

    return *this;
}

Image& Image::invert_color(uint8_t channel) {
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] = 255 - data[(i * w + j) * channels + channel];
        }
    }
    return *this;
}

Image& Image::gamma(uint8_t channel, double gamma_delta) {
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            data[(i * w + j) * channels + channel] =
                255.0 * pow(data[(i * w + j) * channels + channel] / 255.0, 1 / gamma_delta);
        }
    }

    return *this;
}

Image& Image::color_reduce(bool error_diffusion) {
    int r_diff, g_diff, b_diff, a[] = {0, 0, 0};
    std::pair<int, int*> min(0x3f3f3f3f, a);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            min.first = 0x3f3f3f3f;
            for (int color = 0; color < 8; color++) {
                r_diff = data[(i * w + j) * channels] - (((1 << 2) & color) ? 255 : 0);
                g_diff = data[(i * w + j) * channels + (uint8_t)fmin(1, channels - 1)] - (((1 << 1) & color) ? 255 : 0);
                b_diff = data[(i * w + j) * channels + (uint8_t)fmin(2, channels - 1)] - ((1 & color) ? 255 : 0);
                if (r_diff * r_diff + g_diff * g_diff + b_diff * b_diff < min.first) {
                    min.first = r_diff * r_diff + g_diff * g_diff + b_diff * b_diff;
                    min.second[0] = r_diff;
                    min.second[1] = g_diff;
                    min.second[2] = b_diff;
                }
            }
            if (channels < 3) {
                data[(i * w + j) * channels] -= min.second[0];
            } else {
                data[(i * w + j) * channels] -= min.second[0];
                data[(i * w + j) * channels + (uint8_t)fmin(1, channels - 1)] -= min.second[1];
                data[(i * w + j) * channels + (uint8_t)fmin(2, channels - 1)] -= min.second[2];
            }
            if (error_diffusion) {
                // Floyd-Steinberg Error Diffusion
                for (int c_cn = 0; c_cn < fmin(3, channels); c_cn++) {
                    if (j + 1 < w)
                        data[(i * w + (j + 1)) * channels + c_cn] =
                            BYTE_BOUND(data[(i * w + (j + 1)) * channels + c_cn] + 7.0 / 16 * min.second[c_cn]);
                    if (i + 1 < h && j + 1 < w)
                        data[((i + 1) * w + (j + 1)) * channels + c_cn] =
                            BYTE_BOUND(data[((i + 1) * w + (j + 1)) * channels + c_cn] + 1.0 / 16 * min.second[c_cn]);
                    if (i + 1 < h)
                        data[((i + 1) * w + j) * channels + c_cn] =
                            BYTE_BOUND(data[((i + 1) * w + j) * channels + c_cn] + 5.0 / 16 * min.second[c_cn]);
                    if (i + 1 < h && j > 0)
                        data[((i + 1) * w + (j - 1)) * channels + c_cn] =
                            BYTE_BOUND(data[((i + 1) * w + (j - 1)) * channels + c_cn] + 3.0 / 16 * min.second[c_cn]);
                }
            }
        }
    }
    return *this;
}
Image& Image::color_ramp(std::vector<std::pair<double, Color>> points, OneDimInterp method) {
    // assuming all points are sorted
    // std::sort(points.begin(), points.end());
    grayscale_avg();
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            int ceil_idx =
                std::upper_bound(points.begin(), points.end(),
                                 std::make_pair((double)data[(i * w + j) * channels] / 255, Color(-1, -1, -1))) -
                points.begin();
            if (ceil_idx == 0 && (method == OneDimInterp::Constant || method == OneDimInterp::Linear)) {
                if (channels < 3) {
                    data[(i * w + j) * channels] = points[0].second.luminance();
                } else {
                    for (int clr = 0; clr < 3; clr++) {
                        data[(i * w + j) * channels + clr] = points[0].second.get(clr);
                    }
                }
            } else if (method == OneDimInterp::Constant ||
                       (method == OneDimInterp::Linear && ceil_idx == points.size())) {
                if (channels < 3) {
                    data[(i * w + j) * channels] = points[ceil_idx - 1].second.luminance();
                } else {
                    for (int clr = 0; clr < 3; clr++) {
                        data[(i * w + j) * channels + clr] = points[ceil_idx - 1].second.get(clr);
                    }
                }
            } else if (method == OneDimInterp::Linear) {
                std::pair<double, Color> x1 = points[ceil_idx - 1], x2 = points[ceil_idx];
                double x = (double)data[(i * w + j) * channels] / 255.0;
                if (channels < 3) {
                    data[(i * w + j) * channels] = (x - x1.first) / (x2.first - x1.first) * x2.second.luminance() +
                                                   (x2.first - x) / (x2.first - x1.first) * x1.second.luminance();
                } else {
                    for (int clr = 0; clr < 3; clr++) {
                        x = (double)data[(i * w + j) * channels + clr] / 255.0;
                        data[(i * w + j) * channels + clr] =
                            (x - x1.first) / (x2.first - x1.first) * x2.second.get(clr) +
                            (x2.first - x) / (x2.first - x1.first) * x1.second.get(clr);
                    }
                }
            } else {
                throw std::invalid_argument("The scale method specified is not yet supported\n");
            }
        }
    }
    return *this;
}
Image& Image::preview_color_ramp(std::vector<std::pair<double, Color>> points, OneDimInterp method) const {
    Image* ret = new Image(256, 20, 3);
    for (int i = 0; i < ret->h; i++) {
        for (int j = 0; j < ret->w; j++) {
            for (int cd = 0; cd < 3; cd++) {
                ret->data[(i * ret->w + j) * ret->channels + cd] = j;
            }
        }
    }
    ret->write("images/preview_start.png");
    ret->color_ramp(points, method);
    return *ret;
}

// Notes: will not add alpha channel even if fill does have alpha channel != 255
Image& Image::translate(int x, int y, Color fill) {
    uint8_t* new_data = new uint8_t[size];
    for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
            for (uint8_t cd = 0; cd < channels; cd++) {
                if (i - y < 0 || j - y < 0 || i - y >= h || j - x >= w) {
                    if (cd < 3 && channels >= 3)
                        new_data[(i * w + j) * channels + cd] = fill.get(cd);
                    else if (cd < 3)
                        new_data[(i * w + j) * channels + cd] = fill.luminance();
                    else
                        new_data[(i * w + j) * channels + cd] = fill.a;
                } else
                    new_data[(i * w + j) * channels + cd] =
                        data[((uint32_t)round(i - y) * w + (uint32_t)round(j - x)) * channels + cd];
            }
        }
    }
    delete[] data;
    data = new_data;
    return *this;
}

std::vector<Image*> Image::seperate_channels() {
    std::vector<Image*> v;
    for (int8_t cd = 0; cd < channels; cd++) {
        Image* sep_chn = new Image(w, h, 1);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                sep_chn->data[i * w + j] = data[(i * w + j) * channels + cd];
            }
        }
        v.push_back(sep_chn);
    }
    v.shrink_to_fit();
    return v;
}

// Notes: assuming imgs is in order of {r_img, b_img, g_img, a_img, useless channel...} and default values are 0 if not
// the same size;
Image& Image::combine_channels(std::vector<Image*> imgs, bool resize_to_fit, TwoDimInterp method) {
    uint32_t max_w = 0, max_h = 0;
    for (int i = 0; i < imgs.size(); i++) {
        max_w = fmax(max_w, imgs[i]->w);
        max_h = fmax(max_h, imgs[i]->h);
    }
    uint8_t* new_data = new uint8_t[max_w * max_h * imgs.size()];
    if (resize_to_fit) {
        for (int i = 0; i < imgs.size(); i++) {
            imgs[i]->f_scale(max_w, max_h, false, method);
        }
    }
    for (int img = 0; img < imgs.size(); img++) {
        for (int i = 0; i < max_h; i++) {
            for (int j = 0; j < max_w; j++) {
                new_data[(i * max_w + j) * imgs.size() + img] =
                    (i >= imgs[img]->h || j >= imgs[img]->w)
                        ? 0
                        : imgs[img]->data[(i * imgs[img]->w + j) * imgs[img]->channels];
            }
        }
    }
    w = max_w;
    h = max_h;
    channels = imgs.size();
    delete[] data;
    data = new_data;
    return *this;
}

// Notes:
// - default values are opaque (255)
// - requires testing for single-channel grayscale images (try 2 channel images)
Image& Image::set_alpha(Image& alph, bool resize_to_fit, TwoDimInterp method) {
    if (resize_to_fit) {
        alph.f_scale(w, h, false, method);
    }
    if (channels >= 4) {
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                data[(i * w + j) * channels + 3] =
                    (i >= alph.h || j >= alph.w) ? 255 : alph.data[(i * alph.w + j) * alph.channels];
            }
        }
    } else {
        uint8_t* new_data = new uint8_t[w * h * (channels + 1)];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                memcpy(&new_data[(i * w + j) * (channels + 1)], &data[(i * w + j) * channels], (size_t)channels);
                new_data[(i * w + j) * (channels + 1) + 3] =
                    (i >= alph.h || j >= alph.w) ? 255 : alph.data[(i * alph.w + j) * alph.channels];
            }
        }
        channels++;
        delete[] data;
        data = new_data;
    }
    return *this;
}
// Notes: alternative formula available
Image& Image::color_balance(Color lift, Color gamma, Color gain) {
    for (int i = 0; i < 3; i++) {
        lift.set(i, lift.get(i) / 255.0);
        gamma.set(i, gamma.get(i) / 255.0);
        gain.set(i, gain.get(i) / 255.0);
    }
    double x;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            for (int cd = 0; cd < fmin(channels, 3); cd++) {
                double x = data[(i * w + j) * channels + cd] / 255.0;
                data[(i * w + j) * channels + cd] =
                    BYTE_BOUND(round(pow(gain.get(cd) * (x + lift.get(cd) * (1 - x)), 1.0 / gamma.get(cd)) * 255.0));
                // data[(i * w + j)*channels + cd] = BYTE_BOUND(round(pow(x*(gain.get(cd) - lift.get(cd)) +
                // lift.get(cd), 1.0/gamma.get(cd)) * 255.0));
            }
        }
    }
    return *this;
}
// Notes: channels < 3 images require testing
Image& Image::histogram(bool inc_lum, int channel) {
    Image* hist = new Image(256, 256, 3);
    uint32_t cnt_clr[4][256] = {0};
    uint64_t max_cnt = 0;
    Color color(0, 0, 0);

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (channel < 0) {
                for (int cn = 0; cn < fmin(channels, 3); cn++) {
                    max_cnt = fmax(++cnt_clr[cn][data[(i * w + j) * channels + cn]], max_cnt);
                }
                if (inc_lum && channels >= 3) {
                    max_cnt = fmax(++cnt_clr[3][(uint32_t)round(color.luminance(data[(i * w + j) * channels],
                                                                                data[(i * w + j) * channels + 1],
                                                                                data[(i * w + j) * channels + 2]))],
                                   max_cnt);
                } else if (inc_lum) {
                    printf("Unable to calculate luminance for image with %d channels.\n", channels);
                }
            } else {
                max_cnt = fmax(++cnt_clr[0][data[(i * w + j) * channels + channel]], max_cnt);
            }
        }
    }

    for (int c = 0; c < 256; c++) {
        for (int cn = 0; cn < fmin(channels, channel < 0 ? 3.0 : 1.0); cn++) {
            for (int r = 0; r <= (double)cnt_clr[(channel < 0 ? cn : 0)][c] * 255.0 / (double)max_cnt; r++) {
                if (channels < 3 || channel > 0) {
                    hist->data[((255 - r) * 256 + c) * 3] = hist->data[((255 - r) * 256 + c) * 3 + 1] =
                        hist->data[((255 - r) * 256 + c) * 3 + 2] = 155;
                } else {
                    hist->data[((255 - r) * 256 + c) * 3 + cn] += 155;
                }
            }
        }
        if (inc_lum && channels >= 3) {
            for (int r = 0; r <= (double)cnt_clr[3][c] * 255.0 / (double)max_cnt; r++) {
                hist->set((255 - r), c, 0, std::clamp(hist->get((255 - r), c, 0) + 100, 0, 255));
                hist->set((255 - r), c, 1, std::clamp(hist->get((255 - r), c, 1) + 100, 0, 255));
                hist->set((255 - r), c, 2, std::clamp(hist->get((255 - r), c, 2) + 100, 0, 255));
            }
        }
    }
    return *hist;
}

// Notes: channels < 3 images require testing
Image& Image::histogram_lum() {
    Image* hist = new Image(256, 256, 3);
    uint32_t cnt_clr[256] = {0};
    uint64_t max_cnt = 0;
    Color color(0, 0, 0);

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (channels < 3) {
                max_cnt = fmax(++cnt_clr[data[(i * w + j) * channels]], max_cnt);
            } else {
                max_cnt = fmax(++cnt_clr[(uint32_t)round(color.luminance(data[(i * w + j) * channels],
                                                                         data[(i * w + j) * channels + 1],
                                                                         data[(i * w + j) * channels + 2]))],
                               max_cnt);
            }
        }
    }

    for (int c = 0; c < 256; c++) {
        for (int r = 0; r <= (double)cnt_clr[c] * 255.0 / (double)max_cnt; r++) {
            hist->data[((255 - r) * 256 + c) * 3] = hist->data[((255 - r) * 256 + c) * 3 + 1] =
                hist->data[((255 - r) * 256 + c) * 3 + 2] = 155;
        }
    }
    return *hist;
}
Image& Image::HSV(double hue_delta, double saturation_delta, double value_delta) {
    if (channels < 3) {
        printf("Given there are %d channels, hue and saturation adjustments will not take any effect.\n", channels);
    }
    Color c(0, 0, 0);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (channels < 3) {
                c = c.rgb_to_hsv(data[(i * w + j) * channels], data[(i * w + j) * channels],
                                 data[(i * w + j) * channels]);
                c.hsv_to_rgb(c.r + hue_delta, c.g + saturation_delta, c.b + value_delta);
                data[(i * w + j) * channels] = c.r;
            } else {
                c = c.rgb_to_hsv(data[(i * w + j) * channels], data[(i * w + j) * channels + 1],
                                 data[(i * w + j) * channels + 2]);
                c.hsv_to_rgb(c.r + hue_delta, std::clamp(c.g + saturation_delta, 0.0, 1.0),
                             std::clamp(c.b + value_delta, 0.0, 1.0));
                data[(i * w + j) * channels] = c.r;
                data[(i * w + j) * channels + 1] = c.g;
                data[(i * w + j) * channels + 2] = c.b;
            }
        }
    }
    return *this;
}
Image& Image::false_color(bool overwrite) {
    Color c(0, 0, 0);
    Color lookup[] = {Color(0, 0, 0),     Color(0, 0, 255),     Color(0, 127, 255),  Color(0, 255, 255),
                      Color(0, 255, 127), Color(127, 127, 127), Color(127, 255, 0),  Color(255, 255, 0),
                      Color(255, 127, 0), Color(255, 0, 0),     Color(255, 255, 255)};
    int endpoints[] = {0, 15, 58, 102, 117, 137, 153, 196, 239, 254, 255};
    Image* ret = this;
    if (!overwrite || channels < 3) {
        printf("yes\n");
        Image* tmp = new Image(w, h, 3);
        ret = tmp;
    }
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            double lum = data[(i * w + j) * channels];
            if (channels >= 3) {
                lum = c.luminance(data[(i * w + j) * channels], data[(i * w + j) * channels + 1],
                                  data[(i * w + j) * channels + 2]);
            }
            if (lum <= 0) {
                ret->data[(i * w + j) * 3] = lookup[0].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[0].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[0].b;
            } else if (lum <= 15) {
                ret->data[(i * w + j) * 3] = lookup[1].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[1].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[1].b;
            } else if (lum <= 58) {
                ret->data[(i * w + j) * 3] = lookup[2].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[2].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[2].b;
            } else if (lum <= 102) {
                ret->data[(i * w + j) * 3] = lookup[3].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[3].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[3].b;
            } else if (lum <= 117) {
                ret->data[(i * w + j) * 3] = lookup[4].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[4].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[4].b;
            } else if (lum <= 137) {
                ret->data[(i * w + j) * 3] = lookup[5].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[5].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[5].b;
            } else if (lum <= 153) {
                ret->data[(i * w + j) * 3] = lookup[6].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[6].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[6].b;
            } else if (lum <= 196) {
                ret->data[(i * w + j) * 3] = lookup[7].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[7].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[7].b;
            } else if (lum <= 239) {
                ret->data[(i * w + j) * 3] = lookup[8].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[8].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[8].b;
            } else if (lum <= 254) {
                ret->data[(i * w + j) * 3] = lookup[9].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[9].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[9].b;
            } else {
                ret->data[(i * w + j) * 3] = lookup[10].r;
                ret->data[(i * w + j) * 3 + 1] = lookup[10].g;
                ret->data[(i * w + j) * 3 + 2] = lookup[10].b;
            }
        }
    }
    if (overwrite) {
        data = ret->data;
        size = w * h * 3;
    }
    return *ret;
}
Image& Image::tone_correct(uint8_t midtones_start, uint8_t midtones_end, Adjustment shadow, Adjustment midtone,
                           Adjustment highlight) {
    if (channels < 3) {
        printf("Given there are %d channels, hue and saturation adjustments will not take any effect.\n", channels);
    }
    Color c(0, 0, 0);
    double s_fac, m_fac, h_fac, lum, l, hh, fl, fh;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (channels < 3) {
                c.r = c.g = c.b = data[(i * w + j) * channels];
            } else {
                c.r = data[(i * w + j) * channels];
                c.g = data[(i * w + j) * channels + 1];
                c.b = data[(i * w + j) * channels + 2];
            }
            lum = c.luminance();
            if (lum < midtones_start) {
                s_fac = lum / midtones_start;
                m_fac = 0.1 * lum / midtones_start;
                h_fac = 0.1 * lum / ((midtones_start + midtones_end) / 2);
            } else if (lum < midtones_end) {
                if (lum <= ((double)midtones_end - midtones_start) / 2.0) {
                    l = midtones_start;
                    hh = ((double)midtones_end - midtones_start) / 2.0;
                    fl = 1;
                    fh = 0.1;
                } else {
                    l = ((double)midtones_end - midtones_start) / 2.0;
                    hh = 255;
                    fl = 0.1;
                    fh = 0;
                }
                s_fac = ((lum - l) * fh + (hh - lum) * fl) / (hh - l);
                if (lum < (midtones_start + midtones_end) / 2.0) {
                    l = midtones_start;
                    hh = (midtones_start + midtones_end) / 2.0;
                    fl = 0.1;
                    fh = 1;
                } else {
                    l = (midtones_start + midtones_end) / 2.0;
                    hh = midtones_end;
                    fl = 1;
                    fh = 0.1;
                }
                m_fac = ((lum - l) * fh + (hh - lum) * fl) / (hh - l);

                if (lum <= ((double)midtones_start + midtones_end) / 2.0) {
                    l = 0;
                    hh = ((double)midtones_start + midtones_end) / 2.0;
                    fl = 0;
                    fh = 0.1;
                } else {
                    l = ((double)midtones_start + midtones_end) / 2.0;
                    hh = midtones_end;
                    fl = 0.1;
                    fh = 0.75;
                }
                h_fac = ((lum - l) * fh + (hh - lum) * fl) / (hh - l);
            } else {
                s_fac = ((255 - lum) * 0.1) / (255 - ((double)midtones_start + midtones_end) / 2.0);
                m_fac = ((255 - lum) * 0.1) / (255 - midtones_end);
                h_fac = ((lum - midtones_end) * 1 + (255.0 - lum) * 0.8) / (255 - midtones_end);
            }
            // printf("%f, %f, %f\n", s_fac, m_fac, h_fac);
            c = c.apply_adj_rgb(shadow, s_fac);
            c = c.apply_adj_rgb(midtone, m_fac);
            c = c.apply_adj_rgb(highlight, h_fac);
            if (channels < 3) {
                data[(i * w + j) * channels] = c.r;
            } else {
                data[(i * w + j) * channels] = c.r;
                data[(i * w + j) * channels + 1] = c.g;
                data[(i * w + j) * channels + 2] = c.b;
            }
        }
    }
    return *this;
}
Image& Image::rotate(double origin_x, double origin_y, double angle, TwoDimInterp method, Color fill) {
    uint8_t* new_data = new uint8_t[w * h * size];

    double x_old, y_old;
    Color ret(0, 0, 0);
    Interpolation I;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            x_old = (j - origin_x) * cos(-angle * M_PI / 180) - (i - origin_y) * sin(-angle * M_PI / 180) + origin_x;
            y_old = (j - origin_x) * sin(-angle * M_PI / 180) + (i - origin_y) * cos(-angle * M_PI / 180) + origin_y;
            if (method == TwoDimInterp::Nearest) {
                if (x_old < 0 || x_old >= w || y_old < 0 || y_old >= h) {
                    for (int cd = 0; cd < fmin(4, channels); cd++) {
                        new_data[(i * w + j) * channels + cd] = fill.get(cd);
                    }
                }
                for (int cd = 0; cd < fmin(4, channels); cd++) {
                    new_data[(i * w + j) * channels + cd] =
                        data[((uint32_t)y_old * w + (uint32_t)x_old) * channels + cd];
                }
            } else if (method == TwoDimInterp::Bilinear) {
                if (round(x_old) <= -1 || round(x_old) > w || round(y_old) <= -1 || round(y_old) > h) {
                    for (int cd = 0; cd < fmin(4, channels); cd++) {
                        new_data[(i * w + j) * channels + cd] = fill.get(cd);
                    }
                } else {
                    ret = I.bilinear(*this, y_old, x_old, true);
                    for (int cd = 0; cd < fmin(4, channels); cd++) {
                        new_data[(i * w + j) * channels + cd] = ret.get(cd);
                    }
                }
            } else {
                printf("The interpolation method is not yet supported");
            }
        }
    }
    delete[] data;
    data = new_data;
    return *this;
}
