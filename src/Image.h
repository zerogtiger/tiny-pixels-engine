#include <cstdio>
#include <stdint.h>

enum ImageType { PNG, JPG, BMP, TGA };

struct Image {
    uint8_t *data = NULL; // Butes of image
    size_t size = 0;      // size of data
    int w;
    int h;
    int channels;

    Image(const char *filename);
    Image(int w, int h, int channels);
    Image(const Image &img);
    ~Image(); // destructor

    bool read(const char *filename);
    bool write(const char *filename);

    ImageType getFileType(const char *filename);

    Image &grayscale_avg();
    Image &grayscale_lum();

    Image &color_mask(float r, float g, float b);

    Image &encodemessage(const char *message);

    Image &decodemessage(char *buffer, size_t *messageLength);

    Image& diffmap(Image& img);
    Image& diffmap_scale(Image& img, uint8_t scl = 0);

    Image& std_convolve_clamp_to_zero(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& std_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image& std_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    Image& flip_x();
    Image& flip_y();
};
