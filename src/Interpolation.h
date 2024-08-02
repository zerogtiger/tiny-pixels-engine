#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "Image.h"
#include <algorithm>
class Interpolation {
  private:
  public:
    static Color& bilinear(Image image, double r, double c, bool interp_edge_with_fill, Color fill = Color(0, 0, 0)) {
        Color* ret = new Color(0, 0, 0);
        if (round(r) < 0 || round(r) > image.h || round(c) < 0 || round(c) > image.w) {
            ret->r = -1;
            ret->g = -1;
            ret->b = -1;
            ret->a = -1;
            return *ret;
        }

        if (!interp_edge_with_fill) {
            r = std::clamp(r, 0.0, image.h - 1.0);
            c = std::clamp(c, 0.0, image.w - 1.0);
        }
        if (r == floor(r) && c == floor(c)) {
            ret->set(image.get_color_or_default((uint32_t)floor(r), (uint32_t)floor(c), fill));
        } else if (r == floor(r)) {
            ret->set(image.get_color_or_default((uint32_t)r, (uint32_t)floor(c), fill) * (ceil(c) - c) +
                     image.get_color_or_default((uint32_t)r, (uint32_t)ceil(c), fill) * (c - floor(c)));
        } else if (c == floor(c)) {
            ret->set(image.get_color_or_default((uint32_t)floor(r), (uint32_t)c, fill) * (ceil(r) - r) +
                     image.get_color_or_default((uint32_t)ceil(r), (uint32_t)c, fill) * (r - floor(r)));
        } else {
            ret->set(image.get_color_or_default((uint32_t)floor(r), (uint32_t)floor(c), fill) * ((double)ceil(c) - c) *
                         (ceil(r) - r) +
                     image.get_color_or_default((uint32_t)ceil(r), (uint32_t)floor(c), fill) * ((double)ceil(c) - c) *
                         (r - floor(r)) +
                     image.get_color_or_default((uint32_t)floor(r), (uint32_t)ceil(c), fill) * ((double)r - floor(r)) *
                         (ceil(r) - r) +
                     image.get_color_or_default((uint32_t)ceil(r), (uint32_t)ceil(c), fill) * ((double)c - floor(c)) *
                         (r - floor(r)));
        }
        return *ret;
    }
};

#endif
