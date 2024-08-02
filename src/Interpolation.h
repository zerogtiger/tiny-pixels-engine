#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "Image.h"
#include <algorithm>
class Interpolation {
  private:
  public:
    static Color bilinear(Image& image, double r, double c, bool interp_edge_with_fill, Color fill = Color(0, 0, 0)) {
        Color ret(0, 0, 0);
        if (round(r) < 0 || round(r) > image.h || round(c) < 0 || round(c) > image.w) {
            ret.r = -1;
            ret.g = -1;
            ret.b = -1;
            ret.a = -1;
            return ret;
        }

        if (!interp_edge_with_fill) {
            r = std::clamp(r, 0.0, image.h - 1.0);
            c = std::clamp(c, 0.0, image.w - 1.0);
        }
        int fr = floor(r), cr = ceil(r), fc = floor(c), cc = ceil(c);
        if (r == fr && c == fc) {
            ret.set(image.get_color_or_default(fr, fc, fill));
        } else if (r == fr) {
            ret.set(image.get_color_or_default(r, fc, fill) * (cc - c) +
                     image.get_color_or_default(r, c, fill) * (c - fc));
        } else if (c == fc) {
            ret.set(image.get_color_or_default(fr, c, fill) * (cr - r) +
                     image.get_color_or_default(cr, c, fill) * (r - fr));
        } else {
            ret.set(image.get_color_or_default(fr, fc, fill) * (cc - c) * (cr - r) +
                     image.get_color_or_default(cr, fc, fill) * (cc - c) * (r - fr) +
                     image.get_color_or_default(fr, cc, fill) * (c - fc) * (cr - r) +
                     image.get_color_or_default(cr, cc, fill) * (c - fc) * (r - fr));
        }
        return ret;
    }
};

#endif
