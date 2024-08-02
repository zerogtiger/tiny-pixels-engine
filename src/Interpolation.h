#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "Image.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
class Interpolation {
  private:
    static double eval_cubic_bezier(double t, double p0, double p1, double p2, double p3) {
        return pow(1 - t, 3) * p0 + 3 * pow(1 - t, 2) * t * p1 + 3 * (1 - t) * pow(t, 2) * p2 + pow(t, 3) * p3;
    }

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

    // Notes: assuming control_points is sorted in the order it wish to be interpolated
    static std::vector<double> cubic_bezier(std::vector<std::pair<double, double>> control_points,
                                            std::vector<double> interp_points, double error_bound = 0.0001) {
        std::vector<double> ret;
        double l = 0, r = 1, mid, rslt;
        for (double x : interp_points) {
            l = 0;
            r = 1;
            while (abs(l - r) >= error_bound) {
                mid = (l + r) / 2.0;
                rslt = eval_cubic_bezier(mid, control_points[0].first, control_points[1].first, control_points[2].first,
                                         control_points[3].first);
                if (abs(x - rslt) < error_bound) {
                    break;
                } else if (rslt < x) {
                    l = mid;
                } else {
                    r = mid;
                }
            }
            mid = (l + r) / 2.0;
            ret.push_back(eval_cubic_bezier(mid, control_points[0].second, control_points[1].second,
                                            control_points[2].second, control_points[3].second));
        }
        return ret;
    }
};

#endif
