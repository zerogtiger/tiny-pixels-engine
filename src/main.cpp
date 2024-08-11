#include "Image.h"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>

void make_random_complex_arr(uint32_t len, std::complex<double>* z) {
    for (uint32_t i = 0; i < len; i++) {
        z[i] = std::complex<double>(rand() % 100, rand() % 100);
    }
}

void print_complex_arr(uint32_t len, std::complex<double> z[]) {
    for (uint32_t i = 0; i < len; i++) {
        printf("%f + %fi, \n", z[i].real(), z[i].imag());
    }
    printf("\n");
}

void matrix_scalar(double* arr, int size, double scalar) {
    for (int i = 0; i < size; i++) {
        arr[i] /= scalar;
    }
}

int main(int argc, char** argv) {

    Image colorful("images/colorful.jpg");
    colorful.shade_h();
    colorful.write("images/colorful_shade_h.png");
}

void test1() {

    Image test("images/sobel_test.jpg");

    // grayscale
    test.grayscale_avg();
    int test_size = test.w * test.h;

    Image gray_test(test.w, test.h, 1);
    for (uint64_t k = 0; k < test_size; k++) {
        gray_test.data[k] = test.data[k * test.channels];
    }
    gray_test.write("images/test_gray.png");

    double gaussian_blur[] = {1 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 4 / 16.0, 2 / 16.0, 1 / 16.0, 2 / 16.0, 1 / 16.0};

    gray_test.convolve_linear(0, 3, 3, gaussian_blur, 1, 1);
    gray_test.write("images/test_blur.png");
    Image blur_test(gray_test);

    Image blur_test_x(blur_test);
    Image blur_test_y(blur_test);

    double scharr_x[] = {47, 0, -47, 162, 0, -162, 47, 0, -47};
    double scharr_y[] = {47, 162, 47, 0, 0, 0, -47, -162, -47};
    blur_test_x.fd_convolve_clamp_to_border(0, 3, 3, scharr_x, 1, 1, false);
    blur_test_y.fd_convolve_clamp_to_border(0, 3, 3, scharr_y, 1, 1, true);
    blur_test_x.write("images/blur_test_x.png");
    blur_test_y.write("images/blur_test_y.png");

    // make test image
    double mxx = -INFINITY, mxy = -INFINITY, mnx = INFINITY, mny = INFINITY;

    // gets rid of insignificant edges
    double threshold = 0.09;
    double* g = new double[test_size];     // percent
    double* theta = new double[test_size]; // direction
    double x, y;
    for (uint64_t k = 0; k < test_size; k++) {
        x = blur_test_x.data[k*blur_test_x.channels];
        y = blur_test_y.data[k*blur_test_y.channels];
        g[k] = sqrt(x * x/255/255 + y * y/255/255);
        theta[k] = atan2(y, x);
    }

    // makes image hsl -> rgb
    double mx = -INFINITY, mn = INFINITY;

    for (uint64_t k = 0; k < test_size; k++) {
        mx = fmax(mx, g[k]);
        mn = fmin(mn, g[k]);
    }
    Image G(test.w, test.h, 1);
    Image GT(test.w, test.h, 3);

    double h, s, l;
    double v;
    for (uint64_t k = 0; k < test_size; k++) {
        h = theta[k] * 180.0 / M_PI + 180;

        if (mx == mn) {
            v = 0;
        } else {
            v = (g[k] - mn) / (mx - mn) > threshold ? (g[k] - mn) / (mx - mn) : 0;
        }
        s = l = v;

        // hsl -> rgb
        double c = (1 - abs(2 * l - 1)) * s;
        double x = c * (1 - abs(fmod((h / 60), 2) - 1));
        double m = l - c / 2;

        double rt, gt, bt;
        rt = gt = bt = 0;
        if (h < 60) {
            rt = c;
            gt = x;
        } else if (h < 120) {
            rt = x;
            gt = c;
        } else if (h < 180) {
            gt = c;
            bt = x;
        } else if (h < 240) {
            gt = x;
            bt = c;
        } else if (h < 300) {
            bt = c;
            rt = x;
        } else {
            bt = x;
            rt = c;
        }

        uint8_t red = (uint8_t)(255 * (rt + m)), green = (uint8_t)(255 * (gt + m)), blue = (uint8_t)(255 * (bt + m));

        GT.data[k * 3] = red;
        GT.data[k * 3 + 1] = green;
        GT.data[k * 3 + 2] = blue;
        G.data[k] = (uint8_t)(255 * v);
    }

    G.write("images/G.png");
    GT.write("images/GT.png");
    delete[] g;
    delete[] theta;

    //
    //     const uint32_t len = 16;
    //     std::complex<double>* a = new std::complex<double>[len];
    //     make_random_complex_arr(len, a);
    //     print_complex_arr(len, a);
    //
    //     std::complex<double>* A = new std::complex<double>[len];
    //     std::complex<double>* a_recovered = new std::complex<double>[len];
    //
    //     auto fft_2D_start = std::chrono::system_clock::now();
    //     Image::dft_2D(4, 4, a, A);
    //     auto fft_2D_end = std::chrono::system_clock::now();
    //
    //     print_complex_arr(len, A);
    //
    //     auto ifft_2D_start = std::chrono::system_clock::now();
    //     Image::idft_2D(4, 4, A, a_recovered);
    //     auto ifft_2D_end = std::chrono::system_clock::now();
    //     print_complex_arr(len, a_recovered);
    //
    //     printf("2D fft took %lldns\n", std::chrono::duration_cast<std::chrono::nanoseconds>(fft_2D_end - fft_2D_start).count());
    //     printf("2D ifft took %lldns\n", std::chrono::duration_cast<std::chrono::nanoseconds>(ifft_2D_end - ifft_2D_start).count());
    //
    //     delete[] a;
    //     delete[] A;
    //     delete[] a_recovered;

    // Image test("test1.jpg");

    // test.crop(100, 100, 3000, 4000);
    // test.write("cropped.png");

    // Font iosevka("Iosevka Term Nerd Font Complete.ttf", 100);
    // test.overlay_text("TEST", iosevka, 300, 300, 255, 0, 0, 100); // seg fault when 100, 100 for position
    // test.write("overlay_text.png");
    //
    // Image overlay("clear_bg.png");
    // test.overlay(overlay, 400, 500);
    // test.write("overlayed.png");

    // test.flip_y();
    // test.flip_x();
    // test.write("flipped_x.png");

    // embossing
    double ridge[] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    double edge[] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};
    double sharpen[] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
    double emboss[] = {-2 / 3.0, -1 / 3.0, 0, -1 / 3.0, 1 / 3.0, 1 / 3.0, 0, 1 / 3.0, 2 / 3.0};
    double g_blur[] = {1, 4, 6, 4, 1, 4, 16, 24, 16, 4, 6, 24, 36, 24, 6, 4, 16, 24, 16, 4, 1, 4, 6, 4, 1};
    double unsharp[] = {1, 4, 6, 4, 1, 4, 16, 24, 16, 4, 6, 24, -476, 24, 6, 4, 16, 24, 16, 4, 1, 4, 6, 4, 1};
    matrix_scalar(g_blur, 25, 256);
    matrix_scalar(unsharp, 25, -256.0);
    // double emboss[] = {-2 / 3.0, -1 / 3.0, 0, -1 / 3.0, 1 / 3.0, 1 / 3.0, 0, 1 / 3.0, 2 / 3.0};
    // double gaussian_blur[] = {
    //     1/16.0, 2/16.0, 1/16.0,
    //     2/16.0, 4/16.0, 2/16.0,
    //     1/16.0, 2/16.0, 1/16.0
    // };

    // Image t1("colorful2.jpg");
    // Image t2 = t1;
    // Image t0 = t1;
    // t1.std_convolve_clamp_to_border(0, 5, 5, unsharp, 2, 2);
    // t1.std_convolve_clamp_to_border(1, 5, 5, unsharp, 2, 2);
    // t1.std_convolve_clamp_to_border(2, 5, 5, unsharp, 2, 2);
    // t1.std_convolve_clamp_to_border(0, 3, 3, gaussian_blur, 1, 1);
    // t1.std_convolve_clamp_to_border(1, 3, 3, gaussian_blur, 1, 1);
    // t1.std_convolve_clamp_to_border(2, 3, 3, gaussian_blur, 1, 1);
    // t1.std_convolve_clamp_to_border(0, 5, 5, g_blur, 2, 2);
    // t1.std_convolve_clamp_to_border(1, 5, 5, g_blur, 2, 2);
    // t1.std_convolve_clamp_to_border(2, 5, 5, g_blur, 2, 2);
    // t1.fd_convolve_clamp_to_border(0, 5, 5, g_blur, 2, 2);
    // t1.fd_convolve_clamp_to_border(1, 5, 5, g_blur, 2, 2);
    // t1.fd_convolve_clamp_to_border(2, 5, 5, g_blur, 2, 2);
    // t1.write("dft_g_blur_5x5.png");
    // t2.std_convolve_clamp_to_zero(0, 3, 3, gaussian_blur, 1, 1);
    // t2.std_convolve_clamp_to_zero(1, 3, 3, gaussian_blur, 1, 1);
    // t2.std_convolve_clamp_to_zero(2, 3, 3, gaussian_blur, 1, 1);
    //
    // t1.diffmap_scale(t0);
    // t1.write("con.png");

    // Image orig("test1.jpg");
    //
    // test1.diffmap_scale(orig);
    // test1.write("cov_diff.png");

    // Image test("test.jpg");
    // Image test2("test1.jpg");
    // Image out2("out2.png");
    // Image sec("secret.png");
    // test1.diffmap_scale(sec);
    // test1.write("secret_dif_s.png");

    // test.encodemessage("Secret message");
    //
    // test.write("secret.png");
    //
    // Image secret("secret.png");
    //
    // char buffer[256] = {0};
    // size_t len = 0;
    // secret.decodemessage(buffer, &len);
    //
    // printf("Message: %s (%zu)\n", buffer, len);

    // test.color_mask(1, 0.8, 0.9);
    // test.write("out2.png");

    // test.grayscale_lum();
    // test.write("out1.png");

    // test.write("new.png");
    // Image copy = test;
    // copy.write("out.png");
    // Image blank(100, 100, 3);
    // for (int i = 0; i < blank.w * blank.channels; ++i) {
    //     blank.data[i] = 255;
    // }
    // blank.write("blank.jpg");
}
