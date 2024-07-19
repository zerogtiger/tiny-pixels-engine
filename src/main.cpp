#include "Image.h"

int main(int argc, char **argv) {

    Image test("test1.jpg");

    // test.flip_y();
    // test.flip_x();
    // test.write("flipped_x.png");

    //embossing
    double emboss[] = {
        -2/3.0, -1/3.0, 0,
        -1/3.0, 1/3.0, 1/3.0,
        0, 1/3.0, 2/3.0}; 
    // double gaussian_blur[] = {
    //     1/16.0, 2/16.0, 1/16.0,
    //     2/16.0, 4/16.0, 2/16.0, 
    //     1/16.0, 2/16.0, 1/16.0};
    // 
    Image t1("colorful2.jpg");
    // Image t2 = t1;
    // Image t0 = t1;
    t1.std_convolve_clamp_to_border(0, 3, 3, emboss, 1, 1);
    t1.std_convolve_clamp_to_border(1, 3, 3, emboss, 1, 1);
    t1.std_convolve_clamp_to_border(2, 3, 3, emboss, 1, 1);
    // t1.write("convolution.png");
    // t2.std_convolve_clamp_to_zero(0, 3, 3, gaussian_blur, 1, 1);
    // t2.std_convolve_clamp_to_zero(1, 3, 3, gaussian_blur, 1, 1);
    // t2.std_convolve_clamp_to_zero(2, 3, 3, gaussian_blur, 1, 1);
    //
    // t1.diffmap_scale(t0);
    t1.write("con.png");
    

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
