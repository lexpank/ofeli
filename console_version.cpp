#include <iostream>
#include <QImage>
#include <QFileInfo>

#include "ac_withoutedges.hpp"
#include "ac_withoutedges_yuv.hpp"
#include "geodesic_ac.hpp"
#include "filters.hpp"

static const int list_end = -9999999;

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "Please check arguments!\n";
        exit(1);
    }

    string mode = argv[1];
    string filename_in = argv[2];
    string filename_out = argv[3];

    QImage img = QImage(filename_in.c_str());
    QImage image;

    int img_format = img.format();
    int img_width = img.width();
    int img_height = img.height();
    int img_size = img_width * img_height;

    int positionX = 0, positionY = 0, x, y;

    unsigned char* img1 = NULL;
    bool isRgb1;

    if (img.isGrayscale()) {
        isRgb1 = false;
        img1 = new unsigned char[img_size];

        QRgb pix;

        for( int offset = 0; offset < img_size; offset++ ) {
            // offset = x+y*img_width <==>
            y = offset/img_width;
            x = offset-y*img_width;

            pix = img.pixel(x,y);

            img1[offset] = (unsigned char)(qRed(pix));
        }
    } else {
        isRgb1 = true;
        img1 = new unsigned char[3*img_size];

        QRgb pix;

        for( int offset = 0; offset < img_size; offset++ ) {
            // offset = x+y*img_width <==>
            y = offset/img_width;
            x = offset-y*img_width;

            pix = img.pixel(x,y);

            img1[3*offset] = (unsigned char)(qRed(pix));
            img1[3*offset+1] = (unsigned char)(qGreen(pix));
            img1[3*offset+2] = (unsigned char)(qBlue(pix));
        }
    }

    unsigned char* image_uchar = new unsigned char[3*img_size];
    unsigned char* image_result_uchar = new unsigned char[3*img_size];
    image = QImage(image_uchar, img_width, img_height, 3*img_width, QImage::Format_RGB888);

    if (isRgb1) {
        std::memcpy(image_uchar,img1,3*img_size);
    } else {
        unsigned char I;

        for( int i = 0; i < img_size; i++ ) {
            I = img1[i];

            image_uchar[3*i] = I;
            image_uchar[3*i+1] = I;
            image_uchar[3*i+2] = I;
        }
    }

    bool hasSmoothingCycle1 = false;

    char* shape = new char[img_size];;
    int* shape_points = new int[2*img_size+1];
    shape_points[0] = list_end;
    int* Lout_shape1 = new int[img_size+1];
    Lout_shape1[0] = list_end;
    int* Lin_shape1 = new int[img_size+1];
    Lin_shape1[0] = list_end;
    int* Lout_2 = new int[img_size+1];
    Lout_2[0] = list_end;
    int* Lin_2 = new int[img_size+1];
    Lin_2[0] = list_end;
    char* phi_init1 = new char[img_size];
    char* phi_init2 = new char[img_size];
    char* phi_init1_clean = new char[img_size];
    char* phi_init2_clean = new char[img_size];

    int lambda_out1 = 1;
    int lambda_in1 = 1;
    int kernel_curve1 = 7;
    double std_curve1 = 2;
    int Na1 = 20, Ns1 = 20;

    ofeli::Filters* filters = new ofeli::Filters(img1, img_width, img_height, isRgb1 ? 3 : 1);

    const unsigned char* img1_filtered = NULL;

    if (mode == "1") {
        img1_filtered = filters->get_filtered();
    } else {
        filters->morphological_gradient(3);
        img1_filtered = filters->get_filtered();
    }

    ofeli::ActiveContour* ac = NULL;

    if (mode == "1") {
        ac = new ofeli::ACwithoutEdges(img1_filtered, img_width, img_height, phi_init1_clean, hasSmoothingCycle1, kernel_curve1, std_curve1, Na1, Ns1, lambda_out1, lambda_in1);
    } else {
        ac = new ofeli::GeodesicAC(img1_filtered, img_width, img_height, phi_init1_clean, hasSmoothingCycle1, kernel_curve1, std_curve1, Na1, Ns1);
    }

    const ofeli::list<int>* Lout1;
    const ofeli::list<int>* Lin1;

    Lout1 = &ac->get_Lout();
    Lin1 = &ac->get_Lin();

    int max, min, offset;
    unsigned char I;

    if( !isRgb1 || mode == "2" ) {
        if( mode == "2") {    // && hasHistoNormaliz )
            max = 0;
            min = 255;

            for( offset = 0; offset < img_size; offset++ )
            {
                if( img1_filtered[offset] > max )
                {
                    max = img1_filtered[offset];
                }
                if( img1_filtered[offset] < min )
                {
                    min = img1_filtered[offset];
                }
            }

            for( offset = 0; offset < img_size; offset++ )
            {
                I = (unsigned char)(255.0*double(img1_filtered[offset]-min)/double(max-min));

                image_result_uchar[3*offset] = I;
                image_result_uchar[3*offset+1] = I;
                image_result_uchar[3*offset+2] = I;
            }
        } else {
            for( offset = 0; offset < img_size; offset++ ) {
                I = img1_filtered[offset];

                image_result_uchar[3*offset] = I;
                image_result_uchar[3*offset+1] = I;
                image_result_uchar[3*offset+2] = I;
            }
        }
    } else {
        std::memcpy(image_result_uchar,img1_filtered,3*img_size);
    }

    ac->evolve();
    Lout1 = &ac->get_Lout();
    Lin1 = &ac->get_Lin();

    /**
     * Segmentation
     **/

    const char* phi = ac->get_phi();

    for (int i = 0; i < img_width*img_height; ++i) {
        offset = i*3;

        image_result_uchar[offset] = phi[i] > 0 ? 0 : 255;
        image_result_uchar[offset+1] = phi[i] > 0 ? 0 : 255;
        image_result_uchar[offset+2] = phi[i] > 0 ? 0 : 255;
    }

    QImage output(image_result_uchar, img_width, img_height, 3*img_width, QImage::Format_RGB888);
    output.save(filename_out.c_str());

    return 0;
}
