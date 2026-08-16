#include "qt_image_io.hpp"

#include <QFileInfo>
#include <QTemporaryDir>

#include <stdexcept>

namespace
{

template <typename Function>
bool expects_invalid(Function&& function)
{
    try
    {
        function();
    }
    catch(const std::invalid_argument&)
    {
        return true;
    }
    return false;
}

bool pixel_is(const QImage& image, int x, int y, int r, int g, int b)
{
    const QRgb pixel = image.pixel(x, y);
    return qRed(pixel) == r && qGreen(pixel) == g && qBlue(pixel) == b;
}

}

int main()
{
    QImage padded_gray(3, 2, QImage::Format_Grayscale8);
    for(int y = 0; y < padded_gray.height(); ++y)
    {
        unsigned char* row = padded_gray.scanLine(y);
        row[0] = static_cast<unsigned char>(y + 1);
        row[1] = static_cast<unsigned char>(y + 11);
        row[2] = static_cast<unsigned char>(y + 21);
        for(int x = 3; x < padded_gray.bytesPerLine(); ++x)
        {
            row[x] = 250;
        }
    }
    const ofeli::Image gray = ofeli::qimage_to_image(padded_gray);
    if(gray.channels != 1 || gray.pixels != std::vector<unsigned char>{1, 11, 21, 2, 12, 22})
    {
        return 1;
    }

    QImage rgb(2, 2, QImage::Format_RGB888);
    const unsigned char rgb_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    for(int y = 0; y < rgb.height(); ++y)
    {
        unsigned char* row = rgb.scanLine(y);
        for(int x = 0; x < 6; ++x)
        {
            row[x] = rgb_values[y * 6 + x];
        }
    }
    const ofeli::Image rgb_image = ofeli::qimage_to_image(rgb);
    if(rgb_image.channels != 3 || rgb_image.pixels != std::vector<unsigned char>(rgb_values, rgb_values + 12))
    {
        return 1;
    }

    if(ofeli::qimage_to_image(ofeli::image_to_qimage(gray)).pixels != gray.pixels
       || ofeli::qimage_to_image(ofeli::image_to_qimage(rgb_image)).pixels != rgb_image.pixels)
    {
        return 1;
    }

    ofeli::SegmentationResult result;
    result.width = 2;
    result.height = 2;
    result.mask = {0, 1, 255, 2};
    result.phi = {1, -1, -3, -1};
    const QImage mask = ofeli::binary_mask_image(result);
    if(mask.format() != QImage::Format_Grayscale8
       || mask.constScanLine(0)[0] != 0 || mask.constScanLine(0)[1] != 255
       || mask.constScanLine(1)[0] != 255 || mask.constScanLine(1)[1] != 255)
    {
        return 1;
    }

    ofeli::Image overlay_input;
    overlay_input.width = 2;
    overlay_input.height = 2;
    overlay_input.channels = 1;
    overlay_input.pixels = {10, 20, 30, 40};
    result.Lout = {0};
    result.Lin = {3};
    const QImage overlay = ofeli::overlay_image(overlay_input, result);
    if(overlay.format() != QImage::Format_RGB888
       || !pixel_is(overlay, 0, 0, 0, 0, 255)
       || !pixel_is(overlay, 1, 0, 20, 20, 20)
       || !pixel_is(overlay, 1, 1, 255, 0, 0))
    {
        return 1;
    }

    QTemporaryDir temporary;
    if(!temporary.isValid())
    {
        return 1;
    }
    const QString output = temporary.path() + "/mask.png";
    ofeli::save_image_atomic(mask, output, false);
    if(QFileInfo(output).size() == 0
       || !expects_invalid([&] { ofeli::save_image_atomic(mask, output, false); })
       || !expects_invalid([&] { ofeli::save_image_atomic(mask, temporary.path() + "/mask.unknown", true); }))
    {
        return 1;
    }

    return 0;
}
