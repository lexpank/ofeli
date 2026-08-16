#include "qt_image_io.hpp"

#include <QDir>
#include <QFileInfo>
#include <QImageWriter>
#include <QSaveFile>

#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{

std::size_t checked_image_size(const ofeli::Image& image)
{
    if(image.width <= 0 || image.height <= 0 || (image.channels != 1 && image.channels != 3))
    {
        throw std::invalid_argument("image dimensions or channels are invalid");
    }
    const std::size_t width = static_cast<std::size_t>(image.width);
    const std::size_t height = static_cast<std::size_t>(image.height);
    const std::size_t channels = static_cast<std::size_t>(image.channels);
    if(width > std::numeric_limits<std::size_t>::max() / height)
    {
        throw std::invalid_argument("image size overflows");
    }
    const std::size_t pixels = width * height;
    if(pixels > std::numeric_limits<std::size_t>::max() / channels)
    {
        throw std::invalid_argument("image size overflows");
    }
    const std::size_t bytes = pixels * channels;
    if(bytes != image.pixels.size())
    {
        throw std::invalid_argument("image pixel buffer size does not match dimensions");
    }
    return pixels;
}

std::runtime_error image_error(const QString& message)
{
    return std::runtime_error(message.toStdString());
}

void validate_offsets(const ofeli::SegmentationResult& result, std::size_t pixel_count)
{
    for(const std::vector<int>* offsets : {&result.Lout, &result.Lin})
    {
        for(int offset : *offsets)
        {
            if(offset < 0 || static_cast<std::size_t>(offset) >= pixel_count)
            {
                throw std::invalid_argument("contour offset is outside the image");
            }
        }
    }
}

}

namespace ofeli
{

Image qimage_to_image(const QImage& image)
{
    if(image.isNull())
    {
        throw std::invalid_argument("input image is null");
    }

    const bool grayscale = image.isGrayscale();
    const QImage normalized = image.convertToFormat(
        grayscale ? QImage::Format_Grayscale8 : QImage::Format_RGB888);
    if(normalized.isNull())
    {
        throw image_error("could not convert image format");
    }

    Image result;
    result.width = normalized.width();
    result.height = normalized.height();
    result.channels = grayscale ? 1 : 3;
    const std::size_t row_bytes = static_cast<std::size_t>(result.width)
        * static_cast<std::size_t>(result.channels);
    result.pixels.resize(static_cast<std::size_t>(result.height) * row_bytes);
    for(int y = 0; y < result.height; ++y)
    {
        std::memcpy(result.pixels.data() + static_cast<std::size_t>(y) * row_bytes,
                    normalized.constScanLine(y), row_bytes);
    }
    return result;
}

QImage image_to_qimage(const Image& image)
{
    checked_image_size(image);
    const QImage::Format format = image.channels == 1
        ? QImage::Format_Grayscale8 : QImage::Format_RGB888;
    QImage result(image.width, image.height, format);
    if(result.isNull())
    {
        throw image_error("could not allocate image");
    }

    const std::size_t row_bytes = static_cast<std::size_t>(image.width)
        * static_cast<std::size_t>(image.channels);
    for(int y = 0; y < image.height; ++y)
    {
        std::memcpy(result.scanLine(y), image.pixels.data() + static_cast<std::size_t>(y) * row_bytes,
                    row_bytes);
    }
    return result;
}

QImage binary_mask_image(const SegmentationResult& result)
{
    if(result.width <= 0 || result.height <= 0)
    {
        throw std::invalid_argument("mask dimensions are invalid");
    }
    const std::size_t pixel_count = static_cast<std::size_t>(result.width)
        * static_cast<std::size_t>(result.height);
    if(result.mask.size() != pixel_count && result.phi.size() != pixel_count)
    {
        throw std::invalid_argument("mask or level-set size does not match dimensions");
    }

    QImage mask(result.width, result.height, QImage::Format_Grayscale8);
    if(mask.isNull())
    {
        throw image_error("could not allocate mask image");
    }
    for(int y = 0; y < result.height; ++y)
    {
        unsigned char* row = mask.scanLine(y);
        for(int x = 0; x < result.width; ++x)
        {
            const std::size_t offset = static_cast<std::size_t>(y) * result.width + x;
            row[x] = result.phi.size() == pixel_count
                ? (result.phi[offset] < 0 ? 255 : 0)
                : (result.mask[offset] == 0 ? 0 : 255);
        }
    }
    return mask;
}

QImage overlay_image(const Image& image, const SegmentationResult& result)
{
    const std::size_t pixel_count = checked_image_size(image);
    if(result.width != image.width || result.height != image.height)
    {
        throw std::invalid_argument("overlay dimensions do not match the image");
    }
    validate_offsets(result, pixel_count);

    QImage overlay = image_to_qimage(image).convertToFormat(QImage::Format_RGB888);
    for(int offset : result.Lout)
    {
        unsigned char* pixel = overlay.scanLine(offset / image.width) + 3 * (offset % image.width);
        pixel[0] = 0;
        pixel[1] = 0;
        pixel[2] = 255;
    }
    for(int offset : result.Lin)
    {
        unsigned char* pixel = overlay.scanLine(offset / image.width) + 3 * (offset % image.width);
        pixel[0] = 255;
        pixel[1] = 0;
        pixel[2] = 0;
    }
    return overlay;
}

QByteArray image_format_for_path(const QString& path)
{
    const QByteArray format = QFileInfo(path).suffix().toLatin1().toLower();
    if(format.isEmpty())
    {
        throw std::invalid_argument("output filename must have an image suffix");
    }
    for(const QByteArray& supported : QImageWriter::supportedImageFormats())
    {
        if(supported.toLower() == format)
        {
            return format;
        }
    }
    throw std::invalid_argument("unsupported output image suffix: " + format.toStdString());
}

void validate_output_path(const QString& path, bool force)
{
    image_format_for_path(path);
    const QFileInfo info(path);
    if(!info.absoluteDir().exists())
    {
        throw std::invalid_argument("output directory does not exist");
    }
    if(!force && info.exists())
    {
        throw std::invalid_argument("output exists; use --force to overwrite");
    }
}

void save_image_atomic(const QImage& image, const QString& path, bool force)
{
    const QByteArray format = image_format_for_path(path);
    validate_output_path(path, force);

    QSaveFile file(path);
    if(!file.open(QIODevice::WriteOnly))
    {
        throw image_error("could not open output: " + file.errorString());
    }
    if(!image.save(&file, format.constData()))
    {
        file.cancelWriting();
        throw image_error("could not encode output: " + file.errorString());
    }
    if(!file.commit())
    {
        throw image_error("could not commit output: " + file.errorString());
    }
}

}
