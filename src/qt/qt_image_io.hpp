#ifndef OFELI_QT_IMAGE_IO_HPP
#define OFELI_QT_IMAGE_IO_HPP

#include <ofeli/segmentation.hpp>

#include <QByteArray>
#include <QImage>
#include <QString>

namespace ofeli
{

Image qimage_to_image(const QImage& image);
QImage image_to_qimage(const Image& image);
QImage binary_mask_image(const SegmentationResult& result);
QImage overlay_image(const Image& image, const SegmentationResult& result);

QByteArray image_format_for_path(const QString& path);
void validate_output_path(const QString& path, bool force);
void save_image_atomic(const QImage& image, const QString& path, bool force);

}

#endif
