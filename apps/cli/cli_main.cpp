#include "qt_image_io.hpp"

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QFileInfo>
#include <QImageReader>

#include <chrono>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{

class UsageError : public std::runtime_error
{
public:
    explicit UsageError(const std::string& message) : std::runtime_error(message) {}
};

[[noreturn]] void usage_error(const QString& message)
{
    throw UsageError(message.toStdString());
}

ofeli::Image load_image(const QString& path)
{
    QImageReader reader(path);
    const QImage image = reader.read();
    if(image.isNull())
    {
        const QString detail = reader.errorString().isEmpty()
            ? QStringLiteral("unknown image error") : reader.errorString();
        throw std::runtime_error("could not read image '" + path.toStdString() + "': "
                                 + detail.toStdString());
    }
    try
    {
        return ofeli::qimage_to_image(image);
    }
    catch(const std::exception& error)
    {
        throw std::runtime_error("could not decode image '" + path.toStdString() + "': "
                                 + error.what());
    }
}

void validate_output(const QString& path, bool force)
{
    try
    {
        ofeli::validate_output_path(path, force);
    }
    catch(const std::invalid_argument& error)
    {
        throw UsageError(error.what());
    }
}

void save_output(const QImage& image, const QString& path, bool force)
{
    try
    {
        ofeli::save_image_atomic(image, path, force);
    }
    catch(const std::invalid_argument& error)
    {
        throw UsageError(error.what());
    }
}

const char* model_name(ofeli::SegmentationModel model)
{
    return model == ofeli::SegmentationModel::Geodesic ? "geodesic" : "chan-vese";
}

const char* stop_name(ofeli::StopReason reason)
{
    switch(reason)
    {
    case ofeli::StopReason::Converged:
        return "converged";
    case ofeli::StopReason::IterationLimit:
        return "iteration-limit";
    case ofeli::StopReason::Cancelled:
        return "cancelled";
    }
    return "unknown";
}

int run_cli(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName(QStringLiteral("ofeli-cli"));
#ifdef OFELI_VERSION
    QCoreApplication::setApplicationVersion(QStringLiteral(OFELI_VERSION));
#else
    QCoreApplication::setApplicationVersion(QStringLiteral("unknown"));
#endif

    QCommandLineParser parser;
    parser.setApplicationDescription(QStringLiteral("Batch active-contour segmentation"));
    parser.setSingleDashWordOptionMode(QCommandLineParser::ParseAsLongOptions);
    const QCommandLineOption help_option = parser.addHelpOption();
    const QCommandLineOption version_option = parser.addVersionOption();
    const QCommandLineOption mask_option({QStringLiteral("mask")},
                                          QStringLiteral("write the binary mask"),
                                          QStringLiteral("output"));
    const QCommandLineOption overlay_option({QStringLiteral("overlay")},
                                             QStringLiteral("write an RGB contour overlay"),
                                             QStringLiteral("output"));
    const QCommandLineOption model_option({QStringLiteral("model")},
                                           QStringLiteral("segmentation model (chan-vese or geodesic)"),
                                           QStringLiteral("name"));
    const QCommandLineOption initial_mask_option({QStringLiteral("initial-mask")},
                                                  QStringLiteral("use a grayscale image as the initial mask"),
                                                  QStringLiteral("input"));
    const QCommandLineOption force_option({QStringLiteral("force")},
                                           QStringLiteral("allow output overwrite"));
    parser.addOption(mask_option);
    parser.addOption(overlay_option);
    parser.addOption(model_option);
    parser.addOption(initial_mask_option);
    parser.addOption(force_option);

    if(!parser.parse(app.arguments()))
    {
        usage_error(parser.errorText());
    }
    if(parser.isSet(help_option))
    {
        std::cout << parser.helpText().toStdString();
        return 0;
    }
    if(parser.isSet(version_option))
    {
        std::cout << app.applicationVersion().toStdString() << '\n';
        return 0;
    }

    const QStringList positional = parser.positionalArguments();
    if(positional.size() != 2 || positional[0] != QStringLiteral("segment"))
    {
        usage_error(QStringLiteral("usage: ofeli-cli segment INPUT --mask OUTPUT [options]"));
    }
    if(!parser.isSet(mask_option) || parser.value(mask_option).isEmpty())
    {
        usage_error(QStringLiteral("--mask is required"));
    }

    ofeli::SegmentationModel model = ofeli::SegmentationModel::ChanVese;
    const QString requested_model = parser.isSet(model_option)
        ? parser.value(model_option) : QStringLiteral("chan-vese");
    if(requested_model == QStringLiteral("geodesic"))
    {
        model = ofeli::SegmentationModel::Geodesic;
    }
    else if(requested_model != QStringLiteral("chan-vese"))
    {
        usage_error(QStringLiteral("--model must be chan-vese or geodesic"));
    }

    const bool force = parser.isSet(force_option);
    const QString mask_path = parser.value(mask_option);
    const QString overlay_path = parser.value(overlay_option);
    validate_output(mask_path, force);
    if(parser.isSet(overlay_option))
    {
        if(overlay_path.isEmpty())
        {
            usage_error(QStringLiteral("--overlay requires an output path"));
        }
        if(QFileInfo(mask_path).absoluteFilePath() == QFileInfo(overlay_path).absoluteFilePath())
        {
            usage_error(QStringLiteral("--mask and --overlay must be different files"));
        }
        validate_output(overlay_path, force);
    }

    const auto start = std::chrono::steady_clock::now();
    const ofeli::Image input = load_image(positional[1]);

    ofeli::SegmentationOptions options;
    options.model = model;
    if(parser.isSet(initial_mask_option))
    {
        const ofeli::Image initial_mask = load_image(parser.value(initial_mask_option));
        if(initial_mask.channels != 1 || initial_mask.width != input.width
           || initial_mask.height != input.height)
        {
            usage_error(QStringLiteral("initial mask must be a matching grayscale image"));
        }
        options.initialization.mask = initial_mask.pixels;
    }

    ofeli::SegmentationResult result;
    try
    {
        ofeli::SegmentationSession session(input, options);
        result = session.run();
    }
    catch(const std::invalid_argument& error)
    {
        throw UsageError(error.what());
    }

    save_output(ofeli::binary_mask_image(result), mask_path, force);
    if(parser.isSet(overlay_option))
    {
        save_output(ofeli::overlay_image(input, result), overlay_path, force);
    }

    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << "model=" << model_name(model)
              << " width=" << result.width
              << " height=" << result.height
              << " iterations=" << result.iteration
              << " stop=" << stop_name(result.stop_reason)
              << " elapsed_ms=" << elapsed << '\n';
    return 0;
}

}

int main(int argc, char* argv[])
{
    try
    {
        return run_cli(argc, argv);
    }
    catch(const UsageError& error)
    {
        std::cerr << error.what() << '\n';
        return 2;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
