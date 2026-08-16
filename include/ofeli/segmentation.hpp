#ifndef OFELI_SEGMENTATION_HPP
#define OFELI_SEGMENTATION_HPP

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace ofeli
{

struct Image
{
    int width = 0;
    int height = 0;
    int channels = 0;
    std::vector<unsigned char> pixels;
};

enum class SegmentationModel
{
    ChanVese,
    Geodesic
};

enum class InitializationShape
{
    Ellipse,
    Rectangle
};

struct Initialization
{
    InitializationShape shape = InitializationShape::Ellipse;
    double width = 0.65;
    double height = 0.65;
    double center_x = 0.0;
    double center_y = 0.0;
    std::optional<std::vector<unsigned char>> mask;
};

struct PreprocessingOptions
{
    bool enabled = false;

    bool gaussian_noise = false;
    double gaussian_noise_sigma = 20.0;
    bool salt_pepper_noise = false;
    double salt_pepper_probability = 0.05;
    bool speckle_noise = false;
    double speckle_sigma = 0.16;

    bool mean = false;
    int mean_kernel = 5;
    bool gaussian = false;
    int gaussian_kernel = 5;
    double gaussian_sigma = 2.0;
    bool median = false;
    int median_kernel = 5;
    bool optimized_median = true;
    bool anisotropic = false;
    int anisotropic_option = 1;
    int anisotropic_iterations = 10;
    double anisotropic_lambda = 1.0 / 7.0;
    double anisotropic_kappa = 30.0;
    bool opening = false;
    int opening_kernel = 5;
    bool closing = false;
    int closing_kernel = 5;
    bool top_hat = false;
    int top_hat_kernel = 5;
    bool optimized_morphology = true;
    bool white_top_hat = true;
};

struct SmoothingOptions
{
    bool enabled = true;
    int iterations = 3;
    int kernel = 5;
    double sigma = 2.0;
};

struct GradientOptions
{
    int kernel = 5;
};

struct ChanVeseOptions
{
    int external_iterations = 30;
    int outside_weight = 1;
    int inside_weight = 1;
};

struct YuvOptions
{
    int y_weight = 1;
    int u_weight = 10;
    int v_weight = 10;
};

struct SegmentationOptions
{
    Initialization initialization;
    PreprocessingOptions preprocessing;
    SegmentationModel model = SegmentationModel::ChanVese;
    SmoothingOptions smoothing;
    GradientOptions gradient;
    ChanVeseOptions chan_vese;
    YuvOptions yuv;
};

enum class StopReason
{
    Converged,
    IterationLimit,
    Cancelled
};

struct SegmentationResult
{
    int width = 0;
    int height = 0;
    std::vector<unsigned char> mask;
    std::vector<char> phi;
    std::vector<int> Lout;
    std::vector<int> Lin;
    int iteration = 0;
    int iteration_max = 0;
    bool has_lists_changes = false;
    bool has_oscillation = false;
    StopReason stop_reason = StopReason::Converged;
    Image preprocessed;
    std::optional<int> outside_mean;
    std::optional<int> inside_mean;
    std::optional<std::array<int, 3>> outside_rgb_mean;
    std::optional<std::array<int, 3>> inside_rgb_mean;
    std::optional<unsigned char> otsu_threshold;
};

class Filters;
class ActiveContour;

class SegmentationSession
{
public:
    explicit SegmentationSession(const Image& image,
                                 const SegmentationOptions& options = {});
    ~SegmentationSession();

    SegmentationSession(const SegmentationSession&) = delete;
    SegmentationSession& operator=(const SegmentationSession&) = delete;
    SegmentationSession(SegmentationSession&&) = delete;
    SegmentationSession& operator=(SegmentationSession&&) = delete;

    void step();
    bool finished() const noexcept;
    void cancel() noexcept;
    SegmentationResult snapshot() const;
    SegmentationResult run();

private:
    Image input_;
    SegmentationOptions options_;
    Image preprocessed_;
    std::vector<char> initial_phi_;
    std::unique_ptr<Filters> filters_;
    std::unique_ptr<ActiveContour> contour_;
    bool cancelled_ = false;
};

}

#endif
