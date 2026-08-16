#include <ofeli/segmentation.hpp>

#include <ofeli/ac_withoutedges.hpp>
#include <ofeli/ac_withoutedges_yuv.hpp>
#include <ofeli/activecontour.hpp>
#include <ofeli/filters.hpp>
#include <ofeli/geodesic_ac.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{

using ofeli::PreprocessingOptions;

std::size_t checked_pixel_count(const ofeli::Image& image)
{
    if(image.width <= 0 || image.height <= 0)
    {
        throw std::invalid_argument("image dimensions must be positive");
    }
    if(image.channels != 1 && image.channels != 3)
    {
        throw std::invalid_argument("image channels must be 1 or 3");
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
    if(bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("image is too large");
    }
    if(bytes != image.pixels.size())
    {
        throw std::invalid_argument("image pixel buffer size does not match dimensions");
    }
    return pixels;
}

void validate_kernel(int kernel, const char* name)
{
    if(kernel <= 0 || kernel % 2 == 0 || kernel > 499)
    {
        throw std::invalid_argument(std::string(name) + " must be positive, odd, and at most 499");
    }
}

void validate_finite_nonnegative(double value, const char* name)
{
    if(!std::isfinite(value) || value < 0.0)
    {
        throw std::invalid_argument(std::string(name) + " must be finite and non-negative");
    }
}

void validate_preprocessing(const PreprocessingOptions& preprocessing)
{
    validate_finite_nonnegative(preprocessing.gaussian_noise_sigma, "gaussian noise");
    if(!std::isfinite(preprocessing.salt_pepper_probability)
       || preprocessing.salt_pepper_probability < 0.0
       || preprocessing.salt_pepper_probability > 1.0)
    {
        throw std::invalid_argument("salt-and-pepper probability must be in [0, 1]");
    }
    validate_finite_nonnegative(preprocessing.speckle_sigma, "speckle noise");

    validate_kernel(preprocessing.mean_kernel, "mean kernel");
    validate_kernel(preprocessing.gaussian_kernel, "Gaussian kernel");
    validate_finite_nonnegative(preprocessing.gaussian_sigma, "Gaussian sigma");
    validate_kernel(preprocessing.median_kernel, "median kernel");
    if(preprocessing.anisotropic_iterations < 0)
    {
        throw std::invalid_argument("anisotropic iterations must be non-negative");
    }
    validate_finite_nonnegative(preprocessing.anisotropic_lambda, "anisotropic lambda");
    if(preprocessing.anisotropic_lambda > 0.25)
    {
        throw std::invalid_argument("anisotropic lambda must be in [0, 0.25]");
    }
    validate_finite_nonnegative(preprocessing.anisotropic_kappa, "anisotropic kappa");
    if(preprocessing.anisotropic_option != 1 && preprocessing.anisotropic_option != 2)
    {
        throw std::invalid_argument("anisotropic option must be 1 or 2");
    }
    validate_kernel(preprocessing.opening_kernel, "opening kernel");
    validate_kernel(preprocessing.closing_kernel, "closing kernel");
    validate_kernel(preprocessing.top_hat_kernel, "top-hat kernel");
}

void apply_preprocessing(ofeli::Filters& filters, const PreprocessingOptions& preprocessing)
{
    if(!preprocessing.enabled)
    {
        return;
    }

    if(preprocessing.gaussian_noise)
    {
        filters.gaussian_white_noise(preprocessing.gaussian_noise_sigma);
    }
    if(preprocessing.salt_pepper_noise)
    {
        filters.salt_pepper_noise(preprocessing.salt_pepper_probability);
    }
    if(preprocessing.speckle_noise)
    {
        filters.speckle(preprocessing.speckle_sigma);
    }
    if(preprocessing.mean)
    {
        filters.mean_filtering(preprocessing.mean_kernel);
    }
    if(preprocessing.gaussian)
    {
        filters.gaussian_filtering(preprocessing.gaussian_kernel, preprocessing.gaussian_sigma);
    }
    if(preprocessing.median)
    {
        if(preprocessing.optimized_median)
        {
            filters.median_filtering_o1(preprocessing.median_kernel);
        }
        else
        {
            filters.median_filtering_oNlogN(preprocessing.median_kernel);
        }
    }
    if(preprocessing.anisotropic)
    {
        filters.anisotropic_diffusion(preprocessing.anisotropic_iterations,
                                      preprocessing.anisotropic_lambda,
                                      preprocessing.anisotropic_kappa,
                                      preprocessing.anisotropic_option);
    }
    if(preprocessing.opening)
    {
        if(preprocessing.optimized_morphology)
        {
            filters.opening_o1(preprocessing.opening_kernel);
        }
        else
        {
            filters.opening(preprocessing.opening_kernel);
        }
    }
    if(preprocessing.closing)
    {
        if(preprocessing.optimized_morphology)
        {
            filters.closing_o1(preprocessing.closing_kernel);
        }
        else
        {
            filters.closing(preprocessing.closing_kernel);
        }
    }
    if(preprocessing.top_hat)
    {
        if(preprocessing.white_top_hat)
        {
            if(preprocessing.optimized_morphology)
            {
                filters.white_top_hat_o1(preprocessing.top_hat_kernel);
            }
            else
            {
                filters.white_top_hat(preprocessing.top_hat_kernel);
            }
        }
        else if(preprocessing.optimized_morphology)
        {
            filters.black_top_hat_o1(preprocessing.top_hat_kernel);
        }
        else
        {
            filters.black_top_hat(preprocessing.top_hat_kernel);
        }
    }
}

void validate_options(const ofeli::SegmentationOptions& options, std::size_t pixel_count)
{
    if(options.model != ofeli::SegmentationModel::ChanVese
       && options.model != ofeli::SegmentationModel::Geodesic)
    {
        throw std::invalid_argument("unknown segmentation model");
    }
    const ofeli::Initialization& initialization = options.initialization;
    if(initialization.shape != ofeli::InitializationShape::Ellipse
       && initialization.shape != ofeli::InitializationShape::Rectangle)
    {
        throw std::invalid_argument("unknown initialization shape");
    }
    if(!std::isfinite(initialization.width) || initialization.width <= 0.0
       || !std::isfinite(initialization.height) || initialization.height <= 0.0)
    {
        throw std::invalid_argument("initialization shape dimensions must be finite and positive");
    }
    if(!std::isfinite(initialization.center_x) || !std::isfinite(initialization.center_y))
    {
        throw std::invalid_argument("initialization shape center must be finite");
    }
    if(initialization.mask && initialization.mask->size() != pixel_count)
    {
        throw std::invalid_argument("initialization mask size does not match the image");
    }

    const ofeli::SmoothingOptions& smoothing = options.smoothing;
    if(smoothing.iterations < 1 || smoothing.iterations > 999)
    {
        throw std::invalid_argument("smoothing iterations must be in 1..999");
    }
    validate_kernel(smoothing.kernel, "smoothing kernel");
    validate_finite_nonnegative(smoothing.sigma, "smoothing sigma");

    if(options.chan_vese.external_iterations < 1 || options.chan_vese.external_iterations > 999)
    {
        throw std::invalid_argument("external iterations must be in 1..999");
    }
    if(options.chan_vese.outside_weight < 0 || options.chan_vese.inside_weight < 0)
    {
        throw std::invalid_argument("Chan-Vese weights must be non-negative");
    }
    if(options.yuv.y_weight < 0 || options.yuv.u_weight < 0 || options.yuv.v_weight < 0)
    {
        throw std::invalid_argument("YUV weights must be non-negative");
    }
    validate_kernel(options.gradient.kernel, "gradient kernel");
    validate_preprocessing(options.preprocessing);

}

void copy_mask(const char* phi, std::size_t pixel_count, std::vector<unsigned char>& mask)
{
    mask.resize(pixel_count);
    for(std::size_t offset = 0; offset < pixel_count; ++offset)
    {
        mask[offset] = phi[offset] < 0 ? 255 : 0;
    }
}

void copy_offsets(const ofeli::list<int>& source, std::vector<int>& destination)
{
    destination.clear();
    for(ofeli::list<int>::const_iterator it = source.begin(); !it.end(); ++it)
    {
        destination.push_back(*it);
    }
}

}

namespace ofeli
{

SegmentationSession::SegmentationSession(const Image& image, const SegmentationOptions& options)
    : input_(image), options_(options)
{
    const std::size_t pixel_count = checked_pixel_count(input_);
    validate_options(options_, pixel_count);

    try
    {
        filters_ = std::make_unique<Filters>(input_.pixels.data(), input_.width, input_.height,
                                             input_.channels);
        apply_preprocessing(*filters_, options_.preprocessing);

        const unsigned char* contour_image = nullptr;
        if(options_.model == SegmentationModel::Geodesic)
        {
            if(input_.channels == 3)
            {
                filters_->morphological_gradient_yuv(options_.gradient.kernel,
                                                     options_.yuv.y_weight,
                                                     options_.yuv.u_weight,
                                                     options_.yuv.v_weight);
            }
            else
            {
                filters_->morphological_gradient(options_.gradient.kernel);
            }
            contour_image = filters_->get_gradient();
            if(input_.channels == 1)
            {
                contour_image = filters_->get_filtered();
            }
            preprocessed_.width = input_.width;
            preprocessed_.height = input_.height;
            preprocessed_.channels = 1;
            preprocessed_.pixels.assign(contour_image, contour_image + pixel_count);
        }
        else
        {
            contour_image = filters_->get_filtered();
            preprocessed_ = input_;
            preprocessed_.pixels.assign(contour_image, contour_image + input_.pixels.size());
        }

        if(options_.initialization.mask)
        {
            initial_phi_.resize(pixel_count);
            for(std::size_t offset = 0; offset < pixel_count; ++offset)
            {
                initial_phi_[offset] = (*options_.initialization.mask)[offset] == 0 ? 1 : -1;
            }
        }

        const bool ellipse = options_.initialization.shape == InitializationShape::Ellipse;
        if(options_.model == SegmentationModel::Geodesic)
        {
            if(options_.initialization.mask)
            {
                contour_ = std::make_unique<GeodesicAC>(
                    preprocessed_.pixels.data(), input_.width, input_.height, initial_phi_.data(),
                    options_.smoothing.enabled, options_.smoothing.kernel, options_.smoothing.sigma,
                    options_.chan_vese.external_iterations, options_.smoothing.iterations);
            }
            else
            {
                contour_ = std::make_unique<GeodesicAC>(
                    preprocessed_.pixels.data(), input_.width, input_.height, ellipse,
                    options_.initialization.width, options_.initialization.height,
                    options_.initialization.center_x, options_.initialization.center_y,
                    options_.smoothing.enabled, options_.smoothing.kernel, options_.smoothing.sigma,
                    options_.chan_vese.external_iterations, options_.smoothing.iterations);
            }
        }
        else if(input_.channels == 3)
        {
            if(options_.initialization.mask)
            {
                contour_ = std::make_unique<ACwithoutEdgesYUV>(
                    preprocessed_.pixels.data(), input_.width, input_.height, initial_phi_.data(),
                    options_.smoothing.enabled, options_.smoothing.kernel, options_.smoothing.sigma,
                    options_.chan_vese.external_iterations, options_.smoothing.iterations,
                    options_.chan_vese.outside_weight, options_.chan_vese.inside_weight,
                    options_.yuv.y_weight, options_.yuv.u_weight, options_.yuv.v_weight);
            }
            else
            {
                contour_ = std::make_unique<ACwithoutEdgesYUV>(
                    preprocessed_.pixels.data(), input_.width, input_.height, ellipse,
                    options_.initialization.width, options_.initialization.height,
                    options_.initialization.center_x, options_.initialization.center_y,
                    options_.smoothing.enabled, options_.smoothing.kernel, options_.smoothing.sigma,
                    options_.chan_vese.external_iterations, options_.smoothing.iterations,
                    options_.chan_vese.outside_weight, options_.chan_vese.inside_weight,
                    options_.yuv.y_weight, options_.yuv.u_weight, options_.yuv.v_weight);
            }
        }
        else if(options_.initialization.mask)
        {
            contour_ = std::make_unique<ACwithoutEdges>(
                preprocessed_.pixels.data(), input_.width, input_.height, initial_phi_.data(),
                options_.smoothing.enabled, options_.smoothing.kernel, options_.smoothing.sigma,
                options_.chan_vese.external_iterations, options_.smoothing.iterations,
                options_.chan_vese.outside_weight, options_.chan_vese.inside_weight);
        }
        else
        {
            contour_ = std::make_unique<ACwithoutEdges>(
                preprocessed_.pixels.data(), input_.width, input_.height, ellipse,
                options_.initialization.width, options_.initialization.height,
                options_.initialization.center_x, options_.initialization.center_y,
                options_.smoothing.enabled, options_.smoothing.kernel, options_.smoothing.sigma,
                options_.chan_vese.external_iterations, options_.smoothing.iterations,
                options_.chan_vese.outside_weight, options_.chan_vese.inside_weight);
        }
    }
    catch(const std::invalid_argument&)
    {
        throw;
    }
    catch(const std::exception& error)
    {
        throw std::runtime_error(std::string("segmentation setup failed: ") + error.what());
    }
    if(!contour_)
    {
        throw std::runtime_error("segmentation setup failed");
    }
}

SegmentationSession::~SegmentationSession() = default;

void SegmentationSession::step()
{
    if(!finished())
    {
        ++(*contour_);
    }
}

bool SegmentationSession::finished() const noexcept
{
    return cancelled_ || contour_->get_isStopped()
        || contour_->get_iteration() >= contour_->get_iteration_max();
}

void SegmentationSession::cancel() noexcept
{
    cancelled_ = true;
}

SegmentationResult SegmentationSession::snapshot() const
{
    if(!contour_)
    {
        throw std::runtime_error("segmentation session is not initialized");
    }

    SegmentationResult result;
    result.width = input_.width;
    result.height = input_.height;
    result.preprocessed = preprocessed_;
    result.iteration = contour_->get_iteration();
    result.iteration_max = contour_->get_iteration_max();
    result.has_lists_changes = contour_->get_hasListsChanges();
    result.has_oscillation = contour_->get_hasOscillation();
    result.stop_reason = cancelled_ ? StopReason::Cancelled
                                    : (result.iteration >= result.iteration_max
                                           ? StopReason::IterationLimit
                                           : StopReason::Converged);

    const std::size_t pixel_count = static_cast<std::size_t>(input_.width)
        * static_cast<std::size_t>(input_.height);
    result.phi.assign(contour_->get_phi(), contour_->get_phi() + pixel_count);
    copy_mask(contour_->get_phi(), pixel_count, result.mask);
    copy_offsets(contour_->get_Lout(), result.Lout);
    copy_offsets(contour_->get_Lin(), result.Lin);

    if(options_.model == SegmentationModel::Geodesic)
    {
        const GeodesicAC* geodesic = dynamic_cast<const GeodesicAC*>(contour_.get());
        if(!geodesic)
        {
            throw std::runtime_error("segmentation model state is invalid");
        }
        result.otsu_threshold = geodesic->get_otsu_threshold();
    }
    else if(input_.channels == 3)
    {
        const ACwithoutEdgesYUV* rgb = dynamic_cast<const ACwithoutEdgesYUV*>(contour_.get());
        if(!rgb)
        {
            throw std::runtime_error("segmentation model state is invalid");
        }
        result.outside_rgb_mean = std::array<int, 3>{rgb->get_CoutR(), rgb->get_CoutG(), rgb->get_CoutB()};
        result.inside_rgb_mean = std::array<int, 3>{rgb->get_CinR(), rgb->get_CinG(), rgb->get_CinB()};
    }
    else
    {
        const ACwithoutEdges* grayscale = dynamic_cast<const ACwithoutEdges*>(contour_.get());
        if(!grayscale)
        {
            throw std::runtime_error("segmentation model state is invalid");
        }
        result.outside_mean = grayscale->get_Cout();
        result.inside_mean = grayscale->get_Cin();
    }

    return result;
}

SegmentationResult SegmentationSession::run()
{
    while(!finished())
    {
        step();
    }
    return snapshot();
}

}
