#include <ofeli/segmentation.hpp>
#include <ofeli/filters.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

ofeli::Image grayscale_image()
{
    ofeli::Image image;
    image.width = 8;
    image.height = 8;
    image.channels = 1;
    image.pixels.resize(64);
    for(int offset = 0; offset < 64; ++offset)
    {
        image.pixels[offset] = static_cast<unsigned char>((offset * 37 + 11) % 256);
    }
    return image;
}

ofeli::Image rgb_image()
{
    ofeli::Image image;
    image.width = 8;
    image.height = 8;
    image.channels = 3;
    image.pixels.resize(3 * 64);
    for(int offset = 0; offset < 64; ++offset)
    {
        image.pixels[3 * offset] = static_cast<unsigned char>((offset * 17 + 3) % 256);
        image.pixels[3 * offset + 1] = static_cast<unsigned char>((offset * 29 + 7) % 256);
        image.pixels[3 * offset + 2] = static_cast<unsigned char>((offset * 43 + 13) % 256);
    }
    return image;
}

ofeli::SegmentationOptions quick_options()
{
    ofeli::SegmentationOptions options;
    options.smoothing.enabled = false;
    options.chan_vese.external_iterations = 3;
    options.smoothing.iterations = 1;
    return options;
}

ofeli::Image filtered_image(const ofeli::Image& image,
                            const ofeli::PreprocessingOptions& preprocessing)
{
    ofeli::Filters filters(image.pixels.data(), image.width, image.height, image.channels);
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

    ofeli::Image result = image;
    result.pixels.assign(filters.get_filtered(), filters.get_filtered() + image.pixels.size());
    return result;
}

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

}

int main()
{
    const ofeli::Image gray = grayscale_image();
    const ofeli::Image rgb = rgb_image();
    const ofeli::SegmentationOptions options = quick_options();

    ofeli::SegmentationSession gray_session(gray, options);
    const ofeli::SegmentationResult gray_result = gray_session.run();
    if(gray_result.preprocessed.channels != 1 || gray_result.mask.size() != 64
       || !gray_result.outside_mean || !gray_result.inside_mean
       || gray_result.iteration <= 0 || !gray_session.finished())
    {
        return 1;
    }

    ofeli::SegmentationOptions rgb_options = options;
    rgb_options.model = ofeli::SegmentationModel::ChanVese;
    ofeli::SegmentationSession rgb_session(rgb, rgb_options);
    const ofeli::SegmentationResult rgb_result = rgb_session.run();
    if(rgb_result.preprocessed.channels != 3 || !rgb_result.outside_rgb_mean
       || !rgb_result.inside_rgb_mean)
    {
        return 1;
    }

    ofeli::SegmentationOptions geodesic_options = options;
    geodesic_options.model = ofeli::SegmentationModel::Geodesic;
    ofeli::SegmentationSession geodesic_session(gray, geodesic_options);
    const ofeli::SegmentationResult geodesic_result = geodesic_session.run();
    if(geodesic_result.preprocessed.channels != 1 || !geodesic_result.otsu_threshold)
    {
        return 1;
    }

    ofeli::SegmentationOptions rgb_geodesic_options = geodesic_options;
    ofeli::SegmentationSession rgb_geodesic_session(rgb, rgb_geodesic_options);
    if(rgb_geodesic_session.run().preprocessed.channels != 1)
    {
        return 1;
    }

    ofeli::SegmentationOptions mask_options = options;
    mask_options.initialization.mask = std::vector<unsigned char>(64, 0);
    for(int offset = 0; offset < 64; offset += 3)
    {
        (*mask_options.initialization.mask)[offset] = 7;
    }
    ofeli::SegmentationSession mask_session(gray, mask_options);
    const ofeli::SegmentationResult mask_result = mask_session.snapshot();
    for(int offset = 0; offset < 64; ++offset)
    {
        const unsigned char expected = offset % 3 == 0 ? 255 : 0;
        if(mask_result.mask[offset] != expected)
        {
            return 1;
        }
    }

    ofeli::SegmentationSession cancelled_session(gray, options);
    cancelled_session.cancel();
    if(!cancelled_session.finished()
       || cancelled_session.snapshot().stop_reason != ofeli::StopReason::Cancelled)
    {
        return 1;
    }

    if(!expects_invalid([&] {
        ofeli::Image invalid = gray;
        invalid.width = 0;
        ofeli::SegmentationSession session(invalid, options);
    }) || !expects_invalid([&] {
        ofeli::Image invalid = gray;
        invalid.channels = 2;
        ofeli::SegmentationSession session(invalid, options);
    }) || !expects_invalid([&] {
        ofeli::Image invalid = gray;
        invalid.pixels.pop_back();
        ofeli::SegmentationSession session(invalid, options);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.initialization.mask = std::vector<unsigned char>(63, 1);
        ofeli::SegmentationSession session(gray, invalid);
    }))
    {
        return 1;
    }

    if(!expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.chan_vese.external_iterations = 0;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.smoothing.iterations = 1000;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.smoothing.kernel = 4;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.gradient.kernel = 501;
        ofeli::SegmentationSession session(gray, invalid);
    }))
    {
        return 1;
    }

    if(!expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.initialization.width = 0.0;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.initialization.center_x = std::numeric_limits<double>::quiet_NaN();
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.smoothing.sigma = std::numeric_limits<double>::infinity();
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.chan_vese.outside_weight = -1;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.preprocessing.gaussian_noise_sigma = -1.0;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.preprocessing.salt_pepper_probability = 1.1;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.preprocessing.anisotropic_lambda = 0.3;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.preprocessing.anisotropic_option = 3;
        ofeli::SegmentationSession session(gray, invalid);
    }) || !expects_invalid([&] {
        ofeli::SegmentationOptions invalid = options;
        invalid.preprocessing.anisotropic_kappa = -1.0;
        ofeli::SegmentationSession session(gray, invalid);
    }))
    {
        return 1;
    }

    ofeli::SegmentationOptions preprocess_options = options;
    preprocess_options.preprocessing.enabled = true;
    preprocess_options.preprocessing.mean = true;
    preprocess_options.preprocessing.mean_kernel = 3;
    preprocess_options.preprocessing.gaussian = true;
    preprocess_options.preprocessing.gaussian_kernel = 3;
    preprocess_options.preprocessing.gaussian_sigma = 1.0;
    preprocess_options.preprocessing.median = true;
    preprocess_options.preprocessing.median_kernel = 3;
    preprocess_options.preprocessing.optimized_median = false;
    preprocess_options.preprocessing.anisotropic = true;
    preprocess_options.preprocessing.anisotropic_iterations = 1;
    preprocess_options.preprocessing.anisotropic_lambda = 0.1;
    preprocess_options.preprocessing.anisotropic_kappa = 30.0;
    preprocess_options.preprocessing.anisotropic_option = 2;
    preprocess_options.preprocessing.opening = true;
    preprocess_options.preprocessing.opening_kernel = 3;
    preprocess_options.preprocessing.closing = true;
    preprocess_options.preprocessing.closing_kernel = 3;
    preprocess_options.preprocessing.top_hat = true;
    preprocess_options.preprocessing.top_hat_kernel = 3;
    preprocess_options.preprocessing.optimized_morphology = false;
    preprocess_options.preprocessing.white_top_hat = false;

    const ofeli::Image expected_gray = filtered_image(gray, preprocess_options.preprocessing);
    ofeli::SegmentationSession preprocessed_gray(gray, preprocess_options);
    if(preprocessed_gray.snapshot().preprocessed.pixels != expected_gray.pixels)
    {
        return 1;
    }

    ofeli::SegmentationOptions rgb_preprocess_options = options;
    rgb_preprocess_options.preprocessing.enabled = true;
    rgb_preprocess_options.preprocessing.mean = true;
    rgb_preprocess_options.preprocessing.mean_kernel = 3;
    rgb_preprocess_options.preprocessing.median = true;
    rgb_preprocess_options.preprocessing.median_kernel = 3;
    rgb_preprocess_options.preprocessing.optimized_median = true;
    rgb_preprocess_options.preprocessing.opening = true;
    rgb_preprocess_options.preprocessing.opening_kernel = 3;
    rgb_preprocess_options.preprocessing.closing = true;
    rgb_preprocess_options.preprocessing.closing_kernel = 3;
    rgb_preprocess_options.preprocessing.top_hat = true;
    rgb_preprocess_options.preprocessing.top_hat_kernel = 3;
    rgb_preprocess_options.preprocessing.optimized_morphology = true;
    rgb_preprocess_options.preprocessing.white_top_hat = true;

    const ofeli::Image expected_rgb = filtered_image(rgb, rgb_preprocess_options.preprocessing);
    ofeli::SegmentationSession preprocessed_rgb(rgb, rgb_preprocess_options);
    if(preprocessed_rgb.snapshot().preprocessed.pixels != expected_rgb.pixels)
    {
        return 1;
    }

    return 0;
}
