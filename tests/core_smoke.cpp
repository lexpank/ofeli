#include "ac_withoutedges.hpp"
#include "ac_withoutedges_yuv.hpp"
#include "filters.hpp"
#include "geodesic_ac.hpp"
#include "linked_list.hpp"

int main()
{
    ofeli::list<int> values(8);
    values.push_front(3);
    values.push_front(1);
    values.push_front(2);
    values.sort();
    if(values.size() != 3)
    {
        return 1;
    }
    ofeli::list<int>::iterator value = values.begin();
    if(*value != 1)
    {
        return 1;
    }
    ++value;
    if(*value != 2)
    {
        return 1;
    }
    ++value;
    if(*value != 3)
    {
        return 1;
    }

    unsigned char image[25] = {};
    image[12] = 90;

    ofeli::Filters filters(image, 5, 5, 1);
    filters.mean_filtering(3);
    const unsigned char* filtered = filters.get_filtered();
    if(filtered[12] != 10)
    {
        return 1;
    }

    filters.gaussian_white_noise(10.0);
    filters.salt_pepper_noise(0.5);
    filters.speckle(0.5);
    for(int offset = 0; offset < 25; ++offset)
    {
        if(filtered[offset] > 255)
        {
            return 1;
        }
    }

    unsigned char noise_image[64 * 64];
    for(int offset = 0; offset < 64 * 64; ++offset)
    {
        noise_image[offset] = 128;
    }

    ofeli::Filters gaussian_noise(noise_image, 64, 64, 1);
    gaussian_noise.gaussian_white_noise(20.0);
    const unsigned char* gaussian_result = gaussian_noise.get_filtered();
    bool gaussian_changed = false;
    for(int offset = 0; offset < 64 * 64; ++offset)
    {
        if(gaussian_result[offset] != 128)
        {
            gaussian_changed = true;
        }
    }
    if(!gaussian_changed)
    {
        return 1;
    }

    ofeli::Filters salt_pepper_noise(noise_image, 64, 64, 1);
    salt_pepper_noise.salt_pepper_noise(1.0);
    const unsigned char* salt_pepper_result = salt_pepper_noise.get_filtered();
    bool has_black = false;
    bool has_white = false;
    for(int offset = 0; offset < 64 * 64; ++offset)
    {
        if(salt_pepper_result[offset] == 0)
        {
            has_black = true;
        }
        else if(salt_pepper_result[offset] == 255)
        {
            has_white = true;
        }
        else
        {
            return 1;
        }
    }
    if(!has_black || !has_white)
    {
        return 1;
    }

    ofeli::Filters speckle_noise(noise_image, 64, 64, 1);
    speckle_noise.speckle(0.5);
    const unsigned char* speckle_result = speckle_noise.get_filtered();
    bool speckle_changed = false;
    for(int offset = 0; offset < 64 * 64; ++offset)
    {
        if(speckle_result[offset] != 128)
        {
            speckle_changed = true;
        }
    }
    if(!speckle_changed)
    {
        return 1;
    }

    ofeli::ACwithoutEdges contour(image, 5, 5,
                                  true, 0.5, 0.5, 0.0, 0.0,
                                  false, 1, 1.0, 1, 1, 1, 1);
    if(contour.get_Lout().empty() || contour.get_Lin().empty())
    {
        return 1;
    }
    ++contour;
    if(contour.get_iteration() != 1)
    {
        return 1;
    }

    ofeli::GeodesicAC geodesic(image, 5, 5,
                               true, 0.5, 0.5, 0.0, 0.0,
                               false, 1, 1.0, 1, 1);
    if(geodesic.get_Lout().empty() || geodesic.get_Lin().empty())
    {
        return 1;
    }
    ++geodesic;
    if(geodesic.get_iteration() != 1)
    {
        return 1;
    }

    unsigned char rgb_image[3 * 64 * 64];
    for(int offset = 0; offset < 64 * 64; ++offset)
    {
        rgb_image[3 * offset] = static_cast<unsigned char>(offset % 256);
        rgb_image[3 * offset + 1] = 128;
        rgb_image[3 * offset + 2] = static_cast<unsigned char>((255 - offset) % 256);
    }

    ofeli::ACwithoutEdgesYUV rgb_contour(rgb_image, 64, 64,
                                         true, 0.5, 0.5, 0.0, 0.0,
                                         false, 1, 1.0, 1, 1,
                                         1, 1, 1, 10, 10);
    if(rgb_contour.get_Lout().empty() || rgb_contour.get_Lin().empty())
    {
        return 1;
    }
    ++rgb_contour;
    if(rgb_contour.get_iteration() != 1)
    {
        return 1;
    }

    return 0;
}
