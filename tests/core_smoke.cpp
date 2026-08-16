#include <ofeli/ac_withoutedges.hpp>
#include <ofeli/ac_withoutedges_yuv.hpp>
#include <ofeli/activecontour.hpp>
#include <ofeli/filters.hpp>
#include <ofeli/geodesic_ac.hpp>
#include <ofeli/linked_list.hpp>

#include <cstdint>
#include <iostream>

namespace
{

std::uint64_t phi_hash(const ofeli::ActiveContour& contour, int image_size, int& inside)
{
    constexpr std::uint64_t offset_basis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    const char* phi = contour.get_phi();
    std::uint64_t hash = offset_basis;
    inside = 0;

    for(int offset = 0; offset < image_size; ++offset)
    {
        if(phi[offset] != -3 && phi[offset] != -1 && phi[offset] != 1 && phi[offset] != 3)
        {
            return 0;
        }
        if(phi[offset] < 0)
        {
            ++inside;
        }
        hash ^= static_cast<unsigned char>(phi[offset]);
        hash *= prime;
    }
    return hash;
}

bool check_contour(const char* name, const ofeli::ActiveContour& contour,
                   int width, int height, int expected_inside, std::uint64_t expected_hash)
{
    const int image_size = width * height;
    const int iterations = contour.get_iteration();
    if(iterations <= 0 || iterations > contour.get_iteration_max() || !contour.get_isStopped())
    {
        return false;
    }

    int inside = 0;
    const std::uint64_t hash = phi_hash(contour, image_size, inside);
    if(hash == 0)
    {
        return false;
    }

    for(const auto* list : {&contour.get_Lout(), &contour.get_Lin()})
    {
        for(auto it = list->begin(); !it.end(); ++it)
        {
            if(*it < 0 || *it >= image_size)
            {
                return false;
            }
        }
    }

    if(inside != expected_inside || hash != expected_hash)
    {
        std::cerr << name << " iteration=" << iterations
                  << " inside=" << inside << " hash=" << hash << '\n';
        return false;
    }
    return true;
}

}

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

    constexpr int width = 8;
    constexpr int height = 8;
    unsigned char filter_image[25] = {};
    filter_image[12] = 90;
    ofeli::Filters filters(filter_image, 5, 5, 1);
    filters.mean_filtering(3);
    if(filters.get_filtered()[12] != 10)
    {
        return 1;
    }

    unsigned char grayscale_image[width * height];
    for(int offset = 0; offset < width * height; ++offset)
    {
        grayscale_image[offset] = static_cast<unsigned char>((offset * 37 + 11) % 256);
    }

    unsigned char rgb_image[3 * width * height];
    for(int offset = 0; offset < width * height; ++offset)
    {
        rgb_image[3 * offset] = static_cast<unsigned char>((offset * 17 + 3) % 256);
        rgb_image[3 * offset + 1] = static_cast<unsigned char>((offset * 29 + 7) % 256);
        rgb_image[3 * offset + 2] = static_cast<unsigned char>((offset * 43 + 13) % 256);
    }

    ofeli::ACwithoutEdges contour(grayscale_image, width, height,
                                  true, 0.5, 0.5, 0.0, 0.0,
                                  true, 5, 2.0, 30, 3, 1, 1);
    contour.evolve();
    bool ok = check_contour("chan-vese", contour, width, height,
                            5, 7436660801394953399ULL);

    ofeli::ACwithoutEdgesYUV rgb_contour(rgb_image, width, height,
                                         true, 0.5, 0.5, 0.0, 0.0,
                                         true, 5, 2.0, 30, 3,
                                         1, 1, 1, 10, 10);
    rgb_contour.evolve();
    ok = check_contour("chan-vese-yuv", rgb_contour, width, height,
                       36, 8158441133111810505ULL) && ok;

    ofeli::GeodesicAC geodesic(grayscale_image, width, height,
                               true, 0.5, 0.5, 0.0, 0.0,
                               true, 5, 2.0, 30, 3);
    geodesic.evolve();
    ok = check_contour("geodesic", geodesic, width, height,
                       37, 17170649059186662957ULL) && ok;

    return ok ? 0 : 1;
}
