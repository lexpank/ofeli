#include <ofeli/segmentation.hpp>
#include <utility>

int main()
{
    ofeli::Image image{4, 4, 1,
                       {0, 32, 64, 96,
                        32, 64, 96, 128,
                        64, 96, 128, 160,
                        96, 128, 160, 192}};
    ofeli::SegmentationSession session(std::move(image), ofeli::SegmentationOptions{});
    const ofeli::SegmentationResult result = session.run();
    return result.width == 4 && result.height == 4 && !result.mask.empty() ? 0 : 1;
}
