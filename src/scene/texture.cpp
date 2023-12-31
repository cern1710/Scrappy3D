
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	//A1T6: sample_bilinear
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	// The pixel with the nearest center is the pixel that contains (x,y)
    int32_t i = std::floor(x - 0.5f);
    int32_t j = std::floor(y - 0.5f);

	// Subpixel offsets of texels
	float s = x - (i + 0.5f);
	float t = y - (j + 0.5f);

    // Bound all retrieved coordinates by image width and length
    int32_t im_width = image.w - 1;
    int32_t im_height = image.h - 1;

	// Sample four nearest texels
    int32_t cx00 = std::min(i, im_width);
    int32_t cy00 = std::min(j, im_height);
    int32_t cx01 = std::min(i, im_width);
    int32_t cy01 = std::min(j + 1, im_height);
    int32_t cx10 = std::min(i + 1, im_width);
    int32_t cy10 = std::min(j, im_height);
    int32_t cx11 = std::min(i + 1, im_width);
    int32_t cy11 = std::min(j + 1, im_height);

    // Interpolate coordinates
    Spectrum f00 = image.at(cx00, cy00);
    Spectrum f01 = image.at(cx01, cy01);
    Spectrum f10 = image.at(cx10, cy10);
    Spectrum f11 = image.at(cx11, cy11);

    // Interpolate along x and y axis
    return (1 - t) * ((1 - s) * f00 + s * f10) + t * ((1 - s) * f01 + s * f11);
}


Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	//A1T6: sample_trilinear
	if (lod <= 0) // Use base level when lod (-inf, 0]
		return sample_bilinear(base, uv);
	if (lod >= levels.size()) // Use best mipmap level (last vec value) when lod > number of levels
		return sample_bilinear(levels.back(), uv);
    if (lod <= 1) // Interpolate between base level and level 0 if lod is (0, 1]
        return sample_bilinear(base, uv) * (1 - lod) + sample_bilinear(levels[0], uv) * lod;

    lod -= 1.0f;  // Adjust lod for accessing mipmap levels; lod 0 has been handled (lod <= 1 check)

    int32_t lod_level = std::floor(lod);
    float lod_fraction = lod - lod_level; // Interpolation fraction

    // Calculate the two mip levels to interpolate between
    lod_level = std::clamp(lod_level, 0, int32_t(levels.size()) - 1);
    int32_t next_lod_level = std::min(lod_level + 1, int32_t(levels.size()) - 1);

	// Sample the two mip levels
    Spectrum color1 = sample_bilinear(levels[lod_level], uv);
    Spectrum color2 = sample_bilinear(levels[next_lod_level], uv);

	// Triliner interpolation
    return color1 * (1 - lod_fraction) + color2 * lod_fraction;
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
	auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
		//dst is half the size of src in each dimension:
		assert(std::max(1u, src.w / 2u) == dst.w);
		assert(std::max(1u, src.h / 2u) == dst.h);

		//A1T6: Sample from a 2x2 block in the src image and average them
		for (uint32_t y = 0; y < dst.h; y++)
			for (uint32_t x = 0; x < dst.w; x++)
				dst.at(x, y) = (src.at(2 * x, 2 * y) + src.at(2 * x + 1, 2 * y) +
								src.at(2 * x, 2 * y + 1) + src.at(2 * x + 1, 2 * y + 1)) * 0.25f;

		//Be aware that the alignment of the samples in dst and src will be different depending on whether the image is even or odd.

	};

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;

}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (sampler == Sampler::nearest) {
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		return sample_bilinear(image, uv);
	} else {
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}
