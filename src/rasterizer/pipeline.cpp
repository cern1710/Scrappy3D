#include "pipeline.h"

#include "framebuffer.h"
#include "sample_pattern.h"

#include "../lib/log.h"
#include "../lib/mathlib.h"

#include <iostream>

template< PrimitiveType primitive_type, class Program, uint32_t flags >
void Pipeline< primitive_type, Program, flags >::run(
		std::vector< Vertex > const &vertices,
		typename Program::Parameters const &parameters,
		Framebuffer *framebuffer_) {
	//Framebuffer must be non-null:
	assert(framebuffer_);
	auto &framebuffer = *framebuffer_;

	//A1T7: sample loop
	//TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	// This will probably involve inserting a loop of the form:
	//     std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//     for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//  around some subset of the code.
	// You will also need to transform the input and output of the rasterize_* functions to
	//   account for the fact they deal with pixels centered at (0.5,0.5).

	std::vector< ShadedVertex > shaded_vertices;
	shaded_vertices.reserve(vertices.size());

	//--------------------------
	//shade vertices:
	for (auto const &v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex( parameters, v.attributes, &sv.clip_position, &sv.attributes );
		shaded_vertices.emplace_back(sv);
	}

	//--------------------------
	//assemble + clip + homogeneous divide vertices:
	std::vector< ClippedVertex > clipped_vertices;

	//reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		//clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		//clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	//helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const &sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	//actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line( shaded_vertices[i], shaded_vertices[i+1], emit_vertex );
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle( shaded_vertices[i], shaded_vertices[i+1], shaded_vertices[i+2], emit_vertex );
		}
	} else {
		static_assert( primitive_type == PrimitiveType::Lines, "Unsupported primitive type." );
	}


	//--------------------------
	//rasterize primitives:

	std::vector< Fragment > fragments;

	//helper used to put output of rasterization functions into fragments:
	auto emit_fragment = [&](Fragment const &f) {
		fragments.emplace_back(f);
	};
	//actually do rasterization:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
			rasterize_line( clipped_vertices[i], clipped_vertices[i+1], emit_fragment );
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
			rasterize_triangle( clipped_vertices[i], clipped_vertices[i+1], clipped_vertices[i+2], emit_fragment );
		}
	} else {
		static_assert( primitive_type == PrimitiveType::Lines, "Unsupported primitive type." );
	}

	//--------------------------
	//depth test + shade + blend fragments:
	uint32_t out_of_range = 0; //check if rasterization produced fragments outside framebuffer (indicates something is wrong with clipping)
	for (auto const &f : fragments) {

		//fragment location (in pixels):
		int32_t x = (int32_t)std::floor(f.fb_position.x);
		int32_t y = (int32_t)std::floor(f.fb_position.y);

		//if clipping is working properly, this condition shouldn't be needed;
		//however, it prevents crashes while you are working on your clipping functions,
		//so we suggest leaving it in place:
		if (x < 0 || (uint32_t)x >= framebuffer.width || y < 0 || (uint32_t)y >= framebuffer.height) {
			++out_of_range;
			continue;
		}

		//local names that refer to destination sample in framebuffer:
		float &fb_depth = framebuffer.depth_at(x,y,0);
		Spectrum &fb_color = framebuffer.color_at(x,y,0);

		//depth test:
		if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
			//"Always" means the depth test always passes.
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
			//"Never" means the depth test never passes.
			continue; //discard this fragment
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
			//"Less" means the depth test passes when the new fragment has depth less than the stored depth.
			//A1T4: Depth_Less
			// We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less"
			if (f.fb_position.z >= fb_depth)
				continue; // Depth test failed. Discard this fragment.
		} else {
			static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
		}

		//if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
		if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
			fb_depth = f.fb_position.z;
		}

		//shade fragment:
		ShadedFragment sf;
		sf.fb_position = f.fb_position;
		Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

		//write color to framebuffer if color writes aren't disabled:
		if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {

			//blend fragment:
			if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
				fb_color = sf.color;
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
				//A1T4: Blend_Add
				//framebuffer color should have fragment color multiplied by fragment opacity added to it.
				fb_color += sf.color * sf.opacity;
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
				//A1T4: Blend_Over
				// set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
				// You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
				fb_color = sf.color * sf.opacity + fb_color * (1.0f - sf.opacity);
			} else {
				static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
			}
		}
	}

	if (out_of_range > 0) {
		if constexpr (primitive_type == PrimitiveType::Lines) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely wrong with the clip_line function.", out_of_range);
		} else if constexpr (primitive_type == PrimitiveType::Triangles) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely wrong with the clip_triangle function.", out_of_range);
		}
	}



}

//-------------------------------------------------------------------------
//clipping functions

//helper to interpolate between vertices:
template< PrimitiveType p, class P, uint32_t F >
auto Pipeline< p, P, F >::lerp(ShadedVertex const &a, ShadedVertex const &b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  va, vb: endpoints of line
 *  emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 *
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::clip_line(
		ShadedVertex const &va, ShadedVertex const &vb,
		std::function< void(ShadedVertex const &) > const &emit_vertex
	) {
	//Determine portion of line over which:
	// pt = (b-a) * t + a
	// -pt.w <= pt.x <= pt.w
	// -pt.w <= pt.y <= pt.w
	// -pt.w <= pt.z <= pt.w

	//... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		//restrict range such that:
		//l + t * dl <= r + t * dr
		//re-arranging:
		// l - r <= t * (dr - dl)
		if (dr == dl) {
			//want: l - r <= 0
			if (l - r > 0.0f) {
				//works for none of range, so make range empty:
				min_t = 1.0f; max_t = 0.0f;
			}
		} else if (dr > dl) {
			//since dr - dl is positive:
			//want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { //dr < dl
			//since dr - dl is negative:
			//want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	//local names for clip positions and their difference:
	Vec4 const &a = va.clip_position;
	Vec4 const &b = vb.clip_position;
	Vec4 const ba = b-a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w,-ba.w, a.x, ba.x);
	clip_range( a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w,-ba.w, a.y, ba.y);
	clip_range( a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w,-ba.w, a.z, ba.z);
	clip_range( a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va,vb,min_t);
			//don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) out.attributes = va.attributes;
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va,vb,max_t);
			//don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) out.attributes = va.attributes;
			emit_vertex(out);
		}
	}
}


/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  va, vb, vc: vertices of triangle
 *  emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 *
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::clip_triangle(
		ShadedVertex const &va, ShadedVertex const &vb, ShadedVertex const &vc,
		std::function< void(ShadedVertex const &) > const &emit_vertex
	) {
	//A1EC: clip_triangle
	//TODO: correct code!
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

//-------------------------------------------------------------------------
//rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 *
 *        (x+0.5,y+1)
 *        /        \
 *    (x,y+0.5)  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points.
 *
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces
 *    a different rasterization result.
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 *
 */

template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::rasterize_line(
		ClippedVertex const &va, ClippedVertex const &vb,
		std::function< void(Fragment const &) > const &emit_fragment
	) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}

	// Bresenham's algorithm
	float x0 = std::round(va.fb_position.x) + 0.5f;
	float y0 = std::round(va.fb_position.y) + 0.5f;
	float x1 = std::round(vb.fb_position.x) + 0.5f;
	float y1 = std::round(vb.fb_position.y) + 0.5f;
    float dx = x0-x1;
    float dy = y0-y1;
	float z0 = va.fb_position.z;
	float z1 = vb.fb_position.z;
	float dz = z1 - z0;
	bool steep = std::abs(dy) > std::abs(dx);

	if (steep) {
		std::swap(x0, y0);
		std::swap(x1, y1);
		std::swap(dx, dy);
	}
	if (x0 > x1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
		std::swap(z0, z1);
	}

	// move through the line
	for (float x = x0; x < x1; x++) {
		float t = (x - x0) / dx;
		float y = y0 + dy * t;
		float z = z0 + dz * t;

		Fragment frag;
		if (steep) {
			frag.fb_position.x = y;
			frag.fb_position.y = x;
		} else {
			frag.fb_position.x = x;
			frag.fb_position.y = y;
		}
		frag.fb_position.z = z;
		frag.attributes = va.attributes;
		frag.derivatives.fill(Vec2(0.0f, 0.0f));

		// diamond-exit or if slope = (0 or infinity)
		if (std::abs(std::round(frag.fb_position.y) - frag.fb_position.y) < 0.5 || dx == 0 ||
			std::abs(std::round(frag.fb_position.x) - frag.fb_position.x) < 0.5 || dy == 0)
			emit_fragment(frag);

		// emit fragment only if it is on the top or left of a triangle
		else if (((std::abs(std::round(frag.fb_position.y) - frag.fb_position.y) <= 0.5) && (dy >= 0)) ||
				 ((std::abs(std::round(frag.fb_position.x) - frag.fb_position.x) <= 0.5) && (dx <= 0)))
			emit_fragment(frag);
	}
}


/*
 *
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  (x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Screen: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 *  The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 *
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::rasterize_triangle(
		ClippedVertex const &va, ClippedVertex const &vb, ClippedVertex const &vc,
		std::function< void(Fragment const &) > const &emit_fragment
	) {
	//NOTE: it is okay to restructure this function to allow these tasks to use the
	// same code paths. Be aware, however, that all of them need to remain working!
	// (e.g., if you break Flat while implementing Correct, you won't get points
	//  for Flat.)

	// sort vertices such that v1 < v2 < v3
	std::array<ClippedVertex, 3> vertices = {va, vb, vc};
	std::sort(vertices.begin(), vertices.end(), [](const ClippedVertex& a, const ClippedVertex& b) {
		return a.fb_position.y < b.fb_position.y || (a.fb_position.y == b.fb_position.y && a.fb_position.x < b.fb_position.x);
	});
	ClippedVertex& v1 = vertices[0];
	ClippedVertex& v2 = vertices[1];
	ClippedVertex& v3 = vertices[2];
	ClippedVertex a, b;

	// Calculate bounding box
	float minX = std::floor(std::min({va.fb_position.x, vb.fb_position.x, vc.fb_position.x}));
	float maxX = std::ceil(std::max({va.fb_position.x, vb.fb_position.x, vc.fb_position.x}));
	float minY = std::floor(std::min({va.fb_position.y, vb.fb_position.y, vc.fb_position.y}));
	float maxY = std::ceil(std::max({va.fb_position.y, vb.fb_position.y, vc.fb_position.y}));

	if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
		//A1T3: flat triangles

		float centerX, centerY, d1, d2, d3;
		float area, w1, w2, w3;
		bool has_neg, has_pos;
		// Edge walking
		for (float y = minY; y <= maxY; y++) {
			for (float x = minX; x <= maxX; x++) {
				centerX = x + 0.5f;
        		centerY = y + 0.5f;
				d1 = (centerX - v2.fb_position.x) * (v1.fb_position.y - v2.fb_position.y) - (v1.fb_position.x - v2.fb_position.x) * (centerY - v2.fb_position.y);
				d2 = (centerX - v3.fb_position.x) * (v2.fb_position.y - v3.fb_position.y) - (v2.fb_position.x - v3.fb_position.x) * (centerY - v3.fb_position.y);
				d3 = (centerX - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) - (v3.fb_position.x - v1.fb_position.x) * (centerY - v1.fb_position.y);

				has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
				has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

				if (!(has_neg && has_pos)) {
					// Interpolate Z
					area = (v1.fb_position.x - v3.fb_position.x) * (v2.fb_position.y - v3.fb_position.y) - (v2.fb_position.x - v3.fb_position.x) * (v1.fb_position.y - v3.fb_position.y);
					w1 = ((centerX - v3.fb_position.x) * (v2.fb_position.y - v3.fb_position.y) - (v2.fb_position.x - v3.fb_position.x) * (centerY - v3.fb_position.y)) / area;
					w2 = ((centerX - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) - (v3.fb_position.x - v1.fb_position.x) * (centerY - v1.fb_position.y)) / area;
					w3 = 1.0f - w1 - w2;

					Fragment frag;
					frag.fb_position = Vec3(centerX, centerY, w1 * v1.fb_position.z + w2 * v2.fb_position.z + w3 * v3.fb_position.z);
					frag.attributes = va.attributes; // Flat shading
					emit_fragment(frag);
				}
			}
		}
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Screen) {
		//A1T5: screen-space smooth triangles
		float lambda1, lambda2, lambda3, area;
		float lambda1_h, lambda2_h, lambda3_h;
		float area1, area2, attribute;
		float centerX, centerY, d1, d2, d3;
		float centerX_h, centerY_h;
		bool has_neg, has_pos;

		for (float y = minY; y <= maxY; y++) {
			for (float x = minX; x <= maxX; x++) {
				centerX = x + 0.5f;
        		centerY = y + 0.5f;

				d1 = (centerX - v2.fb_position.x) * (v1.fb_position.y - v2.fb_position.y) - (v1.fb_position.x - v2.fb_position.x) * (centerY - v2.fb_position.y);
				d2 = (centerX - v3.fb_position.x) * (v2.fb_position.y - v3.fb_position.y) - (v2.fb_position.x - v3.fb_position.x) * (centerY - v3.fb_position.y);
				d3 = (centerX - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) - (v3.fb_position.x - v1.fb_position.x) * (centerY - v1.fb_position.y);

				has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
				has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

				// Calculate barycentric coordinates
				area = (v2.fb_position.x - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
						(v2.fb_position.y - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);
				lambda1 = ((v2.fb_position.x - centerX) * (v3.fb_position.y - centerY) - (v2.fb_position.y - centerY) * (v3.fb_position.x - centerX)) / area;
				lambda2 = ((v3.fb_position.x - centerX) * (v1.fb_position.y - centerY) - (v3.fb_position.y - centerY) * (v1.fb_position.x - centerX)) / area;
				lambda3 = 1.0f - lambda1 - lambda2;

				if (!(has_neg && has_pos)) {
					float z = lambda1 * v1.fb_position.z + lambda2 * v2.fb_position.z + lambda3 * v3.fb_position.z;

					// Interpolate attributes
					std::array<float, FA> interpolatedAttributes;
					for (uint32_t i = 0; i < FA; i++) {
						interpolatedAttributes[i] = lambda1 * v1.attributes[i] + lambda2 * v2.attributes[i] + lambda3 * v3.attributes[i];
					}

					// Calculate derivatives for each attribute
					std::array<Vec2, FD> derivatives;
					float dAttribute_dx_forward, dAttribute_dy_forward;
					float dAttribute_dx_backward, dAttribute_dy_backward;
					float dAttribute_dx, dAttribute_dy;

					for (uint32_t i = 0; i < FD; i++) {
						// Forward differencing
						if (centerX < maxX) {
							centerX_h = centerX + 1.0f;
							centerY_h = centerY;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dx_forward = attribute - interpolatedAttributes[i];
						}
						if (centerY < maxY) {
							centerX_h = centerX;
							centerY_h = centerY + 1.0f;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dy_forward = attribute - interpolatedAttributes[i];
						}

						// Backward differencing
						if (centerX > minX) {
							centerX_h = centerX - 1.0f;
							centerY_h = centerY;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dy_forward = attribute - interpolatedAttributes[i];
							dAttribute_dx_backward = interpolatedAttributes[i] - attribute;
						}
						if (centerY > minY) {
							centerX_h = centerX;
							centerY_h = centerY - 1.0f;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dy_forward = attribute - interpolatedAttributes[i];
							dAttribute_dy_backward = interpolatedAttributes[i] - attribute;
						}

						// Choose the appropriate derivative based on the context
						if (std::fabs(dAttribute_dx_forward) < std::fabs(dAttribute_dx_backward) || centerX == maxX - 1) {
							dAttribute_dx = dAttribute_dx_forward;
						} else {
							dAttribute_dx = dAttribute_dx_backward;
						}

						if (std::fabs(dAttribute_dy_forward) < std::fabs(dAttribute_dy_backward) || centerY == maxY - 1) {
							dAttribute_dy = dAttribute_dy_forward;
						} else {
							dAttribute_dy = dAttribute_dy_backward;
						}
						derivatives[i] = Vec2(dAttribute_dx, dAttribute_dy);
					}
					Fragment frag;
					frag.fb_position = Vec3(centerX, centerY, z);
					frag.attributes = interpolatedAttributes;
					frag.derivatives = derivatives;
					emit_fragment(frag);
				}
			}
		}
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		//A1T5: perspective correct triangles
		float lambda1, lambda2, lambda3, area;
		float lambda1_h, lambda2_h, lambda3_h;
		float area1, area2, attribute;
		float centerX, centerY, d1, d2, d3;
		float centerX_h, centerY_h;
		bool has_neg, has_pos;

		for (float y = minY; y <= maxY; y++) {
			for (float x = minX; x <= maxX; x++) {
				centerX = x + 0.5f;
        		centerY = y + 0.5f;

				d1 = (centerX - v2.fb_position.x) * (v1.fb_position.y - v2.fb_position.y) - (v1.fb_position.x - v2.fb_position.x) * (centerY - v2.fb_position.y);
				d2 = (centerX - v3.fb_position.x) * (v2.fb_position.y - v3.fb_position.y) - (v2.fb_position.x - v3.fb_position.x) * (centerY - v3.fb_position.y);
				d3 = (centerX - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) - (v3.fb_position.x - v1.fb_position.x) * (centerY - v1.fb_position.y);

				has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
				has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

				// Calculate barycentric coordinates
				area = (v2.fb_position.x - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
						(v2.fb_position.y - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);
				lambda1 = ((v2.fb_position.x - centerX) * (v3.fb_position.y - centerY) - (v2.fb_position.y - centerY) * (v3.fb_position.x - centerX)) / area;
				lambda2 = ((v3.fb_position.x - centerX) * (v1.fb_position.y - centerY) - (v3.fb_position.y - centerY) * (v1.fb_position.x - centerX)) / area;
				lambda3 = 1.0f - lambda1 - lambda2;

				if (!(has_neg && has_pos)) {
					float z = lambda1 * v1.fb_position.z + lambda2 * v2.fb_position.z + lambda3 * v3.fb_position.z;

					// Interpolate inverse depth (w)
					float w = lambda1 * (1.0f / v1.inv_w) + lambda2 * (1.0f / v2.inv_w) + lambda3 * (1.0f / v3.inv_w);

					// Interpolate attributes with perspective correction
					std::array<float, FA> interpolatedAttributes;
					for (uint32_t i = 0; i < FA; i++) {
						interpolatedAttributes[i] = (lambda1 * v1.attributes[i] / v1.inv_w +
													lambda2 * v2.attributes[i] / v2.inv_w +
													lambda3 * v3.attributes[i] / v3.inv_w) / w;
					}

					// Calculate derivatives for each attribute
					std::array<Vec2, FD> derivatives;
					float dAttribute_dx_forward, dAttribute_dy_forward;
					float dAttribute_dx_backward, dAttribute_dy_backward;
					float dAttribute_dx, dAttribute_dy;

					for (uint32_t i = 0; i < FD; i++) {
						// Forward differencing
						if (centerX < maxX) {
							centerX_h = centerX + 1.0f;
							centerY_h = centerY;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dx_forward = attribute - interpolatedAttributes[i];
						}
						if (centerY < maxY) {
							centerX_h = centerX;
							centerY_h = centerY + 1.0f;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dy_forward = attribute - interpolatedAttributes[i];
						}

						// Backward differencing
						if (centerX > minX) {
							centerX_h = centerX - 1.0f;
							centerY_h = centerY;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dy_forward = attribute - interpolatedAttributes[i];
							dAttribute_dx_backward = interpolatedAttributes[i] - attribute;
						}
						if (centerY > minY) {
							centerX_h = centerX;
							centerY_h = centerY - 1.0f;

							// Calculate areas for the sub-triangles
							area1 = (v2.fb_position.x - centerX_h) * (v3.fb_position.y - centerY_h) -
									(v2.fb_position.y - centerY_h) * (v3.fb_position.x - centerX_h);
							area2 = (centerX_h - v1.fb_position.x) * (v3.fb_position.y - v1.fb_position.y) -
									(centerY_h - v1.fb_position.y) * (v3.fb_position.x - v1.fb_position.x);

							// Compute the barycentric coordinates
							lambda1_h = area1 / area;
							lambda2_h = area2 / area;
							lambda3_h = 1.0f - lambda1_h - lambda2_h;
							attribute = lambda1_h * v1.attributes[i] + lambda2_h * v2.attributes[i] + lambda3_h * v3.attributes[i];

							dAttribute_dy_forward = attribute - interpolatedAttributes[i];
							dAttribute_dy_backward = interpolatedAttributes[i] - attribute;
						}

						// Choose the appropriate derivative based on the context
						if (std::fabs(dAttribute_dx_forward) < std::fabs(dAttribute_dx_backward) || centerX == maxX - 1) {
							dAttribute_dx = dAttribute_dx_forward;
						} else {
							dAttribute_dx = dAttribute_dx_backward;
						}

						if (std::fabs(dAttribute_dy_forward) < std::fabs(dAttribute_dy_backward) || centerY == maxY - 1) {
							dAttribute_dy = dAttribute_dy_forward;
						} else {
							dAttribute_dy = dAttribute_dy_backward;
						}
						derivatives[i] = Vec2(dAttribute_dx, dAttribute_dy);
					}
					Fragment frag;
					frag.fb_position = Vec3(centerX, centerY, z);
					frag.attributes = interpolatedAttributes;
					frag.derivatives = derivatives;
					emit_fragment(frag);
				}
			}
		}
	}
}


//-------------------------------------------------------------------------
//compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Screen >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct >;
