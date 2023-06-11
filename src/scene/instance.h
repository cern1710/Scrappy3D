
#pragma once

#include <memory>
#include <variant>

#include "../geometry/halfedge.h"
#include "../lib/mathlib.h"

#include "camera.h"
#include "delta_light.h"
#include "env_light.h"
#include "material.h"
#include "particles.h"
#include "shape.h"
#include "skeleton.h"
#include "transform.h"

enum class DrawStyle : uint8_t {
	Wireframe, //lines at edges
	Flat,      //triangles with attributes from first vertex
	Smooth,    //triangles with attributes interpolated on-screen
	Correct,   //triangles with attributes interpolated in 3D
};

namespace Instance {

class Geometry_Settings {
public:
	bool visible = true;
	bool collides = true;
	DrawStyle style = DrawStyle::Correct;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		f("visible", t.visible);
		f("collides", t.collides);
		introspect_enum< I >(f, "style", t.style, std::vector< std::pair< const char *, DrawStyle > >{
			{"Wireframe", DrawStyle::Wireframe},
			{"Flat", DrawStyle::Flat},
			{"Smooth", DrawStyle::Smooth},
			{"Correct", DrawStyle::Correct}
		});
	}
	static inline const char *TYPE = "Instance::Geometry_Settings";
};

class Light_Settings {
public:
	bool visible = true;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		f("visible", t.visible);
	}
	static inline const char *TYPE = "Instance::Light_Settings";
};

class Simulate_Settings {
public:
	bool visible = true;
	bool wireframe = false;
	bool simulate_here = false;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		f("visible", t.visible);
		f("wireframe", t.wireframe);
		f("simulate_here", t.simulate_here);
	}
	static inline const char *TYPE = "Instance::Simulate_Settings";
};

class Mesh {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<Halfedge_Mesh> mesh;
	std::weak_ptr<Material> material;
	Geometry_Settings settings;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("mesh", t.mesh);
		if constexpr (I != Intent::Animate) f("material", t.material);
		Geometry_Settings::introspect< I >(std::forward< F >(f), t.settings);
	}
	static inline const char *TYPE = "Instance::Mesh";
};

class Skinned_Mesh {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<::Skinned_Mesh> mesh;
	std::weak_ptr<Material> material;
	Geometry_Settings settings;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("mesh", t.mesh);
		if constexpr (I != Intent::Animate) f("material", t.material);
		Geometry_Settings::introspect< I >(std::forward< F >(f), t.settings);
	}
	static inline const char *TYPE = "Instance::Skinned_Mesh";
};

class Shape {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<::Shape> shape;
	std::weak_ptr<Material> material;
	Geometry_Settings settings;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("shape", t.shape);
		if constexpr (I != Intent::Animate) f("material", t.material);
		Geometry_Settings::introspect< I >(std::forward< F >(f), t.settings);
	}
	static inline const char *TYPE = "Instance::Shape";
};

class Delta_Light {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<::Delta_Light> light;
	Light_Settings settings;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("light", t.light);
		Light_Settings::introspect< I >(std::forward< F >(f), t.settings);
	}
	static inline const char *TYPE = "Instance::Delta_Light";
};

class Environment_Light {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<::Environment_Light> light;
	Light_Settings settings;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("light", t.light);
		Light_Settings::introspect< I >(std::forward< F >(f), t.settings);
	}
	static inline const char *TYPE = "Instance::Environment_Light";
};

class Particles {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<Halfedge_Mesh> mesh;
	std::weak_ptr<Material> material;
	std::weak_ptr<::Particles> particles;
	Simulate_Settings settings;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("mesh", t.mesh);
		if constexpr (I != Intent::Animate) f("material", t.material);
		if constexpr (I != Intent::Animate) f("particles", t.particles);
		Simulate_Settings::introspect< I >(std::forward< F >(f), t.settings);
	}
	static inline const char *TYPE = "Instance::Particles";
};

class Camera {
public:
	std::weak_ptr<Transform> transform;
	std::weak_ptr<::Camera> camera;

	//- - - - - - - - - - - -
	template< Intent I, typename F, typename T >
	static void introspect(F&& f, T&& t) {
		if constexpr (I != Intent::Animate) f("transform", t.transform);
		if constexpr (I != Intent::Animate) f("camera", t.camera);
	}
	static inline const char *TYPE = "Instance::Camera";
};

} // namespace Instance
