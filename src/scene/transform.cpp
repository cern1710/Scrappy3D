
#include "transform.h"

Mat4 Transform::local_to_parent() const {
	return Mat4::translate(translation) * rotation.to_mat() * Mat4::scale(scale);
}

Mat4 Transform::parent_to_local() const {
	return Mat4::scale(1.0f / scale) * rotation.inverse().to_mat() * Mat4::translate(-translation);
}

Mat4 Transform::local_to_world() const {
	return parent.lock() ? parent.lock()->local_to_world() * local_to_parent() : local_to_parent();
}

Mat4 Transform::world_to_local() const {
	return parent.lock() ? parent_to_local() * parent.lock()->world_to_local() : parent_to_local();
}

bool operator!=(const Transform& a, const Transform& b) {
	return a.parent.lock() != b.parent.lock() || a.translation != b.translation ||
	       a.rotation != b.rotation || a.scale != b.scale;
}
