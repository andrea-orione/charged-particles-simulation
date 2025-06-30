#ifndef VECTORS
#define VECTORS

/**
 * A 3D Vector.
 * Contains the x, y, z components as doubles.
 */
typedef struct {
  double x;
  double y;
  double z;
} Vec3;

/**
 * Sum 2 vectors.
 * @param v The first vector.
 * @param w The second vector.
 * @return The sum of the vectors.
 */
inline Vec3 sum(Vec3 v, Vec3 w) {
  return Vec3{v.x + w.x, v.y + w.y, v.z + w.z};
}

/**
 * Subtract 2 vectors.
 * @param v The first vector.
 * @param w The second vector.
 * @return The difference of the vectors.
 */
inline Vec3 subtract(Vec3 v, Vec3 w) {
  return Vec3{v.x - w.x, v.y - w.y, v.z - w.z};
}

/**
 * Compute the per scalar product of a vector.
 * @param v The vector.
 * @param a The scalar.
 * @return The product.
 */
inline Vec3 scale(Vec3 v, double a) { return Vec3{v.x * a, v.y * a, v.z * a}; }

/**
 * Take the dot product between two vectors.
 * @param v The first vector.
 * @param w The second vector.
 * @return The dot product.
 */
inline double dot(Vec3 v, Vec3 w) { return v.x * w.x + v.y * w.y + v.z * w.z; }

/**
 * Take the cross product between two vectors.
 * @param v The first vector.
 * @param w The second vector.
 * @return The cross product.
 */
inline Vec3 cross(Vec3 v, Vec3 w) {
  return Vec3{
      v.y * w.z - v.z * w.y,
      v.z * w.x - v.x * w.z,
      v.x * w.y - v.y * w.x,
  };
}

/**
 * A 3-dimentional vector field not depending on the position of the particles.
 * It should depend only on
 *  - the 3D position
 *  - the time
 */
typedef Vec3 (*Vec3Field)(Vec3, double);

#endif
