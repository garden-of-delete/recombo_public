/*
 * File:   threevector.h
 * Author: kmo
 *
 * Created on November 27, 2012, 6:49 PM
 */

#ifndef THREEVECTOR_H
#define	THREEVECTOR_H

#include <cmath>
#include <iostream>

using namespace std;

template <class T>
class threevector;

template <class T>
std::ostream& operator <<(std::ostream&, const threevector<T>& v);

/**
 * A generic class implementing a vector in a three-dimensional space.
 */
template <class T>
class threevector
{
public:

   /**
    * Default constructor, which makes no guarantee about the values of coordinates.
    */
   threevector() { }

   /**
    * Constructor which explicitly sets coordinates.
    * @param x value to set first coordinate
    * @param y value to set second coordinate
    * @param z value to set third coordinate
    */
   inline threevector(T x, T y, T z);

   /**
    * Copy constructor.
    * @param orig new object will have same coordinates as orig
    */
   inline threevector(const threevector<T>& orig);

   /**
    * Constructor to copy coordinates from an array.
    * @param orig must be an array of at least length three. Constructor will
    * copy coordinates from first three entries of orig
    */
   inline threevector(const T* orig);

   /**
    * Constructor that will make a scaled copy of original threevector.
    * @param scale to multiply coordinates
    * @param orig coordinates of new threevector will be those of orig
    * multiplied by scale
    */
   inline threevector(T scale, const threevector<T>& orig);

   /**
    * Constructor that will make a scaled copy of coordinates from an array.
    * @param scale to multiply coordinates
    * @param orig expected to be an array of at least length three. The
    * coordinates of new threevector will be those from orig multiplied by
    * scale.
    */
   inline threevector(T scale, const T* orig);

   /**
    * Destructor.
    */
   ~threevector();

   /**
    * Accesses one coordinate.
    * @param i expected to be 0, 1 or 2.
    * @return A reference to coordinate specified by i.
    */
   inline T& operator [](int i);

   /**
    * Accesses one coordinate.
    * @param i expected to be 0, 1 or 2.
    * @return A constant reference to coordinate specified by i.
    */
   inline const T& operator [](int i) const;

   /**
    * Sets all three coordinates of vector.
    * @param x
    * @param y
    * @param z
    */
   inline void set(T x, T y, T z);

   /**
    * Accesses the first coordinate.
    * @return The value of the first coordinate.
    */
   inline T getX() const;

   /**
    * Accesses the second coordinate.
    * @return The value of the second coordinate.
    */
   inline T getY() const;

   /**
    * Accesses the third coordinate.
    * @return The value of the third coordinate.
    */
   inline T getZ() const;

   /**
    * Sets threevector to be a linear combination of two other threevectors.
    * @param a scale applied to u
    * @param u a threevector
    * @param b scale applied to v
    * @param v a threevector
    * @return *this set to a*u + b*v.
    */
   inline threevector<T>& combine(T a, const threevector<T>& u, T b, const threevector<T>& v);

   /**
    * Computes the inner product of threevector with another three vector
    * @param v another threevector
    * @return The inner product of *this and v.
    */
   inline T dot(const threevector<T>& v) const;

   /**
    * Sets threevector to be the cross product of two other threevectors
    * @param u a threevector
    * @param v a second threevector
    * @return *this set to the cross product of u and v.
    */
   inline threevector<T>& cross(const threevector<T>& u, const threevector<T>& v);

   /**
    * Computes the Euclidean norm of threevector.
    * @return The Euclidean norm.
    */
   inline double norm() const;

   /**
    * Computes the square of the Euclidean norm of threevector.
    * @return The Euclidean norm squared.
    */
   inline T norm2() const;

   /**
    * Tests for equality of *this and another threevector.
    * @param v another threevector.
    * @return true if v has same three coordinates as this, false otherwise.
    */
   inline bool operator ==(const threevector<T>& v) const;

   /**
    * Tests for equality of *this and another threevector are within eplsion of
    * each other in the Euclidean metric.
    * @param v another threevector.
    * @param epsilon expected to be nonnegative.
    * @return true if the vector difference of this and v has same three
    * has norm less than epsilon, false otherwise.
    */
   inline bool equalsEpsilon(const threevector<T>& v, double epsilon) const;

   /**
    * Opposite of ==.
    * @param v another threevector.
    * @return false if v has same three coordinates as this, true otherwise.
    */
   inline bool operator !=(const threevector<T>& v) const;

   /**
    * Assignment operator.
    * @param v another threevector.
    * @return *this with coordinates set to those of v.
    */
   inline threevector<T>& operator =(const threevector<T>& v);

   /**
    * Produces the sum with a second threevector
    * @param v another threevector
    * @return *this after adding the coordinates of v to the coordinates of
    * *this.
    */
   inline threevector<T>& operator +=(const threevector<T>& v);

   /**
    * Produces the difference with a second threevector
    * @param v another threevector
    * @return *this after subtracting the coordinates of v to the coordinates of
    * *this.
    */
   inline threevector<T>& operator -=(const threevector<T>& v);

   /**
    * Scales *this, that is multiplies the coordinates by a scale factor.
    * @param scale a scale factor.
    * @return *this after scaling.
    */
   inline threevector<T>& operator *=(T scale);

   /**
    * Output operator.
    */
   friend std::ostream& operator << <>(std::ostream&, const threevector<T>& v);

   /**
    * The type of the coordinates of threevector.
    */
   typedef T coordinate_type;

protected:
   T coords[3];
};

template <class T>
std::ostream& operator <<(std::ostream& out, const threevector<T>& v)
{
   //return out << v.getX() << " " << v.getY() << " " << v.getZ(); preserved original
   return out << v.getX() << " " << v.getY() << " " << v.getZ() << endl;

}

template <class T>
T& threevector<T>::operator [](int i)
{
   return coords[i];
}

template <class T>
const T& threevector<T>::operator [](int i) const
{
   return coords[i];
}

template <class T>
threevector<T>::threevector(T x, T y, T z)
{
   coords[0] = x;
   coords[1] = y;
   coords[2] = z;
}

template <class T>
threevector<T>::threevector(const threevector<T>& orig)
{
   coords[0] = orig.coords[0];
   coords[1] = orig.coords[1];
   coords[2] = orig.coords[2];
}

template <class T>
threevector<T>::threevector(const T* orig)
{
   coords[0] = orig[0];
   coords[1] = orig[1];
   coords[2] = orig[2];
}

template <class T>
threevector<T>::threevector(T scale, const threevector<T>& orig)
{
   coords[0] = scale * orig.coords[0];
   coords[1] = scale * orig.coords[1];
   coords[2] = scale * orig.coords[2];
}

template <class T>
threevector<T>::threevector(T scale, const T* orig)
{
   coords[0] = scale * orig[0];
   coords[1] = scale * orig[1];
   coords[2] = scale * orig[2];
}

template <class T>
threevector<T>::~threevector() { }

template <class T>
T threevector<T>::getX() const
{
   return coords[0];
}

template <class T>
T threevector<T>::getY() const
{
   return coords[1];
}

template <class T>
T threevector<T>::getZ() const
{
   return coords[2];
}

template <class T>
bool threevector<T>::operator ==(const threevector<T>& v) const
{
   return coords[0] == v.coords[0]
           && coords[1] == v.coords[1]
           && coords[2] == v.coords[2];
}

template <class T>
bool threevector<T>::equalsEpsilon(const threevector<T>& v, double epsilon) const
{
   double a = coords[0] - v.coords[0];
   double b = coords[1] - v.coords[1];
   double c = coords[2] - v.coords[2];
   return a * a + b * b + c * c < epsilon * epsilon;
}

template <class T>
bool threevector<T>::operator !=(const threevector<T>& v) const
{
   return !(*this == v);
}

template <class T>
threevector<T>& threevector<T>::operator *=(T scale)
{
   coords[0] *= scale;
   coords[1] *= scale;
   coords[2] *= scale;
   return *this;
}

template <class T>
threevector<T>& threevector<T>::operator =(const threevector<T>& v)
{
   coords[0] = v.coords[0];
   coords[1] = v.coords[1];
   coords[2] = v.coords[2];
   return *this;
}

template <class T>
void threevector<T>::set(T x, T y, T z)
{
   coords[0] = x;
   coords[1] = y;
   coords[2] = z;
}

template <class T>
threevector<T>& threevector<T>::operator +=(const threevector<T>& v)
{
   coords[0] += v.coords[0];
   coords[1] += v.coords[1];
   coords[2] += v.coords[2];
   return *this;
}

template <class T>
threevector<T>& threevector<T>::operator -=(const threevector<T>& v)
{
   coords[0] -= v.coords[0];
   coords[1] -= v.coords[1];
   coords[2] -= v.coords[2];
   return *this;
}

template <class T>
threevector<T>& threevector<T>::combine(T a, const threevector<T>& u, T b, const threevector<T>& v)
{
   coords[0] = a * u.coords[0] + b * v.coords[0];
   coords[1] = a * u.coords[1] + b * v.coords[1];
   coords[2] = a * u.coords[2] + b * v.coords[2];
   return *this;
}

template <class T>
T threevector<T>::dot(const threevector<T>& v) const
{
   return coords[0] * v.coords[0] +
           coords[1] * v.coords[1] +
           coords[2] * v.coords[2];
}

template <class T>
double threevector<T>::norm() const
{
   return std::sqrt(norm2());
}

template <class T>
T threevector<T>::norm2() const
{
   return coords[0] * coords[0] +
           coords[1] * coords[1] +
           coords[2] * coords[2];
}

template <class T>
threevector<T>& threevector<T>::cross(const threevector<T>& u, const threevector<T>& v)
{
   coords[0] = u.coords[1] * v.coords[2] - u.coords[2] * v.coords[1];
   coords[1] = u.coords[2] * v.coords[0] - u.coords[0] * v.coords[2];
   coords[2] = u.coords[0] * v.coords[1] - u.coords[1] * v.coords[0];
   return *this;
}

#endif	/* THREEVECTOR_H */

