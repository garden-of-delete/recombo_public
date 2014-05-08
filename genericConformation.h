/* 
 * File:   genericConformation.h
 * Author: kmo
 *
 * Created on December 26, 2012, 12:51 PM
 */

#ifndef GENERICCONFORMATION_H
#define	GENERICCONFORMATION_H

#include <iostream>
#include <string>
#include <iterator>

#include "threevector.h"

#define TWOPI        6.28318530717958647693

/**
 * Generic algorithm to copy a range of coordinates into a range of 
 * three vectors.
 * 
 * @param start1 the beginning of range of scalars to be copied.
 * @param end1 the end of range of scalars to be copied.
 * @param start2 the beginning of range to three vectors into to which to copy.
 * 
 * @return the end of the range into which data was copied.
 */
template <class scalarIterator, class vectorIterator>
vectorIterator copyScalarToVector(scalarIterator start1, scalarIterator end1, vectorIterator start2)
{
   while(start1 != end1)
   {
      scalarIterator i(start1);
      start1++;
      scalarIterator j(start1);
      start1++;
      scalarIterator k(start1);
      start1++;
      start2->set(*i, *j, *k);
      start2++;
   }

   return start2;
}

/**
 * Generic algorithm to output a single threevector.
 * 
 * <p> For example, the vector to output is v = (1,2,3),
 * then calling outputVertex(v, os, "(", ")", ", ") will
 * output
 * 
 * <p> (1, 2, 3)
 * 
 * @param v a threevector
 * @param os the ostream to which output threevectors
 * @param prevector a string to output at the beginning of each three vector
 * @param postvector a string to output at the end of each three vector
 * @param betweenCoords a string to output between coordinates of v
 * 
 * @return a reference to the ostream
 */
template <class threevector>
std::ostream& outputVertex(const threevector& v,
std::ostream& os,
const std::string& prevector = "", const std::string& postvector = "",
const std::string& betweenCoords = " ")
{
   os << prevector << v.getX() << betweenCoords
           << v.getY() << betweenCoords
           << v.getZ() << postvector;
   return os;
}

/**
 * Generic algorithm to output a range of threevectors.
 * 
 * <p> For example, if the range consists of two vertices, (0,0,0) and (1,0,0),
 * then calling outputVertices(start, end, os, "[", "]", "; ", "(", ")", ", ") will
 * output
 * 
 * <p> [(0, 0, 0); (1, 0, 0)]
 * 
 * @param start1 the beginning of range to be output
 * @param end1 the end of range to be output
 * @param os the ostream to which to output threevectors
 * @param prelist a string to output before any three vectors
 * @param postlist a string to output after all other output
 * @param betweenVectors a string to output between threevectors
 * @param prevector a string to output at the beginning of each three vector
 * @param postvector a string to output at the end of each three vector
 * @param betweenCoords a string to output between coordinates threevectors
 * 
 * @return a reference to the ostream
 */
template <class vectorIterator>
std::ostream& outputVertices(vectorIterator start1, vectorIterator end1,
std::ostream& os,
const std::string& prelist = "", const std::string& postlist = "",
const std::string& betweenVectors = " ",
const std::string& prevector = "", const std::string& postvector = "",
const std::string& betweenCoords = " ")
{
   os << prelist;

   if(start1 != end1)
   {
      outputVertex(*start1, os, prevector, postvector, betweenCoords);
      start1++;
   }
   while(start1 != end1)
   {
      os << betweenVectors;
      outputVertex(*start1, os, prevector, postvector, betweenCoords);
      start1++;
   }

   os << postlist;

   return os;
}

/**
 * Generic algorithm to compute the radius of gyration squared of a range of threevectors.
 * 
 * @param start the beginning of range threevectors
 * @param end the end of range of threevectors
 * 
 * @return radius of gyration squared
 */
template <class vectorIterator> //, class T>
double computeRog(vectorIterator start, vectorIterator end)
{
   // first compute the center of mass
   vectorIterator i = start;
   int n = 0;
   typename std::iterator_traits<vectorIterator>::value_type::coordinate_type x, y, z;
   //   std::cout << std::endl << std::endl << typeid(x).name() << std::endl;
   x = y = z = 0;
   while(i != end)
   {
      x += i->getX();
      y += i->getY();
      z += i->getZ();
      i++;
      n++;
   }

   // next compute the second moment of displacement squared from center of mass
   double xbar, ybar, zbar, a, b, c, r2;
   xbar = double(x) / n;
   ybar = double(y) / n;
   zbar = double(z) / n;
   r2 = 0.0;
   i = start;
   while(i != end)
   {
      a = i->getX() - xbar;
      b = i->getY() - ybar;
      c = i->getZ() - zbar;
      r2 += a * a + b * b + c*c;
      i++;
   }
   return r2 / n;
}

//int count;

// TODO: fix computeWrAcn to work with list<>::iterator

/**
 * Generic algorithm to compute the contribution of Writhe and ACN between one
 * edge and a range of threevectors representing the other edges.
 * 
 * @param start the beginning of range threevectors
 * @param end the end of range of threevectors
 * @param Wr will be set to the computed writhe
 * @param Acn will be set to the computed ACN
 */
template <class vectorIterator>
void computeWrAcn(vectorIterator start, vectorIterator end,
vectorIterator vect1, vectorIterator vect2, double &space_writhe_sum,
double &gauss_xing_number_sum)
{
   if(start == end) return;
   //   int before = count;
   //   std::cout << "before = " << count << std::endl;
   typedef typename std::iterator_traits<vectorIterator>::value_type::coordinate_type T;
   // endpoints of two edges
   threevector<T> V1, V2, V3, V4;
   // displacements between the endpoints of two edges
   threevector<T> R13, R14, R24, R23, R34, R12;
   // vectors orthogonal to planes containing three points
   threevector<T> N1, N2, N3, N4;
   // norms of the N vectors
   double normN1, normN2, normN3, normN4;
   // crossproduct of two displacement vectors
   threevector<T> R34xR12;

   double exact_acn_contrib;

   V1 = *vect1;
   R12 = V2 = *vect2;
   R12 -= V1;

   vectorIterator j, k;
   j = start;
   k = j + 1;
   while(k != end)
   {
      //      std::cout << count++ << std::endl;
      //      count++;
      R13 = R23 = V3 = *j;
      R14 = R24 = R34 = V4 = *k;

      R13 -= V1;
      R14 -= V1;
      R24 -= V2;
      R23 -= V2;
      R34 -= V3;

      N1.cross(R13, R14);
      normN1 = N1.norm();
      N2.cross(R14, R24);
      normN2 = N2.norm();
      N3.cross(R24, R23);
      normN3 = N3.norm();
      N4.cross(R23, R13);
      normN4 = N4.norm();

      exact_acn_contrib = asin(double(N1.dot(N2)) / (normN1 * normN2))
              + asin(double(N2.dot(N3)) / (normN2 * normN3))
              + asin(double(N3.dot(N4)) / (normN3 * normN4))
              + asin(double(N4.dot(N1)) / (normN4 * normN1));

      if(!_isnan(exact_acn_contrib))
      {
         gauss_xing_number_sum += exact_acn_contrib;
         R34xR12.cross(R34, R12);


         if(R34xR12.dot(R13) > 0.0)
            space_writhe_sum += exact_acn_contrib;
         else
            space_writhe_sum -= exact_acn_contrib;
      }

      j = k;
      k++;
   }
   //   std::cout << "After = " << count << ", with difference = " << count - before << std::endl;

}

/**
 * Generic algorithm to compute the Writhe and ACN of a range of threevectors.
 * 
 * @param start the beginning of range threevectors.
 * @param end the end of range of threevectors.
 * @param last vertex in the range. It is expected that end == last + 1.
 * @param Wr will be set to the computed writhe.
 * @param Acn will be set to the computed ACN.
 */
template <class vectorIterator> //, class T>
void computeWrAcn(vectorIterator start, vectorIterator end,
vectorIterator last, double &Wr, double &Acn)
{
   //   count = 0;
   if(start == end || start + 1 == end || start + 2 == end)
   {
      // we are looking at an empty conformation, or a conformation 
      // consisting of only one edge, or a triangle, so there is almost nothing
      // to do
      Wr = Acn = 0;
      return;
   }

   double gauss_xing_number_sum = 0.0;
   double space_writhe_sum = 0.0;

   // we need some special treatment when one of the edges is the one between 
   // the first and the last vertex
   // compute the contribution from when one edge is the last edge
   computeWrAcn(start + 1, last, last, start, space_writhe_sum, gauss_xing_number_sum);

   vectorIterator i, j, k;
   i = start;
   j = i + 1;
   k = j + 1;
   while(k != last)
   {
      computeWrAcn(k, end, i, j, space_writhe_sum, gauss_xing_number_sum);
      i = j;
      j = k;
      k++;
   }

   Wr = space_writhe_sum / (TWOPI);
   Acn = gauss_xing_number_sum / (TWOPI);
   //   std::cout << "Counted " << count << " pairs of edges." << std::endl;
}

/**
 * A generic algorithm to find the last element in a range. This algorithm 
 * will not give sensible results for empty ranges.
 * 
 * @param start the beginning of a range.
 * @param end the end of a range.
 * @return the last element of a range. Will satisfy the condition that
 *  last + 1 == end.
 */
template <class iterator>
iterator findLast(iterator start, iterator end)
{
   iterator i = start;
   iterator j = i + 1;
   while(j != end)
   {
      i = j;
      j++;
   }
   return i;
}

/**
 * Generic algorithm to compute the Writhe and ACN of a range of threevectors.
 * For this algorithm, the position of the last vertex need not be known.
 * 
 * @param start the beginning of range threevectors.
 * @param end the end of range of threevectors.
 * @param Wr will be set to the computed writhe.
 * @param Acn will be set to the computed ACN.
 */
template <class vectorIterator> //, class T>
void computeWrAcn(vectorIterator start, vectorIterator end, double &Wr, double &Acn)
{
   computeWrAcn(start, end, findLast(start, end), Wr, Acn);
}

template <class vectorIterator>
void computeWrAcnTwocomps(vectorIterator start0, vectorIterator end0,
vectorIterator start1, vectorIterator end1, double &Wr, double &Acn)
{
   vectorIterator last0 = findLast(start0, end0);
   vectorIterator last1 = findLast(start1, end1);

   double gauss_xing_number_sum = 0.0;
   double space_writhe_sum = 0.0;

   // we need some special treatment when one of the edges is the one between 
   // the first and the last vertex
   // compute the contribution from when one edge is the last edge
   computeWrAcn(start0, end0, last1, start1, space_writhe_sum, gauss_xing_number_sum);
   computeWrAcn(start1, end1, last0, start0, space_writhe_sum, gauss_xing_number_sum);

   vectorIterator i, j, k;
   i = start0;
   j = i + 1;
   //   k = j+1;
   while(j != end0)
   {
      computeWrAcn(start1, end1, i, j, space_writhe_sum, gauss_xing_number_sum);
      i = j;
      //      j = k;
      j++;
   }

   // It remains only to compute the contributions from the two edges which
   // are between the first and last, but in different components
   typedef typename std::iterator_traits<vectorIterator>::value_type::coordinate_type T;

   // endpoints of two edges
   threevector<T> V1, V2, V3, V4;
   // displacements between the endpoints of two edges
   threevector<T> R13, R14, R24, R23, R34, R12;
   // vectors orthogonal to planes containing three points
   threevector<T> N1, N2, N3, N4;
   // norms of the N vectors
   double normN1, normN2, normN3, normN4;
   // crossproduct of two displacement vectors
   threevector<T> R34xR12;

   double exact_acn_contrib;

   V1 = *last0;
   R12 = V2 = *start0;
   R12 -= V1;

   //      std::cout << count++ << std::endl;
   //      count++;
   R13 = R23 = V3 = *last1;
   R14 = R24 = R34 = V4 = *start1;

   R13 -= V1;
   R14 -= V1;
   R24 -= V2;
   R23 -= V2;
   R34 -= V3;

   N1.cross(R13, R14);
   normN1 = N1.norm();
   N2.cross(R14, R24);
   normN2 = N2.norm();
   N3.cross(R24, R23);
   normN3 = N3.norm();
   N4.cross(R23, R13);
   normN4 = N4.norm();

   exact_acn_contrib = asin(double(N1.dot(N2)) / (normN1 * normN2))
           + asin(double(N2.dot(N3)) / (normN2 * normN3))
           + asin(double(N3.dot(N4)) / (normN3 * normN4))
           + asin(double(N4.dot(N1)) / (normN4 * normN1));

   if(!_isnan(exact_acn_contrib)) //changed for windows compatibility
   {
      gauss_xing_number_sum += exact_acn_contrib;
      R34xR12.cross(R34, R12);


      if(R34xR12.dot(R13) > 0.0)
         space_writhe_sum += exact_acn_contrib;
      else
         space_writhe_sum -= exact_acn_contrib;
   }



   Wr = space_writhe_sum / (TWOPI);
   Acn = gauss_xing_number_sum / (TWOPI);


}
#endif	/* GENERICCONFORMATION_H */

